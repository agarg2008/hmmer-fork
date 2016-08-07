/* nhmmer: search profile HMM(s) against a nucleotide sequence database.
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_scorematrix.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"


#ifdef __cplusplus
#    define _FUNCTION_MACRO_ __PRETTY_FUNCTION__
#    include <cstdlib>
#    include <cstdio>
#    include <cstdarg>
#else
#    define _FUNCTION_MACRO_ __func__
#    include <stdlib.h>
#    include <stdio.h>
#    include <stdarg.h>
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#if !NDEBUG
#    define LDB(...) log_debug(_FUNCTION_MACRO_, __FILE__, __LINE__, ##__VA_ARGS__);
#else
#    define LDB(...)
#endif

static inline void log_debug(const char *func, const char *filename, int line, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "[D:%s:%d:%s] ", func, line, filename);
    vfprintf(stderr, fmt, args);
    va_end(args);
}

#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif /*HMMER_THREADS*/

#include "hmmer.h"

/* set the max residue count to 1/4 meg when reading a block */
#ifdef P7_IMPL_DUMMY_INCLUDED
#include "esl_vectorops.h"
#define NHMMER_MAX_RESIDUE_COUNT (1024 * 100)
#else
#define NHMMER_MAX_RESIDUE_COUNT (1024 * 256)  /* 1/4 Mb */
#endif

typedef struct {
  P7_HMM      **h;
  size_t        n;
  size_t        m;
} HMM_VEC;

HMM_VEC p7_hmmvec_Fill(P7_HMMFILE *hfp, ESL_ALPHABET *abc) {
  int status = eslOK;
  HMM_VEC hv = (HMM_VEC) {NULL, 0uL, 4uL};
  ESL_ALLOC(hv.h, 4uL * sizeof(P7_HMM *));
  if(abc == NULL) p7_Fail("ESL_ALPHABET abc must be set to call %s.\n", __func__);
  while((status = p7_hmmfile_Read(hfp, &abc, hv.h + hv.n)) == eslOK) {
    if(++hv.n > hv.m) {
      hv.m = hv.n; kroundup32(hv.m);
      ESL_REALLOC(hv.h, hv.m * sizeof(P7_HMM *));
    }
  }
  LDB("Successfully loaded %lu hmm profiles.\n", hv.n);
  return hv;
  ERROR:
    p7_Fail("Could not allocate memory.\n");
  return hv;
}

void p7_hmmvec_Destroy(HMM_VEC *hv) {
  unsigned i = hv->n;
  while(i) p7_hmm_Destroy(hv->h[--i]);
  free(hv->h);
}

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  P7_BG            *bg;          /* null model                              */
  P7_PIPELINE      *pli;         /* work pipeline                           */
  P7_TOPHITS       *th;          /* top hit results                         */
  P7_OPROFILE      *om;          /* optimized query profile                 */
  P7_SCOREDATA     *scoredata;   /* hmm-specific data used by nhmmer */
} WORKER_INFO;


typedef struct {
  int    id;         /* internal sequence ID  */
  int    length;     /* length of sequence */
} ID_LENGTH;

typedef struct {
  ID_LENGTH  *id_lengths;
  int        count;
  int        size;
} ID_LENGTH_LIST;


static ID_LENGTH_LIST* init_id_length( int size );
static void            destroy_id_length( ID_LENGTH_LIST *list );
static int             add_id_length(ID_LENGTH_LIST *list, int id, int L);
static int             assign_Lengths(P7_TOPHITS *th, ID_LENGTH_LIST *id_length_list);

#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

#define CPUOPTS     NULL
#define MPIOPTS     NULL


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp              help                                                      docgroup*/
  { "-h",           eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  NULL,             "show brief help on version and usage",                         1 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE,      NULL, NULL, NULL,    NULL,  NULL,  NULL,              "direct output to file <f>, not stdout",                      2 },
  { "--stats",           eslARG_OUTFILE,      NULL, NULL, NULL,    NULL,  NULL,  NULL,              "direct stats output to file <f>, not stderr",           2 },
  { "--noali",      eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                2 },
  { "--acc",      eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                2 },
  { "--noheader",      eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "Don't output header",                                     2 },
  /* Control of scoring system */
  { "--singlemx",   eslARG_NONE,        FALSE,   NULL, NULL,    NULL,  NULL,   "",           "use substitution score matrix w/ single-sequence MSA-format inputs",  3 },
  { "--popen",      eslARG_REAL,       "0.03125",NULL,"0<=x<0.5",NULL, NULL, NULL,           "gap open probability",                                         3 },
  { "--pextend",    eslARG_REAL,       "0.75", NULL,  "0<=x<1",  NULL, NULL, NULL,           "gap extend probability",                                       3 },
  { "--mx",         eslARG_STRING,     "DNA1", NULL, NULL,      NULL,  NULL, "--mxfile",     "substitution score matrix choice (of some built-in matrices)", 3 },
  { "--mxfile",     eslARG_INFILE,       NULL, NULL, NULL,      NULL,  NULL,  "--mx",        "read substitution score matrix from file <f>",                 3 },
  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,       "10.0", NULL, "x>0",   NULL,  NULL,  REPOPTS,         "report sequences <= this E-value threshold in output",         4 },
  { "-T",           eslARG_REAL,        FALSE, NULL, NULL,    NULL,  NULL,  REPOPTS,         "report sequences >= this score threshold in output",           4 },
  /* Control of inclusion (significance) thresholds */
  { "--incE",       eslARG_REAL,       "0.01", NULL, "x>0",   NULL,  NULL,  INCOPTS,         "consider sequences <= this E-value threshold as significant",  5 },
  { "--incT",       eslARG_REAL,        FALSE, NULL, NULL,    NULL,  NULL,  INCOPTS,         "consider sequences >= this score threshold as significant",    5 },
  /* Model-specific thresholding for both reporting and inclusion */
  { "--cut_ga",     eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's GA gathering cutoffs to set all thresholding",   6 },
  { "--cut_nc",     eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's NC noise cutoffs to set all thresholding",       6 },
  { "--cut_tc",     eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's TC trusted cutoffs to set all thresholding",     6 },
  /* Control of acceleration pipeline */
  { "--max",        eslARG_NONE,        FALSE,      NULL, NULL,    NULL,  NULL, "--F1,--F2,--F3", "Turn all heuristic filters off (less speed, more power)",      7 },
  { "--F1",         eslARG_REAL, /*set below*/NULL, NULL, NULL,    NULL,  NULL, "--max",          "Stage 1 (SSV) threshold: promote hits w/ P <= F1",             7 },
  { "--F2",         eslARG_REAL,       "3e-3",      NULL, NULL,    NULL,  NULL, "--max",          "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             7 },
  { "--F3",         eslARG_REAL,       "3e-5",      NULL, NULL,    NULL,  NULL, "--max",          "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             7 },
  { "--nobias",     eslARG_NONE,         NULL,      NULL, NULL,    NULL,  NULL, "--max",          "turn off composition bias filter",                             7 },

  /* Selecting the alphabet rather than autoguessing it */
  { "--dna",        eslARG_NONE,        FALSE, NULL, NULL,   "--rna", NULL,   NULL,          "input alignment is DNA sequence data",                         8 },
  { "--rna",        eslARG_NONE,        FALSE, NULL, NULL,   "--dna",  NULL,  NULL,          "input alignment is RNA sequence data",                         8 },

/* Other options */
  { "--nonull2",    eslARG_NONE,         NULL, NULL, NULL,    NULL,  NULL,           NULL,     "turn off biased composition score corrections",                 12 },
  { "-Z",           eslARG_REAL,        FALSE, NULL, "x>0",   NULL,  NULL,           NULL,     "set database size (Megabases) to <x> for E-value calculations", 12 },
  { "--seed",       eslARG_INT,          "42", NULL, "n>=0",  NULL,  NULL,           NULL,     "set RNG seed to <n> (if 0: one-time arbitrary seed)",           12 },
  { "--w_beta",     eslARG_REAL,         NULL, NULL, NULL,    NULL,  NULL,           NULL,     "tail mass at which window length is determined",                12 },
  { "--w_length",   eslARG_INT,          NULL, NULL, NULL,    NULL,  NULL,           NULL,     "window length - essentially max expected hit length" ,          12 },
  { "--block_length", eslARG_INT,        NULL, NULL, "n>=50000", NULL, NULL,         NULL,     "length of blocks read from target database (threaded) ",        12 },
  { "--watson",     eslARG_NONE,         NULL, NULL, NULL,    NULL,  NULL,       "--crick",  "only search the top strand",                                    12 },
  { "--crick",      eslARG_NONE,         NULL, NULL, NULL,    NULL,  NULL,       "--watson",  "only search the bottom strand",                                 12 },


  /* stage-specific window length used for bias composition estimate,
   * hidden because they are confusing/expert options. May drag them out
   * into the daylight eventually
   */
  { "--B1",         eslARG_INT,         "110", NULL, NULL,    NULL,  NULL, "--max,--nobias", "window length for biased-composition modifier (SSV)",          99 },
  { "--B2",         eslARG_INT,         "240", NULL, NULL,    NULL,  NULL, "--max,--nobias", "window length for biased-composition modifier (Vit)",          99 },
  { "--B3",         eslARG_INT,        "1000", NULL, NULL,    NULL,  NULL, "--max,--nobias", "window length for biased-composition modifier (Fwd)",          99 },


/* Not used, but retained because esl option-handling code errors if it isn't kept here.  Placed in group 99 so it doesn't print to help*/
  { "--domZ",       eslARG_REAL,        FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "Not used",   99 },
  { "--domE",       eslARG_REAL,       "10.0", NULL, "x>0",   NULL,  NULL,  DOMREPOPTS,      "Not used",   99 },
  { "--domT",       eslARG_REAL,        FALSE, NULL, NULL,    NULL,  NULL,  DOMREPOPTS,      "Not used",   99 },
  { "--incdomE",    eslARG_REAL,       "0.01", NULL, "x>0",   NULL,  NULL,  INCDOMOPTS,      "Not used",   99 },
  { "--incdomT",    eslARG_REAL,        FALSE, NULL, NULL,    NULL,  NULL,  INCDOMOPTS,      "Not used",   99 },



#ifdef HMMER_THREADS
  { "--cpu",        eslARG_INT, NULL,"HMMER_NCPU","n>=0",NULL,  NULL,  CPUOPTS,         "number of parallel CPU workers to use for multithreads",      12 },
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};


/* struct cfg_s : "Global" application configuration shared by all threads/processes
 *
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  char            *dbfile;            /* target sequence database file                   */
  char            *queryfile;           /* query HMM file                                  */
  int              qfmt;

  int              do_mpi;            /* TRUE if we're doing MPI parallelization         */
  int              nproc;             /* how many MPI processes, total                   */
  int              my_rank;           /* who am I, in 0..nproc-1                         */

  char             *firstseq_key;     /* name of the first sequence in the restricted db range */
  int              n_targetseq;       /* number of sequences in the restricted range */
};

static char usage[]  = "[options] <query hmmfile|alignfile> <target seqfile>";
static char banner[] = "search a DNA model or alignment against a DNA database";
//static char usage[]  = "[options] <hmmfile> <seqdb>";
//static char banner[] = "search a DNA model against a DNA database";


static int  serial_master  (ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop    (WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_SQFILE *dbfp, char *firstseq_key, int n_targetseqs /*, ESL_STOPWATCH *ssv_watch_master, ESL_STOPWATCH *postssv_watch_master, ESL_STOPWATCH *watch_slave*/);
#ifdef HMMER_THREADS
#define BLOCK_SIZE 1000

static int  thread_loop(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp, char *firstseq_key, int n_targetseqs);
static void pipeline_thread(void *arg);

#endif /*HMMER_THREADS*/


static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_queryfile, char **ret_seqfile)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  int          status;

  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { if (printf("Failed to process environment: %s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { if (printf("Failed to parse command line: %s\n", go->errbuf)  < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { if (printf("Failed to parse command line: %s\n", go->errbuf)  < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE)
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      if (puts("\nBasic options:")                                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 100); /* 1= group; 2 = indentation; 120=textwidth*/

      if (puts("\nOptions directing output:")                                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 100);

      if (puts("\nOptions controlling scoring system:")                      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 100);

      if (puts("\nOptions controlling reporting thresholds:")                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 100);

      if (puts("\nOptions controlling inclusion (significance) thresholds:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 100);

      if (puts("\nOptions controlling model-specific thresholding:")         < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);

      if (puts("\nOptions controlling acceleration heuristics:")             < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 100);

      if (puts("\nOptions for selecting query alphabet rather than guessing it:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80);

//      if (puts("\nOptions for restricting search to a range of target database sequences:")             < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
//      esl_opt_DisplayHelp(stdout, go, 8, 2, 100);

      if (puts("\nOptions controlling seed search heuristic:")               < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 100);

      if (puts("\nOther expert options:")                                    < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 100);
      exit(0);

  }

  if (esl_opt_ArgNumber(go)                  != 2)     { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_queryfile = esl_opt_GetArg(go, 1)) == NULL)  { if (puts("Failed to get <queryfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_seqfile   = esl_opt_GetArg(go, 2)) == NULL)  { if (puts("Failed to get <seqdb> argument on command line")   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* Validate any attempted use of stdin streams */
  if (strcmp(*ret_queryfile, "-") == 0 && strcmp(*ret_seqfile, "-") == 0)
    { if (puts("Either <query hmmfile|alignfile> or <seqdb> may be '-' (to read from stdin), but not both.") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  *ret_go = go;
  return eslOK;

 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  if (puts("\nwhere basic options are:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  if (printf("\nTo see more help on available options, do %s -h\n\n", argv[0]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_getopts_Destroy(go);
  exit(1);

 ERROR:
  if (go) esl_getopts_Destroy(go);
  exit(status);
}

static int
output_header(FILE *ofp, const ESL_GETOPTS *go, char *queryfile, char *seqfile, int ncpus)
{
  if(esl_opt_IsUsed(go, "--noheader")) return eslOK;
  p7_banner(ofp, go->argv[0], banner);
  if (fprintf(ofp, "# query file:                      %s\n", queryfile)                                                                                < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(ofp, "# target sequence database:        %s\n", seqfile)                                                                                < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-o")              && fprintf(ofp, "# output directed to file:         %s\n",            esl_opt_GetString(go, "-o"))             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "--singlemx")   && fprintf(ofp, "# Use score matrix for 1-seq MSAs:  on\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--popen")      && fprintf(ofp, "# gap open probability:            %f\n",             esl_opt_GetReal   (go, "--popen"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--pextend")    && fprintf(ofp, "# gap extend probability:          %f\n",             esl_opt_GetReal   (go, "--pextend"))  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--mx")         && fprintf(ofp, "# subst score matrix (built-in):   %s\n",             esl_opt_GetString (go, "--mx"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--mxfile")     && fprintf(ofp, "# subst score matrix (file):       %s\n",             esl_opt_GetString (go, "--mxfile"))   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-E")           && fprintf(ofp, "# sequence reporting threshold:    E-value <= %g\n",  esl_opt_GetReal   (go, "-E"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-T")           && fprintf(ofp, "# sequence reporting threshold:    score >= %g\n",    esl_opt_GetReal   (go, "-T"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incE")       && fprintf(ofp, "# sequence inclusion threshold:    E-value <= %g\n",  esl_opt_GetReal   (go, "--incE"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incT")       && fprintf(ofp, "# sequence inclusion threshold:    score >= %g\n",    esl_opt_GetReal   (go, "--incT"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_ga")     && fprintf(ofp, "# model-specific thresholding:     GA cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_nc")     && fprintf(ofp, "# model-specific thresholding:     NC cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_tc")     && fprintf(ofp, "# model-specific thresholding:     TC cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--max")        && fprintf(ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n")                       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F1")         && fprintf(ofp, "# SSV filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F1"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F2")         && fprintf(ofp, "# Vit filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F2"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F3")         && fprintf(ofp, "# Fwd filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F3"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--nobias")     && fprintf(ofp, "# biased composition HMM filter:   off\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "--B1")         && fprintf(ofp, "# biased comp SSV window len:      %d\n",             esl_opt_GetInteger(go, "--B1"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--B2")         && fprintf(ofp, "# biased comp Viterbi window len:  %d\n",             esl_opt_GetInteger(go, "--B2"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--B3")         && fprintf(ofp, "# biased comp Forward window len:  %d\n",             esl_opt_GetInteger(go, "--B3"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "--dna")        && fprintf(ofp, "# input query is asserted as:      DNA\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--rna")        && fprintf(ofp, "# input query is asserted as:      RNA\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");


  if (esl_opt_IsUsed(go, "--nonull2")    && fprintf(ofp, "# null2 bias corrections:          off\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "--watson")    && fprintf(ofp, "# search only top strand:          on\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--crick") && fprintf(ofp, "# search only bottom strand:       on\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-Z")           && fprintf(ofp, "# database size is set to:         %.1f Mb\n",        esl_opt_GetReal(go, "-Z"))            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed") == 0 && fprintf(ofp, "# random number seed:              one-time arbitrary\n")                              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    else if                              (  fprintf(ofp, "# random number seed set to:       %d\n",             esl_opt_GetInteger(go, "--seed"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }
  if (esl_opt_IsUsed(go, "--w_beta")     && fprintf(ofp, "# window length beta value:        %g\n",             esl_opt_GetReal(go, "--w_beta"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--w_length")   && fprintf(ofp, "# window length :                  %d\n",             esl_opt_GetInteger(go, "--w_length")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--block_length")&&fprintf(ofp, "# block length :                   %d\n",             esl_opt_GetInteger(go, "--block_length")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
#ifdef HMMER_THREADS
  //if (esl_opt_IsUsed(go, "--cpu")        && fprintf(ofp, "# number of worker threads:        %d\n",             esl_opt_GetInteger(go, "--cpu"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(ofp, "# number of worker threads:        %d\n",             ncpus)      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
#endif
  if (fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}



int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go       = NULL;
  struct cfg_s     cfg;
  int              status   = eslOK;

  impl_Init();             /* processor specific initialization */
  p7_FLogsumInit();        /* we're going to use table-driven Logsum() approximations at times */

  /* Initialize what we can in the config structure (without knowing the alphabet yet)
   */
  cfg.queryfile  = NULL;
  cfg.dbfile     = NULL;
  cfg.qfmt       = eslMSAFILE_UNKNOWN;
  cfg.do_mpi     = FALSE;               /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;                   /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;                   /* this gets reset below, if we init MPI */

  cfg.firstseq_key = NULL;
  cfg.n_targetseq  = -1;

  process_commandline(argc, argv, &go, &cfg.queryfile, &cfg.dbfile);

  status = serial_master(go, &cfg);

  esl_getopts_Destroy(go);
  return status;
}



/* serial_master()
 * The serial version of hmmsearch.
 * For each query HMM in <queryfile> search the database for hits.
 *
 * A master can only return if it's successful. All errors are handled
 * immediately and fatally with p7_Fail().
 */
static int
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp          = stdout;          /* results output file (-o)                              */
  FILE            *statsfp      = stderr;

  /*Some fraction of these will be used, depending on what sort of input is used for the query*/
  P7_HMMFILE      *hfp        = NULL;              /* open input HMM file    */
  P7_HMM          *hmm        = NULL;              /* one HMM query          */
  ESL_SQ          *qsq        = NULL;              /* query sequence         */

  int              dbformat  =  eslSQFILE_FASTA; /* format of dbfile          */
  ESL_SQFILE      *dbfp      = NULL;               /* open input sequence file  */

  ESL_ALPHABET    *abc       = NULL;              /* digital alphabet           */
  ESL_STOPWATCH   *w;
  P7_SCOREDATA    *scoredata = NULL;

  int              nquery    = 0;
  int              status    = eslOK;
  int              qhstatus  = eslOK;
  int              sstatus   = eslOK;
  int              i;
  double           resCnt    = 0;
  /* used to keep track of the lengths of the sequences that are processed */
  ID_LENGTH_LIST  *id_length_list = NULL;


  int              ncpus    = 0;

  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;
#ifdef HMMER_THREADS
  ESL_SQ_BLOCK    *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif
  char   errbuf[eslERRBUFSIZE];
  double window_beta = -1.0 ;
  int window_length  = -1;

  P7_BUILDER       *builder     = NULL;


//  ESL_STOPWATCH *ssv_watch_master = esl_stopwatch_Create();
//  ESL_STOPWATCH *ssv_watch_master_tot = esl_stopwatch_Create();

//  ESL_STOPWATCH *postssv_watch_master = esl_stopwatch_Create();
//  ESL_STOPWATCH *postssv_watch_master_tot = esl_stopwatch_Create();

//  ESL_STOPWATCH *watch_slave  = esl_stopwatch_Create();

  if (esl_opt_IsUsed(go, "--w_beta")) { if (  ( window_beta   = esl_opt_GetReal(go, "--w_beta") )  < 0 || window_beta > 1  ) esl_fatal("Invalid window-length beta value\n"); }
  if (esl_opt_IsUsed(go, "--w_length")) { if (( window_length = esl_opt_GetInteger(go, "--w_length")) < 4  ) esl_fatal("Invalid window length value\n"); }

  w = esl_stopwatch_Create();


  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbformat, p7_SEQDBENV, &dbfp);
  switch(status) {
       default:  p7_Fail("Unexpected error %d opening target sequence database file %s\n", status, cfg->dbfile);
       case eslENOTFOUND: p7_Fail("Failed to open target sequence database %s for reading\n",      cfg->dbfile);
       case eslEFORMAT:   p7_Fail("Target sequence database file %s is empty or misformatted\n",   cfg->dbfile);
       case eslEINVAL:    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
       case eslOK: break;
  }

  if ( esl_opt_IsOn(go, "--dna") )
    abc     = esl_alphabet_Create(eslDNA);
  else if ( esl_opt_IsOn(go, "--rna") )
    abc     = esl_alphabet_Create(eslRNA);
  if(!abc) fprintf(stderr, "Defaulting to DNA alphabet.\n"), abc = esl_alphabet_Create(eslDNA);

  LDB("Trying to open %s.\n", cfg->queryfile);
  status = p7_hmmfile_OpenE(cfg->queryfile, NULL, &hfp, errbuf);

  if      (status == eslENOTFOUND) {
    // File just doesn't exist
    p7_Fail("File existence/permissions problem in  == NULLtrying to open query file %s.\n%s\n", cfg->queryfile, errbuf);
  } else if (status == eslOK) {
    //Successfully read HMM file
    qhstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
    if (qhstatus != eslOK) p7_Fail("reading hmm from file %s (%d)\n", cfg->queryfile, qhstatus);
  }
  fprintf(stderr, "abc: %p. ret: %i.\n", abc, qhstatus);

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))              { if ((ofp      = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "--stats"))              { if ((statsfp      = fopen(esl_opt_GetString(go, "--stats"), "w")) == NULL) p7_Fail("Failed to open stats file %s for writing\n",    esl_opt_GetString(go, "-o")); }
  LDB("out files open.\n");


#ifdef HMMER_THREADS
  /* initialize thread data */
  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                                   esl_threads_CPUCount(&ncpus);


  if (ncpus > 0) {
        threadObj = esl_threads_Create(&pipeline_thread);
        queue = esl_workqueue_Create(ncpus * 2);
  }
#endif
  infocnt = (ncpus) ? ncpus : 1;
  LDB("infocnt: %i.\n", infocnt);
  ESL_ALLOC(info, sizeof(*info) * infocnt);
  fprintf(stderr, "Trying to allocate.\n");

  if (! (abc->type == eslRNA || abc->type == eslDNA))
    p7_Fail("Invalid alphabet type in hmm for nhmmer. Expect DNA or RNA\n");
  LDB("alphabet type: %i.\n", abc->type);
  if (qhstatus == eslOK) {
      /* One-time initializations after alphabet <abc> becomes known */
     output_header(ofp, go, cfg->queryfile, cfg->dbfile, ncpus);

     dbfp->abc = abc;

      for (i = 0; i < infocnt; ++i)    {
          info[i].pli    = NULL;
          info[i].th     = NULL;
          info[i].om     = NULL;
          info[i].bg = p7_bg_Create(abc);

#ifdef HMMER_THREADS
          info[i].queue = queue;
#endif
      }
      LDB("About to make blocks.\n");

#ifdef HMMER_THREADS
      for (i = 0; i < ncpus * 2; ++i) {
        {
          block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, abc);
          if (block == NULL)           esl_fatal("Failed to allocate sequence block");

          status = esl_workqueue_Init(queue, block);
          if (status != eslOK)          esl_fatal("Failed to add block to work queue");
        }
      }
#endif
  }


  /* Outer loop: over each query HMM or alignment in <queryfile>. */
  while (qhstatus == eslOK) {
      LDB("Processing one thing.\n");
      P7_PROFILE      *gm      = NULL;
      P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */

      //fprintf(stderr, "qhstatus: %i.\n", qhstatus);


      // Assign HMM max_length
      if      (window_length > 0)     hmm->max_length = window_length;
      else if (window_beta   > 0)     p7_Builder_MaxLength(hmm, window_beta);
      else if (hmm->max_length == -1 ) p7_Builder_MaxLength(hmm, p7_DEFAULT_WINDOW_BETA);


      nquery++;
      resCnt = 0;
      esl_stopwatch_Start(w);
      fprintf(stderr, "dbfp: %p.\n", dbfp);

      /* seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1) {
        {
          if (! esl_sqfile_IsRewindable(dbfp))
            esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", cfg->dbfile);

        }
      }
      LDB("Checked for rewinding needs.\n");

        if ( cfg->firstseq_key != NULL ) { //it's tempting to want to do this once and capture the offset position for future passes, but ncbi files make this non-trivial, so this keeps it general
          sstatus = esl_sqfile_PositionByKey(dbfp, cfg->firstseq_key);
          if (sstatus != eslOK)
            p7_Fail("Failure setting restrictdb_stkey to %d\n", cfg->firstseq_key);
        }


      /* Convert to an optimized model */
      gm = p7_profile_Create (hmm->M, abc);
      om = p7_oprofile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, info->bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; and MSVFilter requires local mode */
      p7_oprofile_Convert(gm, om);                  /* <om> is now p7_LOCAL, multihit */

      scoredata = p7_hmm_ScoreDataCreate(om, NULL);

      LDB("Building info for queue.\n");
      for (i = 0; i < infocnt; ++i) {
          /* Create processing pipeline and hit list */
          info[i].th  = p7_tophits_Create();
          info[i].om = p7_oprofile_Copy(om);
          info[i].pli = p7_pipeline_Create(go, om->M, 100, TRUE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */

          //set method specific --F1, if it wasn't set at command line
          if (!esl_opt_IsOn(go, "--F1") )
              info[i].pli->F1 = 0.02;

          p7_pli_NewModel(info[i].pli, info[i].om, info[i].bg);

          info[i].pli->do_alignment_score_calc = 0;

          info[i].pli->strands = p7_STRAND_BOTH;


          if (  esl_opt_IsUsed(go, "--block_length") )
            info[i].pli->block_length = esl_opt_GetInteger(go, "--block_length");
          else
            info[i].pli->block_length = NHMMER_MAX_RESIDUE_COUNT;

          info[i].scoredata = p7_hmm_ScoreDataClone(scoredata, om->abc->Kp);

#ifdef HMMER_THREADS
          if (ncpus > 0)
            esl_threads_AddThread(threadObj, &info[i]);
#endif
      }

      /* establish the id_lengths data structutre */
      id_length_list = init_id_length(1000);

#ifdef HMMER_THREADS
      if (ncpus > 0)  sstatus = thread_loop    (info, id_length_list, threadObj, queue, dbfp, cfg->firstseq_key, cfg->n_targetseq);
      else            sstatus = serial_loop    (info, id_length_list, dbfp, cfg->firstseq_key, cfg->n_targetseq/*, ssv_watch_master, postssv_watch_master, watch_slave*/);
#else //HMMER_THREADS
        sstatus = serial_loop    (info, id_length_list, dbfp, cfg->firstseq_key, cfg->n_targetseq/*, ssv_watch_master, postssv_watch_master, watch_slave*/);
#endif //HMMER_THREADS


      switch(sstatus) {
        case eslEFORMAT:
          esl_fatal("Parse failed (sequence file %s):\n%s\n",
                    dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
          break;
        case eslEOF:
        case eslOK:
          /* do nothing */
          break;
        default:
          esl_fatal("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);
      }

      //need to re-compute e-values before merging (when list will be sorted)
      if (esl_opt_IsUsed(go, "-Z")) {
    	  resCnt = 1000000*esl_opt_GetReal(go, "-Z");

    	  if ( info[0].pli->strands == p7_STRAND_BOTH)
    	    resCnt *= 2;

      } else
          for (i = 0; i < infocnt; ++i)
            resCnt += info[i].pli->nres;

      for (i = 0; i < infocnt; ++i)
          p7_tophits_ComputeNhmmerEvalues(info[i].th, resCnt, info[i].om->max_length);

      /* merge the results of the search results */
      for (i = 1; i < infocnt; ++i) {
          p7_tophits_Merge(info[0].th, info[i].th);
          p7_pipeline_Merge(info[0].pli, info[i].pli);

          p7_pipeline_Destroy(info[i].pli);
          p7_tophits_Destroy(info[i].th);
          p7_oprofile_Destroy(info[i].om);
      }

      /* Print the results.  */
      p7_tophits_SortBySeqidxAndAlipos(info->th);
      assign_Lengths(info->th, id_length_list);
      p7_tophits_RemoveDuplicates(info->th, info->pli->use_bit_cutoffs);

      p7_tophits_SortBySortkey(info->th);
      p7_tophits_Threshold(info->th, info->pli);
      p7_tophits_PrintHtsSummary(ofp, info->th, hmm->name, hmm->acc);


      //tally up total number of hits and target coverage
      info->pli->n_output = info->pli->pos_output = 0;
      for (i = 0; i < info->th->N; i++) {
          if (info->th->hit[i]->flags & (p7_IS_REPORTED | p7_IS_INCLUDED)) {
              info->pli->n_output++;
              info->pli->pos_output += abs(info->th->hit[i]->dcl[0].jali - info->th->hit[i]->dcl[0].iali) + 1;
          }
      }

      //p7_tophits_Targets(ofp, info->th, info->pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      //p7_tophits_Domains(ofp, info->th, info->pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      esl_stopwatch_Stop(w);

      //p7_pli_Statistics(statsfp, info->pli, w);

//      esl_stopwatch_Display(stdout, ssv_watch_master,     "# SSV time: ");
//      esl_stopwatch_Display(stdout, postssv_watch_master, "# POSTSSV time: ");

//      esl_stopwatch_Include(ssv_watch_master_tot, ssv_watch_master);
//      esl_stopwatch_Include(postssv_watch_master_tot, postssv_watch_master);

      //reset the per-query master stopwatches
//      esl_stopwatch_Start(ssv_watch_master);      esl_stopwatch_Stop(ssv_watch_master);
//      esl_stopwatch_Start(postssv_watch_master);  esl_stopwatch_Stop(postssv_watch_master);

      for (i = 0; i < infocnt; ++i)
        p7_hmm_ScoreDataDestroy(info[i].scoredata);

      p7_hmm_ScoreDataDestroy(scoredata);
      p7_pipeline_Destroy(info->pli);
      p7_tophits_Destroy(info->th);
      p7_oprofile_Destroy(info->om);
      p7_oprofile_Destroy(om);
      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);
      destroy_id_length(id_length_list);
      if (qsq != NULL) esl_sq_Reuse(qsq);

      if (hfp != NULL) {
        qhstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
      }

  } /* end outer loop over queries */

  if (hfp != NULL) {
    switch(qhstatus) {
      case eslEOD:        p7_Fail("read failed, HMM file %s may be truncated?", cfg->queryfile);      break;
      case eslEFORMAT:    p7_Fail("bad file format in HMM file %s",             cfg->queryfile);      break;
      case eslEINCOMPAT:  p7_Fail("HMM file %s contains different alphabets",   cfg->queryfile);      break;
      case eslEOF:        /* do nothing; EOF is what we expect here */                              break;
      default:            p7_Fail("Unexpected error (%d) in reading HMMs from %s", qhstatus, cfg->queryfile);
    }
  }



//  esl_stopwatch_Display(stdout, ssv_watch_master_tot, "# Total SSV time: ");
//  esl_stopwatch_Display(stdout, postssv_watch_master_tot, "# Total POSTSSV time: ");


 /* Terminate outputs - any last words?
   */
  //if (statsfp)      { if (fprintf(statsfp, "#[ok]\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

  /* Cleanup - prepare for successful exit
   */
  for (i = 0; i < infocnt; ++i)
    p7_bg_Destroy(info[i].bg);

#ifdef HMMER_THREADS
  if (ncpus > 0) {
      esl_workqueue_Reset(queue);
      {
        while (esl_workqueue_Remove(queue, (void **) &block) == eslOK) {
          esl_sq_DestroyBlock(block);
        }
      }
      esl_workqueue_Destroy(queue);
      esl_threads_Destroy(threadObj);
  }
#endif

  free(info);



  if (hfp)     p7_hmmfile_Close(hfp);

  if (builder) p7_builder_Destroy(builder);
  if (qsq)     esl_sq_Destroy(qsq);

  esl_sqfile_Close(dbfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);


  if (ofp != stdout) fclose(ofp);
  if (statsfp != stderr) fclose(statsfp);

  return eslOK;

 ERROR:
   if (hfp)     p7_hmmfile_Close(hfp);

   if (builder) p7_builder_Destroy(builder);
   if (qsq)     esl_sq_Destroy(qsq);

   if (ofp != stdout) fclose(ofp);

   return eslFAIL;
}

//TODO: MPI code needs to be added here
static int
serial_loop(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_SQFILE *dbfp, char *firstseq_key, int n_targetseqs/*, ESL_STOPWATCH *ssv_watch_master, ESL_STOPWATCH *postssv_watch_master, ESL_STOPWATCH *watch_slave*/)
{

  int      wstatus = eslOK;
  int      status = eslOK;
  int seq_id = 0;
  ESL_SQ   *dbsq   =  esl_sq_CreateDigital(info->om->abc);
#ifdef eslAUGMENT_ALPHABET
  ESL_SQ   *dbsq_revcmp = NULL;


  if (dbsq->abc->complement != NULL)
    dbsq_revcmp =  esl_sq_CreateDigital(info->om->abc);
#endif /*eslAUGMENT_ALPHABET*/

  wstatus = esl_sqio_ReadWindow(dbfp, 0, info->pli->block_length, dbsq);

  while (wstatus == eslOK && (n_targetseqs==-1 || seq_id < n_targetseqs) ) {
      dbsq->idx = seq_id;
      p7_pli_NewSeq(info->pli, dbsq);
      ESL_REALLOC(dbsq->acc, dbsq->n + 1);
      if((status = esl_abc_Textize(dbsq->abc, dbsq->dsq, dbsq->n, dbsq->acc) != eslOK))
        p7_Fail("Error writing alphabetic sequence from binary: %i.\n", status);

      if (info->pli->strands != p7_STRAND_BOTTOMONLY) {

        info->pli->nres -= dbsq->C; // to account for overlapping region of windows
        p7_Pipeline_LongTarget(info->pli, info->om, info->scoredata, info->bg, info->th, info->pli->nseqs, dbsq, p7_NOCOMPLEMENT, NULL, NULL, NULL/*, ssv_watch_master, postssv_watch_master, watch_slave*/);
        p7_pipeline_Reuse(info->pli); // prepare for next search

      } else {
        info->pli->nres -= dbsq->n;
      }
#ifdef eslAUGMENT_ALPHABET
      //reverse complement
      if (info->pli->strands != p7_STRAND_TOPONLY && dbsq->abc->complement != NULL )
      {
          esl_sq_Copy(dbsq, dbsq_revcmp);
          esl_sq_ReverseComplement(dbsq_revcmp);
          ESL_REALLOC(dbsq_revcmp->acc, dbsq_revcmp->n + 1);
          if((status = esl_abc_Textize(dbsq_revcmp->abc, dbsq_revcmp->dsq, dbsq_revcmp->n, dbsq_revcmp->acc) != eslOK))
            p7_Fail("Error writing alphabetic sequence from binary: %i.\n", status);
          p7_Pipeline_LongTarget(info->pli, info->om, info->scoredata, info->bg, info->th, info->pli->nseqs, dbsq_revcmp, p7_COMPLEMENT, NULL, NULL, NULL/*, ssv_watch_master, postssv_watch_master, watch_slave*/);
          p7_pipeline_Reuse(info->pli); // prepare for next search

          info->pli->nres += dbsq_revcmp->W;

      }
#endif /*eslAUGMENT_ALPHABET*/

      wstatus = esl_sqio_ReadWindow(dbfp, info->om->max_length, info->pli->block_length, dbsq);
      if (wstatus == eslEOD) { // no more left of this sequence ... move along to the next sequence.
          add_id_length(id_length_list, dbsq->idx, dbsq->L);

          info->pli->nseqs++;
          esl_sq_Reuse(dbsq);
          wstatus = esl_sqio_ReadWindow(dbfp, 0, info->pli->block_length, dbsq);

          seq_id++;

      }
    }


  if (dbsq) esl_sq_Destroy(dbsq);
  if (dbsq_revcmp) esl_sq_Destroy(dbsq_revcmp);

  ERROR: p7_Fail("Could not allocate memory.\n");
  return wstatus;

}


#ifdef HMMER_THREADS
static int
thread_loop(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp, char *firstseq_key, int n_targetseqs)
{

  int          i;
  int          status  = eslOK;
  int          sstatus = eslOK;
  int          eofCount = 0;
  ESL_SQ_BLOCK *block;
  void         *newBlock;
  int          seqid = -1;

  ESL_SQ      *tmpsq = esl_sq_CreateDigital(info->om->abc);
  int          abort = FALSE; // in the case n_targetseqs != -1, a block may get abbreviated


  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
  ((ESL_SQ_BLOCK *)newBlock)->complete = TRUE;

  /* Main loop: */
  while (sstatus == eslOK  ) {
      block = (ESL_SQ_BLOCK *) newBlock;

      if (abort) {
        block->count = 0;
        sstatus = eslEOF;
      } else {
        sstatus = esl_sqio_ReadBlock(dbfp, block, info->pli->block_length, n_targetseqs, TRUE);
      }

      block->first_seqidx = info->pli->nseqs;
      seqid = block->first_seqidx;
      for (i=0; i<block->count; i++) {
        block->list[i].idx = seqid;
        add_id_length(id_length_list, seqid, block->list[i].L);
        seqid++;

        if (   seqid == n_targetseqs // hit the sequence target
            && ( i<block->count-1 ||  block->complete ) // and either it's not the last sequence (so it's complete), or its complete
        ) {
          abort = TRUE;
          block->count = i+1;
          break;
        }
      }
      info->pli->nseqs += block->count  - ((abort || block->complete) ? 0 : 1);// if there's an incomplete sequence read into the block wait to count it until it's complete.


      if (sstatus == eslEOF) {
          if (eofCount < esl_threads_GetWorkerCount(obj)) sstatus = eslOK;
          ++eofCount;
      } else if (!block->complete ) {
          // The final sequence on the block was an incomplete window of the active sequence,
          // so our next read will need a copy of it to correctly deal with overlapping
          // regions. We capture a copy of the sequence here before sending it off to the
          // pipeline to avoid odd race conditions that can occur otherwise.
          // Copying the entire sequence isn't really necessary, and is a bit heavy-
          // handed. Could accelerate if this proves to have any notable impact on speed.
          esl_sq_Copy(block->list + (block->count - 1) , tmpsq);
      }


      if (sstatus == eslOK) {

          status = esl_workqueue_ReaderUpdate(queue, block, &newBlock);
          if (status != eslOK) esl_fatal("Work queue reader failed");

          //newBlock needs all this information so the next ReadBlock call will know what to do
          ((ESL_SQ_BLOCK *)newBlock)->complete = block->complete;
          if (!block->complete) {
              // Push the captured copy of the previously-read sequence into the new block,
              // in preparation for ReadWindow  (double copy ... slower than necessary)
              esl_sq_Copy(tmpsq, ((ESL_SQ_BLOCK *)newBlock)->list);

              if (  ((ESL_SQ_BLOCK *)newBlock)->list->n < info->om->max_length ) {
                //no reason to search the final partial sequence on the block, as the next block will search this whole chunk
                ((ESL_SQ_BLOCK *)newBlock)->list->C = ((ESL_SQ_BLOCK *)newBlock)->list->n;
                (((ESL_SQ_BLOCK *)newBlock)->count)--;
              } else {
                ((ESL_SQ_BLOCK *)newBlock)->list->C = info->om->max_length;
              }

          }
      }
  }


  status = esl_workqueue_ReaderUpdate(queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue reader failed");

  if (sstatus == eslEOF) {
      /* wait for all the threads to complete */
      esl_threads_WaitForFinish(obj);
      esl_workqueue_Complete(queue);
    }

  esl_sq_Destroy(tmpsq);

  return sstatus;
}

static void
pipeline_thread(void *arg)
{
  int i;
  int status;
  int workeridx;
  WORKER_INFO   *info;
  ESL_THREADS   *obj;
  ESL_SQ_BLOCK  *block = NULL;
  void          *newBlock;

  impl_Init();

  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue worker failed");


  /* loop until all blocks have been processed */
  block = (ESL_SQ_BLOCK *) newBlock;

  while (block->count > 0)
  {
      /* Main loop: */
      for (i = 0; i < block->count; ++i)
    {
      ESL_SQ *dbsq = block->list + i;

      p7_pli_NewSeq(info->pli, dbsq);
      ESL_REALLOC(dbsq->acc, dbsq->n + 1);
      if((status = esl_abc_Textize(dbsq->abc, dbsq->dsq, dbsq->n, dbsq->acc) != eslOK))
        p7_Fail("Error writing alphabetic sequence from binary: %i.\n", status);

      if (info->pli->strands != p7_STRAND_BOTTOMONLY) {
        info->pli->nres -= dbsq->C; // to account for overlapping region of windows

        p7_Pipeline_LongTarget(info->pli, info->om, info->scoredata, info->bg, info->th, block->first_seqidx + i, dbsq, p7_NOCOMPLEMENT, NULL, NULL, NULL/*, NULL, NULL, NULL*/);
        p7_pipeline_Reuse(info->pli); // prepare for next search

      } else {
        info->pli->nres -= dbsq->n;
      }

#ifdef eslAUGMENT_ALPHABET
      //reverse complement
      if (info->pli->strands != p7_STRAND_TOPONLY && dbsq->abc->complement != NULL)
      {
          esl_sq_ReverseComplement(dbsq);
          ESL_REALLOC(dbsq->acc, dbsq->n + 1);
          if((status = esl_abc_Textize(dbsq->abc, dbsq->dsq, dbsq->n, dbsq->acc) != eslOK))
            p7_Fail("Error writing alphabetic sequence from binary: %i.\n", status);
          p7_Pipeline_LongTarget(info->pli, info->om, info->scoredata, info->bg, info->th, block->first_seqidx + i, dbsq, p7_COMPLEMENT, NULL, NULL, NULL/*, NULL, NULL, NULL*/);
          p7_pipeline_Reuse(info->pli); // prepare for next search

          info->pli->nres += dbsq->W;
      }

#endif /*eslAUGMENT_ALPHABET*/

    }

      status = esl_workqueue_WorkerUpdate(info->queue, block, &newBlock);
      if (status != eslOK) esl_fatal("Work queue worker failed");

      block = (ESL_SQ_BLOCK *) newBlock;

  }

  status = esl_workqueue_WorkerUpdate(info->queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  esl_threads_Finished(obj, workeridx);
  return;
  ERROR: p7_Fail("Could not allocate memory.\n");
}



#endif   /* HMMER_THREADS */



/* helper functions for tracking id_lengths */

static ID_LENGTH_LIST *
init_id_length( int size )
{
  int status;
  ID_LENGTH_LIST *list;

  ESL_ALLOC (list, sizeof(ID_LENGTH_LIST));
  list->count = 0;
  list->size  = size;
  list->id_lengths = NULL;

  ESL_ALLOC (list->id_lengths, size * sizeof(ID_LENGTH));

  return list;

ERROR:
  return NULL;
}

static void
destroy_id_length( ID_LENGTH_LIST *list )
{

  if (list != NULL) {
    if (list->id_lengths != NULL) free (list->id_lengths);
    free (list);
  }

}


static int
add_id_length(ID_LENGTH_LIST *list, int id, int L)
{
   int status;

   if (list->count > 0 && list->id_lengths[list->count-1].id == id) {
     // the last time this gets updated, it'll have the sequence's actual length
     list->id_lengths[list->count-1].length = L;
   } else {

     if (list->count == list->size) {
       list->size *= 10;
       ESL_REALLOC(list->id_lengths, list->size * sizeof(ID_LENGTH));
     }

     list->id_lengths[list->count].id     = id;
     list->id_lengths[list->count].length = L;

     list->count++;
   }
   return eslOK;

ERROR:
   return status;
}

static int
assign_Lengths(P7_TOPHITS *th, ID_LENGTH_LIST *id_length_list) {

  int i;
  int j = 0;
  for (i=0; i<th->N; i++) {
    while (th->hit[i]->seqidx != id_length_list->id_lengths[j].id) { j++;   }
    th->hit[i]->dcl[0].ad->L = id_length_list->id_lengths[j].length;
  }

  return eslOK;
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/



