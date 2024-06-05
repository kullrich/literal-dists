#include "kseq.h"
#include <getopt.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <zlib.h>
#include <string.h>

#define VERSION "0.0.3"
#define EXENAME "literal-dists"
#define GITHUB_URL "https://github.com/kullrich/literal-dists"

const int MAX_SEQ = 100000;
const char IGNORE_CHAR = '.';
const int PHRED = 33;
const char GAP_CHAR0 = '.';
const char GAP_CHAR1 = '-';
const char GAP_CHAR2 = 'N';
const char GAP_CHAR3 = 'X';

// define score matrix
char *nucnuc[324] = {
  "AA", "AC", "AG", "AT", "AR", "AY", "AS", "AW", "AK", "AM", "AB", "AD", "AH", "AV", "A.", "A-", "AN", "AX",
  "CA", "CC", "CG", "CT", "CR", "CY", "CS", "CW", "CK", "CM", "CB", "CD", "CH", "CV", "C.", "C-", "CN", "CX",
  "GA", "GC", "GG", "GT", "GR", "GY", "GS", "GW", "GK", "GM", "GB", "GD", "GH", "GV", "G.", "G-", "GN", "GX",
  "TA", "TC", "TG", "TT", "TR", "TY", "TS", "TW", "TK", "TM", "TB", "TD", "TH", "TV", "T.", "T-", "TN", "TX",
  "RA", "RC", "RG", "RT", "RR", "RY", "RS", "RW", "RK", "RM", "RB", "RD", "RH", "RV", "R.", "R-", "RN", "RX",
  "YA", "YC", "YG", "YT", "YR", "YY", "YS", "YW", "YK", "YM", "YB", "YD", "YH", "YV", "Y.", "Y-", "YN", "YX",
  "SA", "SC", "SG", "ST", "SR", "SY", "SS", "SW", "SK", "SM", "SB", "SD", "SH", "SV", "S.", "S-", "SN", "SX",
  "WA", "WC", "WG", "WT", "WR", "WY", "WS", "WW", "WK", "WM", "WB", "WD", "WH", "WV", "W.", "W-", "WN", "WX",
  "KA", "KC", "KG", "KT", "KR", "KY", "KS", "KW", "KK", "KM", "KB", "KD", "KH", "KV", "K.", "K-", "KN", "KX",
  "MA", "MC", "MG", "MT", "MR", "MY", "MS", "MW", "MK", "MM", "MB", "MD", "MH", "MV", "M.", "M-", "MN", "MX",
  "BA", "BC", "BG", "BT", "BR", "BY", "BS", "BW", "BK", "BM", "BB", "BD", "BH", "BV", "B.", "B-", "BN", "BX",
  "DA", "DC", "DG", "DT", "DR", "DY", "DS", "DW", "DK", "DM", "DB", "DD", "DH", "DV", "D.", "D-", "DN", "DX",
  "HA", "HC", "HG", "HT", "HR", "HY", "HS", "HW", "HK", "HM", "HB", "HD", "HH", "HV", "H.", "H-", "HN", "HX",
  "VA", "VC", "VG", "VT", "VR", "VY", "VS", "VW", "VK", "VM", "VB", "VD", "VH", "VV", "V.", "V-", "VN", "VX",
  ".A", ".C", ".G", ".T", ".R", ".Y", ".S", ".W", ".K", ".M", ".B", ".D", ".H", ".V", "..", ".-", ".N", ".X",
  "-A", "-C", "-G", "-T", "-R", "-Y", "-S", "-W", "-K", "-M", "-B", "-D", "-H", "-V", "-.", "--", "-N", "-X",
  "NA", "NC", "NG", "NT", "NR", "NY", "NS", "NW", "NK", "NM", "NB", "ND", "NH", "NV", "N.", "N-", "NN", "NX",
  "XA", "XC", "XG", "XT", "XR", "XY", "XS", "XW", "XK", "XM", "XB", "XD", "XH", "XV", "X.", "X-", "XN", "XX"
};
double nucnucscore[324] = {
  0.0, 1.0, 1.0, 1.0, 0.5, 1.0, 1.0, 0.5, 1.0, 0.5, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  1.0, 0.0, 1.0, 1.0, 1.0, 0.5, 0.5, 1.0, 1.0, 0.5, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  1.0, 1.0, 0.0, 1.0, 0.5, 1.0, 0.5, 1.0, 0.5, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  1.0, 1.0, 1.0, 0.0, 1.0, 0.5, 1.0, 0.5, 0.5, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  0.5, 1.0, 0.5, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  1.0, 0.5, 1.0, 0.5, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  1.0, 0.5, 0.5, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  0.5, 1.0, 1.0, 0.5, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  1.0, 1.0, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
};

KSEQ_INIT(gzFile, gzread)

//------------------------------------------------------------------------
size_t distance(const char* restrict a, const char* restrict b, size_t L, int remgaps, int usesite[])
{
  size_t diff = 0;
  for (size_t i = 0; i < L; i++) {
    if (remgaps) {
      if (usesite[i] == 0) {
        continue;
      }
    }
    if (a[i] != b[i] && a[i] != IGNORE_CHAR && b[i] != IGNORE_CHAR) {
      diff++;
    }
  }
  return diff;
}

void literaldistance(const char* restrict a, const char* restrict b, size_t L, double scorematrix[][57], double about[], int remgaps, int usesite[])
{
  double abdist = 0.0;
  double abnsites = 0.0;
  double abscore = 0.0;
  double abscoreij = 0.0;
  int ai, bi;
  abnsites = L;
  for (size_t i = 0; i < L; i++) {
    if (remgaps) {
      if (usesite[i] == 0) {
        abnsites-=1.0;
        continue;
      }
    }
    if (a[i] != IGNORE_CHAR && b[i] != IGNORE_CHAR) {
      ai = a[i];
      bi = b[i];
      abscoreij=scorematrix[ai-PHRED][bi-PHRED];
      //printf("i: %zu a[i]: %c %d b[i]: %c %d score: %f\n", i, a[i], ai, b[i], bi, abscore);
      if (abscoreij>=0) {
        //abnsites+=1.0;
        abscore+=abscoreij;
      } else {
        abnsites-=1.0;
      }
    } else {
      abnsites-=1.0;
    }
  }
  abdist=abscore/abnsites;
  //printf("abdist: %f abscore: %f abnsites: %f\n", abdist, abscore, abnsites);
  about[0]=abdist;
  about[1]=abscore;
  about[2]=abnsites;
}

//------------------------------------------------------------------------
void show_help(int retcode)
{
  FILE* out = (retcode == EXIT_SUCCESS ? stdout : stderr);
  
  static const char str[] = {
      "SYNOPSIS\n  Pairwise literal distance matrix from a FASTA alignment\n"
      "USAGE\n  %s [options] -f alignment.fasta[.gz] > matrix.tsv\n"
      "or\n"
      "USAGE\n  %s [options] < alignment.fasta[.gz] > matrix.tsv\n"
      "OPTIONS\n"
      "  -h\tShow this help\n"
      "  -v\tPrint version and exit\n"
      "  -t\tUse original snp-dists distance\n"
      "  -q\tQuiet mode; do not print progress information\n"
      "  -a\tCount all differences not just [AGTC]\n"
      "  -k\tKeep case, don't uppercase all letters\n"
      "  -m\tOutput MOLTEN instead of TSV\n"
      "  -c\tUse comma instead of tab in output\n"
      "  -b\tBlank top left corner cell\n"
      "  -i\tOutput used sites, scores and gap sites to TSV\n"
      "  -g\tSkip gap sites if gap frequency is met, gap sites [.-NX]\n"
      "  -z\tGap frequency [default: 0.5]\n"
      "  -f\tInput FASTA file\n"
      "URL\n  %s\n"};
  fprintf(out, str, EXENAME, VERSION, GITHUB_URL);
  exit(retcode);
}

//------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  // parse command line parameters
  int opt, quiet = 0, csv = 0, corner = 1, sitestotsv = 0, allchars = 0, remgaps = 0, keepcase = 0, molten = 0, snpdists = 0; fasta_file = 0;
  double gapfreq = 0.5;
  int gapsites = 0;
  while ((opt = getopt(argc, argv, "htqcakmbigz:vf")) != -1) {
    switch (opt) {
      case 'h': show_help(EXIT_SUCCESS); break;
      case 't': snpdists = 1; break;
      case 'q': quiet = 1; break;
      case 'c': csv = 1; break;
      case 'a': allchars = 1; break;
      case 'k': keepcase = 1; break;
      case 'm': molten = 1; break;
      case 'b': corner = 0; break;
      case 'i': sitestotsv = 1; break;
      case 'g': remgaps = 1; break;
      case 'z': gapfreq = strtod(optarg, NULL); break;
      case 'v': printf("%s %s\n", EXENAME, VERSION); exit(EXIT_SUCCESS);
      case 'f': if (strcmp(optarg, "-f") == 0) { fasta_file = 1; } else { fasta_file = 0; } break;
      default: show_help(EXIT_FAILURE);
    }
  }
  // say hello
  if (!quiet) {
    fprintf(stderr, "This is %s %s\n", EXENAME, VERSION);
  }
  // open filename via libz
  gzFile fp;
  if (fasta_file) {
    if (optind >= argc) {
      show_help(EXIT_FAILURE);
      return 0;
    }
    const char* fasta = argv[optind];
    fp = gzopen(fasta, "r");
    if (!fp) {
      fprintf(stderr, "ERROR: Could not open filename '%s'\n", fasta);
      exit(EXIT_FAILURE);
    }
  } else {
    fp = gzdopen(fileno(stdin), "r");
    if (!fp) {
      fprintf(stderr, "ERROR: Could not open stdin\n");
      exit(EXIT_FAILURE);
    }
  }
  // load all the sequences
  char** seq = calloc(MAX_SEQ, sizeof(char*));
  char** name = calloc(MAX_SEQ, sizeof(char*));
  ssize_t l, N=0, L=-1;
  kseq_t* kseq = kseq_init(fp);
  while ((l = kseq_read(kseq)) >= 0) {
    // first sequence
    if (L < 0) {
      L = l;
    }
    // not first sequence - so ensure length matches first one
    if (l != L) {
      fprintf(stderr,
              "ERROR: sequence #%zu '%s' has length %zu but expected %zu\n",
              N + 1, kseq->name.s, l, L);
      exit(EXIT_FAILURE);
    }
    // have we exceeded the number of sequences we can handle?
    if (N >= MAX_SEQ) {
      fprintf(stderr,
              "ERROR: %s can only handle %d sequences at most. Please change MAX_SEQ "
              "and recompile.\n",
              EXENAME, MAX_SEQ);
      exit(EXIT_FAILURE);
    }
    // save the sequence and name
    seq[N] = (char*)calloc(kseq->seq.l + 1, sizeof(char));
    strcpy(seq[N], kseq->seq.s);
    name[N] = (char*)calloc(kseq->name.l + 1, sizeof(char));
    strcpy(name[N], kseq->name.s);
    // uppercase all sequences
    if (!keepcase) {
      for (char* s = seq[N]; *s; s++) {
        *s = toupper(*s);
      }
    }
    // clean the sequence depending on -a option
    if (!allchars) {
      for (char* s = seq[N]; *s; s++) {
        if(snpdists){
          if (*s != 'A' && *s != 'T' && *s != 'C' && *s != 'G') {
            *s = IGNORE_CHAR;
          }
        } else {
          if (*s != 'A' && *s != 'C' && *s != 'G' && *s != 'T' && *s != 'R' && *s != 'Y' && *s != 'S' && *s != 'W' && *s != 'K' && *s != 'M' && *s != 'B' && *s != 'D' && *s != 'H' && *s != 'V' && *s != '-' && *s != 'X' && *s != 'N') {
            *s = IGNORE_CHAR;
          }
        }
      }
    }
    // keep track of how many we have
    N++;
  }
  kseq_destroy(kseq);
  gzclose(fp);
  if (N < 1) {
    fprintf(stderr, "ERROR: file contained no sequences\n");
    exit(EXIT_FAILURE);
  }
  if (!quiet) {
    fprintf(stderr, "Read %zu sequences of length %zu\n", N, L);
  }
  // output TSV or CSV
  char sep = csv ? ',' : '\t';
  // create score matrix
  double scorematrix[57][57] = {{0.0}};
  if (!snpdists) {
    for (int i=0;i < 324;i++){
        //printf("i: %d %d %c %d %c %f\n", i, nucnuc[i][0]-PHRED, nucnuc[i][0], nucnuc[i][1]-PHRED, nucnuc[i][1], nucnucscore[i]);
        scorematrix[nucnuc[i][0]-PHRED][nucnuc[i][1]-PHRED]=nucnucscore[i];
    }
  }
  // create gap frequencies
  double gfreq[L];
  int usesite[L];
  if (remgaps) {
    for (int i = 0; i < L; i++) {
      gfreq[i] = 0.0;
      usesite[i] = 1;
    }
    for (int j = 0; j < N; j++) {
      const char* restrict seq_j;
      seq_j = seq[j];
      for (int i = 0; i < L; i++) {
        if (seq_j[i] == GAP_CHAR0 || seq_j[i] == GAP_CHAR1 || seq_j[i] == GAP_CHAR2 || seq_j[i] == GAP_CHAR3) {
          gfreq[i] += 1;  
        }
      }
    }
    for (int i = 0; i < L; i++) {
      gfreq[i] = gfreq[i]/N;
      // printf("gapfreq site: %d = %f\n", i, gfreq[i]);
      if (gfreq[i] != 0) {
        if (gfreq[i] >= gapfreq) {
          usesite[i] = 0;
          gapsites += 1;
          // printf("skip site: %d\n", i);
        }
      }
    }
  }
  if (molten) {
    // "molten" format, one row per pair
    if (!quiet) {
      if (remgaps) {
        printf("#seq1%cseq2%cdist%cscore%cnsites%cgapsites\n", sep, sep, sep, sep, sep);
      } else {
        printf("#seq1%cseq2%cdist%cscore%cnsites\n", sep, sep, sep, sep);
      }
    }
    for (int j = 0; j < N; j++) {
      for (int i = 0; i < N; i++) {
        if(snpdists){
            size_t d = distance(seq[j], seq[i], L, remgaps, usesite);
            printf("%s%c%s%c%zu\n", name[j], sep, name[i], sep, d);
        } else{
          double about[3];
          literaldistance(seq[j], seq[i], L, scorematrix, about, remgaps, usesite);
          //printf("abdist: %f abscore: %f abnsites: %f\n", about[0], about[1], about[2]);
          if (remgaps) {
            printf("%s%c%s%c%f%c%f%c%f%c%d\n", name[j], sep, name[i], sep, about[0], sep, about[1], sep, about[2], sep, gapsites);
          } else {
            printf("%s%c%s%c%f%c%f%c%f\n", name[j], sep, name[i], sep, about[0], sep, about[1], sep, about[2]);
          }
        }
      }
    }
  }
  else {
    // regular TSV matrix output
    // header seq
    if (corner){
      if(snpdists){
        printf("snp-dists 0.7.0");
      } else{
        printf("%s %s", EXENAME, VERSION);
      }
    }
    for (int j = 0; j < N; j++) {
      printf("%c%s", sep, name[j]);
    }
    printf("\n");
    // Output the distance matrix to stdout
    // (does full matrix, wasted computation i know)
    for (int j = 0; j < N; j++) {
      printf("%s", name[j]);
      for (int i = 0; i < N; i++) {
        if(snpdists){
          size_t d = distance(seq[j], seq[i], L, remgaps, usesite);
          printf("%c%zu", sep, d);
        } else {
          double about[3];
          literaldistance(seq[j], seq[i], L, scorematrix, about, remgaps, usesite);
          if(sitestotsv) {
            printf("%c%f/%f/%f/%d", sep, about[0], about[1], about[2], gapsites);
          } else {
            printf("%c%f", sep, about[0]);
          }
        }
      }
      printf("\n");
    }
  }
  // free memory
  for (int k = 0; k < N; k++) {
    free(seq[k]);
    free(name[k]);
  }
  free(seq);
  free(name);
  return 0;
}

//------------------------------------------------------------------------
