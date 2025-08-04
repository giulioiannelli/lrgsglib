#include <stdio.h>
#include <stddef.h>
#include <inttypes.h>

#ifndef __CUSTOM_H_INC__
#define __CUSTOM_H_INC__
/* preprocessors functions */
#define DUMP(x) (#x)
#define UNUSED(x) (void)(x)
#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))
#define ARRAY_SIZE_2D(arr) (sizeof(arr) / sizeof((arr)[0][0]))
/* files extension */
#define BINX ".bin"
#define TXTX ".txt"
/*basic string formatters*/
#define __ " "
#define _C ","
#define _E "="
#define _H "-"
#define _N "\n"
#define _S "/"
#define _T "\t"
#define _U "_"
#define _V "|"
#define _W "!"
/* composite string formatter */
#define __V __ _V     /*                                                 " |" */
#define _N_ _N __     /*                                               "\\n " */
#define _NN _N _N     /*                                             "\\n\\n" */
#define _NT _N _T     /*                                             "\\n\\t" */
#define _T_ _T __     /*                                               "\\t " */
#define _TT _T _T     /*                                             "\\t\\t" */
#define _TV _T _V     /*                                               "\\t|" */
#define _VT _V _T     /*                                               "|\\t" */
#define _WT _W _T     /*                                               "!\\t" */
#define _N_V _N_ _V   /*                                              "\\n |" */
#define _N_W _N_ _W   /*                                              "\\n |" */
#define _N_VT _N_ _VT /*                                           "\\n |\\t" */
#define _N_WT _N_ _WT /*                                      "\\n !\\t |\\t" */
#define _NT_V _NT __V /*                                           "\\n\\t |" */
#define _NT_V_ _NT __V __ /*                                      "\\n\\t | " */
/* log file formatters */
#define _HLIFI "{+}"
#define _LLPFI "[?]"
#define _LLPWI "[!]"
#define _LLPII "[+]"
#define _LLPEI "[*]"
#define _LLPPI "[•]"
#define _LOGON "| ON |"
#define _LOGOF "| OFF|"
/* stderr messages */
#define __ERR "ERROR: "
#define __WRN "WARNING: "
#define MSG_ERR_FGETS __ERR "(fgets)"
#define MSG_ERR_FREAD __ERR "(fread)"
#define MSG_ERR_POPEN __ERR "(popen)"
#define MSG_ERR_SCNU32 __ERR "strtou32 attempting to scan uinteger > UINT32_MAX"
#define MSG_WRN_SCNU32 __WRN "strtou32 found extra characters after number."
/* string lengths */
#define STRL16 16
#define STRL64 64
#define STRL128 128
#define STRL256 256
#define STRL512 512
#define STRL1024 1024

typedef double *double_p;
typedef size_t *size_tp;
typedef int8_t *spin_tp;
typedef struct {
    uint64_t u, v;
    double w;
} Edge;
typedef Edge* Edges;

typedef struct {
    size_t *neighbors;
    double *weights;
} NodeEdges;
typedef NodeEdges* NodesEdges;


#define LRGSG_ROOT getenv("LRGSG_ROOT")
#define LRGSG_BUILD getenv("LRGSG_BUILD")
#define LRGSG_DATA getenv("LRGSG_DATA")
#define LRGSG_IPYNB getenv("LRGSG_IPYNB")
#define LRGSG_SRC getenv("LRGSG_SRC")
#define LRGSG_LIB getenv("LRGSG_LIB")
#define LRGSG_LIB_CCORE getenv("LRGSG_LIB_CCORE")
#define LRGSG_LIB_CONFIG getenv("LRGSG_LIB_CONFIG")
#define LRGSG_LIB_GT_PATCHES getenv("LRGSG_LIB_GT_PATCHES")
#define LRGSG_LIB_NX_PATCHES getenv("LRGSG_LIB_NX_PATCHES")
#define LRGSG_LIB_STOCPROC getenv("LRGSG_LIB_STOCPROC")
#define LRGSG_TEST getenv("LRGSG_TEST")
#define LRGSG_TOOLS getenv("LRGSG_TOOLS")

#endif /* __CUSTOM_H_INC__ */