//nog enkele definities:
#ifdef PQ

#define __Q_CON

#endif

#ifdef PQG

#define __Q_CON
#define __G_CON

#endif

#ifdef PQGT1

#define __Q_CON
#define __G_CON
#define __T1_CON

#endif

#ifdef PQGT2

#define __Q_CON
#define __G_CON
#define __T2_CON

#endif

#ifdef PQGT

#define __Q_CON
#define __G_CON
#define __T1_CON
#define __T2_CON

#endif

#ifdef PQGT2P

#define __Q_CON
#define __G_CON
#define __T2P_CON

#endif

#ifdef PQGTP

#define __Q_CON
#define __G_CON
#define __T1_CON
#define __T2P_CON

#endif

#include "lapack.h"
#include "Matrix.h"
#include "Vector.h"
#include "TPM.h"
#include "SPM.h"
#include "PHM.h"
#include "DPM.h"
#include "PPHM.h"
#include "T2PM.h"

#include "SUP.h"
#include "EIG.h"

/* vim: set ts=3 sw=3 expandtab :*/
