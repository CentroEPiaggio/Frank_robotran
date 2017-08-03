#ifndef __c1_Frank_lqr_h__
#define __c1_Frank_lqr_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef struct_struct_ZwdsKLYK9S2KrT5tyDKoxE_tag
#define struct_struct_ZwdsKLYK9S2KrT5tyDKoxE_tag

struct struct_ZwdsKLYK9S2KrT5tyDKoxE_tag
{
  real_T dpt[57];
  real_T m[17];
};

#endif                                 /*struct_struct_ZwdsKLYK9S2KrT5tyDKoxE_tag*/

#ifndef typedef_c1_struct_ZwdsKLYK9S2KrT5tyDKoxE
#define typedef_c1_struct_ZwdsKLYK9S2KrT5tyDKoxE

typedef struct struct_ZwdsKLYK9S2KrT5tyDKoxE_tag
  c1_struct_ZwdsKLYK9S2KrT5tyDKoxE;

#endif                                 /*typedef_c1_struct_ZwdsKLYK9S2KrT5tyDKoxE*/

#ifndef struct_s3P2Hgr8IQaO66xIl8CEZeB
#define struct_s3P2Hgr8IQaO66xIl8CEZeB

struct s3P2Hgr8IQaO66xIl8CEZeB
{
  real_T P[3];
  real_T R[9];
  real_T V[3];
  real_T OM[3];
  real_T A[3];
  real_T OMP[3];
  real_T J[102];
};

#endif                                 /*struct_s3P2Hgr8IQaO66xIl8CEZeB*/

#ifndef typedef_c1_s3P2Hgr8IQaO66xIl8CEZeB
#define typedef_c1_s3P2Hgr8IQaO66xIl8CEZeB

typedef struct s3P2Hgr8IQaO66xIl8CEZeB c1_s3P2Hgr8IQaO66xIl8CEZeB;

#endif                                 /*typedef_c1_s3P2Hgr8IQaO66xIl8CEZeB*/

#ifndef typedef_SFc1_Frank_lqrInstanceStruct
#define typedef_SFc1_Frank_lqrInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c1_sfEvent;
  boolean_T c1_isStable;
  boolean_T c1_doneDoubleBufferReInit;
  uint8_T c1_is_active_c1_Frank_lqr;
  c1_struct_ZwdsKLYK9S2KrT5tyDKoxE c1_rob_str;
  real_T c1_Ts;
} SFc1_Frank_lqrInstanceStruct;

#endif                                 /*typedef_SFc1_Frank_lqrInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c1_Frank_lqr_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c1_Frank_lqr_get_check_sum(mxArray *plhs[]);
extern void c1_Frank_lqr_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
