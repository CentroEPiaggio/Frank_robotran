#ifndef __c2_Frank_feedback_h__
#define __c2_Frank_feedback_h__

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

#ifndef typedef_c2_struct_ZwdsKLYK9S2KrT5tyDKoxE
#define typedef_c2_struct_ZwdsKLYK9S2KrT5tyDKoxE

typedef struct struct_ZwdsKLYK9S2KrT5tyDKoxE_tag
  c2_struct_ZwdsKLYK9S2KrT5tyDKoxE;

#endif                                 /*typedef_c2_struct_ZwdsKLYK9S2KrT5tyDKoxE*/

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

#ifndef typedef_c2_s3P2Hgr8IQaO66xIl8CEZeB
#define typedef_c2_s3P2Hgr8IQaO66xIl8CEZeB

typedef struct s3P2Hgr8IQaO66xIl8CEZeB c2_s3P2Hgr8IQaO66xIl8CEZeB;

#endif                                 /*typedef_c2_s3P2Hgr8IQaO66xIl8CEZeB*/

#ifndef typedef_SFc2_Frank_feedbackInstanceStruct
#define typedef_SFc2_Frank_feedbackInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c2_sfEvent;
  boolean_T c2_isStable;
  boolean_T c2_doneDoubleBufferReInit;
  uint8_T c2_is_active_c2_Frank_feedback;
  c2_struct_ZwdsKLYK9S2KrT5tyDKoxE c2_rob_str;
  real_T c2_Ts;
} SFc2_Frank_feedbackInstanceStruct;

#endif                                 /*typedef_SFc2_Frank_feedbackInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c2_Frank_feedback_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c2_Frank_feedback_get_check_sum(mxArray *plhs[]);
extern void c2_Frank_feedback_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
