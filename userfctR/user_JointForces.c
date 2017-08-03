//---------------------------
// UCL-CEREM-MBS
//
// @version MBsysLab_s 1.7.a
//
// Creation : 2006
// Last update : 01/10/2008
//---------------------------

#define _USE_MATH_DEFINES  // to use M_PI in Windows
#include "math.h" 

#include "mbs_data.h"
#include "user_model.h"

#include "set_output.h"

double* user_JointForces(MbsData *mbs_data, double tsim)
{

/*-- Begin of user code --*/

/* 
    // example for a spring-damper in joint number 2
    double K = 10;
    double C = 1;
    double L0 = 0.1;
    mbs_data->Qq[2] = -(K*(mbs_data->q[2]-L0) + C*mbs_data->qd[2]);
*/
    

/*-- End of user code --*/

    return mbs_data->Qq;
}
