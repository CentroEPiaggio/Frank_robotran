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
#include "user_all_id.h"
#include "mbs_data.h"


void user_cons_hJ(double *h, double **Jac, MbsData *mbs_data, double tsim)
{

/*-- Begin of user code --*/

    // declare and define Radius of wheels and width of base
     double R, W;
     R = 0.1;
     W = 0.4570;

    // define the value of the constraint
     
    // x constraint
    h[1] = mbs_data->q[Joint_0_id] - (R/2)*(mbs_data->q[Joint_13_id] + mbs_data->q[Joint_20_id])*cos(mbs_data->q[Joint_3_id]);      
    Jac[1][Joint_0_id] = 1;
    Jac[1][Joint_13_id] = -(R/2)*cos(mbs_data->q[Joint_3_id]);
    Jac[1][Joint_20_id] = -(R/2)*cos(mbs_data->q[Joint_3_id]);
    Jac[1][Joint_3_id] = (R/2)*(mbs_data->q[Joint_13_id] + mbs_data->q[Joint_20_id])*sin(mbs_data->q[Joint_3_id]);
    
    // y constraint
    h[2] = mbs_data->q[Joint_1_id] + (R/2)*(mbs_data->q[Joint_13_id] + mbs_data->q[Joint_20_id])*sin(mbs_data->q[Joint_3_id]);
    Jac[2][Joint_1_id] = 1;
    Jac[2][Joint_13_id] = (R/2)*sin(mbs_data->q[Joint_3_id]);
    Jac[2][Joint_20_id] = (R/2)*sin(mbs_data->q[Joint_3_id]);
    Jac[2][Joint_3_id] = (R/2)*(mbs_data->q[Joint_13_id] + mbs_data->q[Joint_20_id])*cos(mbs_data->q[Joint_3_id]);
    
    // z constraint
    h[3] = mbs_data->q[Joint_2_id] - 0.1;
    Jac[3][Joint_2_id] = 1;
    
    // yaw constraint
    h[4] = mbs_data->q[Joint_3_id] - (R/W)*(mbs_data->q[Joint_13_id] - mbs_data->q[Joint_20_id]);
    Jac[4][Joint_3_id] = 1;
    Jac[4][Joint_13_id] = -(R/W);
    Jac[4][Joint_20_id] = (R/W);
    
    // roll constraint
    h[5] = mbs_data->q[Joint_5_id];
    Jac[5][Joint_5_id] = 1;
    
    
    //...
/* 
    h[1] = ...;      
    Jac[1][1] = ...;
    Jac[1][2] = ...;
    ...
 */
    
/*
    h[2] = 
    Jac[2][1] =
    Jac[2][2] =
    ...
*/

/*-- End of user code --*/

}
