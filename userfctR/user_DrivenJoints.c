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
#include "user_all_id.h"

void user_DrivenJoints(MbsData *mbs_data,double tsim)
{

   // joint id
   int id;
   
//    // Get the Z-prismatic Joint id
//    id = Joint_2_id;
//    
//    // impose the position, velocity and acceleration
//    mbs_data->q[id] = 0.1; // radius of wheel
//    mbs_data->qd[id] = 0;
//    mbs_data->qdd[id] = 0;
   
//    id = Joint_5_id;
//    
//    // impose the position, velocity and acceleration
//    mbs_data->q[id] = 0; 
//    mbs_data->qd[id] = 0;
//    mbs_data->qdd[id] = 0;
   
   // Get the torso Joint id
   id = Joint_6_id;
   
   // impose the position, velocity and acceleration
   mbs_data->q[id] = 0; 
   mbs_data->qd[id] = 0;
   mbs_data->qdd[id] = 0;
   
   // Get the left Joint id
   id = Joint_17_id;
   
   // impose the position, velocity and acceleration
   mbs_data->q[id] = 0; 
   mbs_data->qd[id] = 0;
   mbs_data->qdd[id] = 0;
   
   // Get the torso Joint id
   id = Joint_18_id;
   
   // impose the position, velocity and acceleration
   mbs_data->q[id] = 0; 
   mbs_data->qd[id] = 0;
   mbs_data->qdd[id] = 0;
   // Get the torso Joint id
   id = Joint_14_id;
   
   // impose the position, velocity and acceleration
   mbs_data->q[id] = 0; 
   mbs_data->qd[id] = 0;
   mbs_data->qdd[id] = 0;
   
   
   // Get the right Joint id
   id = Joint_9_id;
   
   // impose the position, velocity and acceleration
   mbs_data->q[id] = 0; 
   mbs_data->qd[id] = 0;
   mbs_data->qdd[id] = 0;
   
   // Get the torso Joint id
   id = Joint_10_id;
   
   // impose the position, velocity and acceleration
   mbs_data->q[id] = 0; 
   mbs_data->qd[id] = 0;
   mbs_data->qdd[id] = 0;
   
   // Get the torso Joint id
   id = Joint_12_id;
   
   // impose the position, velocity and acceleration
   mbs_data->q[id] = 0; 
   mbs_data->qd[id] = 0;
   mbs_data->qdd[id] = 0;
}

 
