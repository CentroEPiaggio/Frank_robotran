//---------------------------
// UCL-CEREM-MBS
//
// @version MBsysLab_s 1.7.a
//
// Creation : 2006
// Last update : 01/10/2008
//
//
//---------------------------

#define _USE_MATH_DEFINES  // to use M_PI in Windows
#include "math.h"

#include "MBSdef.h"
#include "mbs_data.h"
#include "mbs_project_interface.h"

#include "user_all_id.h"

#include "user_IO.h"

// my function
void transpose(void *dest, void *src, int src_h, int src_w)
{
	int i, j;
	double (*d)[src_h] = dest, (*s)[src_w] = src;
	for (i = 0; i < src_h; i++)
		for (j = 0; j < src_w; j++)
			d[j][i] = s[i][j];
}


//
double* user_ExtForces(double PxF[4], double RxF[4][4], 
                       double VxF[4], double OMxF[4], 
                       double AxF[4], double OMPxF[4], 
                       MbsData *mbs_data, double tsim,int ixF)
{
    double Fx=0.0, Fy=0.0, Fz=0.0;
    double Mx=0.0, My=0.0, Mz=0.0;
    double dxF[4] ={0.0, 0.0, 0.0, 0.0};


    double *SWr = mbs_data->SWr[ixF];

    // default application point of the force: anchor point to which it is attached
    int idpt = 0;
    idpt = mbs_data->xfidpt[ixF];
    dxF[1] = mbs_data->dpt[1][idpt];
    dxF[2] = mbs_data->dpt[2][idpt];
    dxF[3] = mbs_data->dpt[3][idpt];

    

    /* Begin of user declaration */
    int i;
    double motorTorque_left, motorTorque_right, Fr[4], Fl[4], Mr[4], Ml[4], Fa, Fas, Mtot, mu; // motors torque 
    double m_wheel, radius_wheel, Inertia_wheel;
    /* End of user declaration */
    
    mu = 0.7;
    motorTorque_left    = mbs_data->user_IO->motorTorque_left;
    motorTorque_right   = mbs_data->user_IO->motorTorque_right;
    
    Mtot            = 8.58;       // Kg;
    m_wheel         = 0.38;       // kg;
    radius_wheel    = 0.1;        // m;
    Inertia_wheel   = 0.0024;     //kgm^2;
    
    switch(ixF){

/* Begin of user code */
        case ExtForce_left_wheel_id:
            
            Fa = motorTorque_left/(radius_wheel*(1 + Inertia_wheel/(m_wheel*radius_wheel*radius_wheel)));
            
            for(i = 1; i<4; i++)
            {
                Fl[i] = RxF[i][1]*Fa;
                Ml[i] = RxF[i][2]*(motorTorque_left - radius_wheel*Fa);
            }
            
            Fx = Fl[1];
            Fy = Fl[2];
            Fz = Fl[3];
            
            Mx = Ml[1];
            My = Ml[2];
            Mz = Ml[3];
            break;

        case ExtForce_right_wheel_id:
            
            Fa = motorTorque_right/(radius_wheel*(1 + Inertia_wheel/(m_wheel*radius_wheel*radius_wheel)));
            
            for(i = 1; i<4; i++)
            {
                Fr[i] = RxF[i][1]*Fa;
                Mr[i] = RxF[i][2]*(motorTorque_left - radius_wheel*Fa);
            }
            
            Fx = Fr[1];
            Fy = Fr[2];
            Fz = Fr[3];
            
            Mx = Mr[1];
            My = Mr[2];
            Mz = Mr[3];
            
            
            break;


/* End of user code */

    }

    SWr[1]=Fx;
    SWr[2]=Fy;
    SWr[3]=Fz;
    SWr[4]=Mx;
    SWr[5]=My;
    SWr[6]=Mz;
    SWr[7]=dxF[1];
    SWr[8]=dxF[2];
    SWr[9]=dxF[3];

    return SWr;
}

 
