//
//-------------------------------------------------------------
//
//	ROBOTRAN - Version 6.6 (build : february 22, 2008)
//
//	Copyright 
//	Universite catholique de Louvain 
//	Departement de Mecanique 
//	Unite de Production Mecanique et Machines 
//	2, Place du Levant 
//	1348 Louvain-la-Neuve 
//	http://www.robotran.be// 
//
//	==> Generation Date : Thu Jul 20 16:41:31 2017
//
//	==> Project name : Frank_segway
//	==> using XML input file 
//
//	==> Number of joints : 17
//
//	==> Function : F 6 : Sensors Kinematical Informations (sens) 
//	==> Flops complexity : 3001
//
//	==> Generation Time :  0.040 seconds
//	==> Post-Processing :  0.060 seconds
//
//-------------------------------------------------------------
//
 
#include <math.h> 

#include "mbs_data.h"
#include "mbs_project_interface.h"
#include "mbs_sensor.h"
 
void  mbs_sensor(MbsSensor *sens, 
              MbsData *s,
              int isens)
{ 
 
#include "mbs_sensor_Frank_segway.h" 
#define q s->q 
#define qd s->qd 
#define qdd s->qdd 
 
 

// === begin imp_aux === 

// === end imp_aux === 

// ===== BEGIN task 0 ===== 
 
// Sensor Kinematics 



// = = Block_0_0_0_0_0_1 = = 
 
// Trigonometric Variables  

  C4 = cos(q[4]);
  S4 = sin(q[4]);
  C5 = cos(q[5]);
  S5 = sin(q[5]);
  C6 = cos(q[6]);
  S6 = sin(q[6]);

// = = Block_0_0_0_0_0_2 = = 
 
// Augmented Joint Position Vectors   

  Dz73 = q[7]+s->dpt[3][2];

// = = Block_0_0_0_0_0_3 = = 
 
// Trigonometric Variables  

  C8 = cos(q[8]);
  S8 = sin(q[8]);
  C9 = cos(q[9]);
  S9 = sin(q[9]);
  C10 = cos(q[10]);
  S10 = sin(q[10]);

// = = Block_0_0_0_0_0_4 = = 
 
// Trigonometric Variables  

  C11 = cos(q[11]);
  S11 = sin(q[11]);
  C12 = cos(q[12]);
  S12 = sin(q[12]);
  C13 = cos(q[13]);
  S13 = sin(q[13]);

// = = Block_0_0_0_0_0_5 = = 
 
// Trigonometric Variables  

  C14 = cos(q[14]);
  S14 = sin(q[14]);

// = = Block_0_0_0_0_0_6 = = 
 
// Trigonometric Variables  

  C16 = cos(q[16]);
  S16 = sin(q[16]);

// = = Block_0_0_0_0_0_7 = = 
 
// Trigonometric Variables  

  C17 = cos(q[17]);
  S17 = sin(q[17]);

// ====== END Task 0 ====== 

// ===== BEGIN task 1 ===== 
 
switch(isens)
{
 
// 
break;
case 1:
 


// = = Block_1_0_0_1_0_1 = = 
 
// Sensor Kinematics 


    ROcp0_15 = C4*C5;
    ROcp0_25 = S4*C5;
    ROcp0_75 = C4*S5;
    ROcp0_85 = S4*S5;
    ROcp0_46 = ROcp0_75*S6-S4*C6;
    ROcp0_56 = ROcp0_85*S6+C4*C6;
    ROcp0_66 = C5*S6;
    ROcp0_76 = ROcp0_75*C6+S4*S6;
    ROcp0_86 = ROcp0_85*C6-C4*S6;
    ROcp0_96 = C5*C6;
    OMcp0_15 = -qd[5]*S4;
    OMcp0_25 = qd[5]*C4;
    OMcp0_16 = OMcp0_15+ROcp0_15*qd[6];
    OMcp0_26 = OMcp0_25+ROcp0_25*qd[6];
    OMcp0_36 = qd[4]-qd[6]*S5;
    OPcp0_16 = ROcp0_15*qdd[6]-qdd[5]*S4-qd[4]*qd[5]*C4-qd[6]*(OMcp0_25*S5+ROcp0_25*qd[4]);
    OPcp0_26 = ROcp0_25*qdd[6]+qdd[5]*C4-qd[4]*qd[5]*S4+qd[6]*(OMcp0_15*S5+ROcp0_15*qd[4]);
    OPcp0_36 = qdd[4]-qdd[6]*S5-qd[5]*qd[6]*C5;

// = = Block_1_0_0_1_1_0 = = 
 
// Symbolic Outputs  

    sens->P[1] = q[1];
    sens->P[2] = q[2];
    sens->P[3] = q[3];
    sens->R[1][1] = ROcp0_15;
    sens->R[1][2] = ROcp0_25;
    sens->R[1][3] = -S5;
    sens->R[2][1] = ROcp0_46;
    sens->R[2][2] = ROcp0_56;
    sens->R[2][3] = ROcp0_66;
    sens->R[3][1] = ROcp0_76;
    sens->R[3][2] = ROcp0_86;
    sens->R[3][3] = ROcp0_96;
    sens->V[1] = qd[1];
    sens->V[2] = qd[2];
    sens->V[3] = qd[3];
    sens->OM[1] = OMcp0_16;
    sens->OM[2] = OMcp0_26;
    sens->OM[3] = OMcp0_36;
    sens->J[1][1] = (1.0);
    sens->J[2][2] = (1.0);
    sens->J[3][3] = (1.0);
    sens->J[4][5] = -S4;
    sens->J[4][6] = ROcp0_15;
    sens->J[5][5] = C4;
    sens->J[5][6] = ROcp0_25;
    sens->J[6][4] = (1.0);
    sens->J[6][6] = -S5;
    sens->A[1] = qdd[1];
    sens->A[2] = qdd[2];
    sens->A[3] = qdd[3];
    sens->OMP[1] = OPcp0_16;
    sens->OMP[2] = OPcp0_26;
    sens->OMP[3] = OPcp0_36;
 
// 
break;
case 2:
 


// = = Block_1_0_0_2_0_1 = = 
 
// Sensor Kinematics 


    ROcp1_15 = C4*C5;
    ROcp1_25 = S4*C5;
    ROcp1_75 = C4*S5;
    ROcp1_85 = S4*S5;
    ROcp1_46 = ROcp1_75*S6-S4*C6;
    ROcp1_56 = ROcp1_85*S6+C4*C6;
    ROcp1_66 = C5*S6;
    ROcp1_76 = ROcp1_75*C6+S4*S6;
    ROcp1_86 = ROcp1_85*C6-C4*S6;
    ROcp1_96 = C5*C6;
    OMcp1_15 = -qd[5]*S4;
    OMcp1_25 = qd[5]*C4;
    OMcp1_16 = OMcp1_15+ROcp1_15*qd[6];
    OMcp1_26 = OMcp1_25+ROcp1_25*qd[6];
    OMcp1_36 = qd[4]-qd[6]*S5;
    OPcp1_16 = ROcp1_15*qdd[6]-qdd[5]*S4-qd[4]*qd[5]*C4-qd[6]*(OMcp1_25*S5+ROcp1_25*qd[4]);
    OPcp1_26 = ROcp1_25*qdd[6]+qdd[5]*C4-qd[4]*qd[5]*S4+qd[6]*(OMcp1_15*S5+ROcp1_15*qd[4]);
    OPcp1_36 = qdd[4]-qdd[6]*S5-qd[5]*qd[6]*C5;

// = = Block_1_0_0_2_0_2 = = 
 
// Sensor Kinematics 


    RLcp1_17 = Dz73*ROcp1_76;
    RLcp1_27 = Dz73*ROcp1_86;
    RLcp1_37 = Dz73*ROcp1_96;
    POcp1_17 = RLcp1_17+q[1];
    POcp1_27 = RLcp1_27+q[2];
    POcp1_37 = RLcp1_37+q[3];
    JTcp1_17_5 = RLcp1_37*C4;
    JTcp1_27_5 = RLcp1_37*S4;
    JTcp1_37_5 = -(RLcp1_17*C4+RLcp1_27*S4);
    JTcp1_17_6 = RLcp1_27*S5+RLcp1_37*ROcp1_25;
    JTcp1_27_6 = -(RLcp1_17*S5+RLcp1_37*ROcp1_15);
    JTcp1_37_6 = -(RLcp1_17*ROcp1_25-RLcp1_27*ROcp1_15);
    ORcp1_17 = OMcp1_26*RLcp1_37-OMcp1_36*RLcp1_27;
    ORcp1_27 = -(OMcp1_16*RLcp1_37-OMcp1_36*RLcp1_17);
    ORcp1_37 = OMcp1_16*RLcp1_27-OMcp1_26*RLcp1_17;
    VIcp1_17 = ORcp1_17+qd[1]+ROcp1_76*qd[7];
    VIcp1_27 = ORcp1_27+qd[2]+ROcp1_86*qd[7];
    VIcp1_37 = ORcp1_37+qd[3]+ROcp1_96*qd[7];
    ACcp1_17 = qdd[1]+OMcp1_26*ORcp1_37-OMcp1_36*ORcp1_27+OPcp1_26*RLcp1_37-OPcp1_36*RLcp1_27+ROcp1_76*qdd[7]+(2.0)*qd[7]*(
 OMcp1_26*ROcp1_96-OMcp1_36*ROcp1_86);
    ACcp1_27 = qdd[2]-OMcp1_16*ORcp1_37+OMcp1_36*ORcp1_17-OPcp1_16*RLcp1_37+OPcp1_36*RLcp1_17+ROcp1_86*qdd[7]-(2.0)*qd[7]*(
 OMcp1_16*ROcp1_96-OMcp1_36*ROcp1_76);
    ACcp1_37 = qdd[3]+OMcp1_16*ORcp1_27-OMcp1_26*ORcp1_17+OPcp1_16*RLcp1_27-OPcp1_26*RLcp1_17+ROcp1_96*qdd[7]+(2.0)*qd[7]*(
 OMcp1_16*ROcp1_86-OMcp1_26*ROcp1_76);

// = = Block_1_0_0_2_1_0 = = 
 
// Symbolic Outputs  

    sens->P[1] = POcp1_17;
    sens->P[2] = POcp1_27;
    sens->P[3] = POcp1_37;
    sens->R[1][1] = ROcp1_15;
    sens->R[1][2] = ROcp1_25;
    sens->R[1][3] = -S5;
    sens->R[2][1] = ROcp1_46;
    sens->R[2][2] = ROcp1_56;
    sens->R[2][3] = ROcp1_66;
    sens->R[3][1] = ROcp1_76;
    sens->R[3][2] = ROcp1_86;
    sens->R[3][3] = ROcp1_96;
    sens->V[1] = VIcp1_17;
    sens->V[2] = VIcp1_27;
    sens->V[3] = VIcp1_37;
    sens->OM[1] = OMcp1_16;
    sens->OM[2] = OMcp1_26;
    sens->OM[3] = OMcp1_36;
    sens->J[1][1] = (1.0);
    sens->J[1][4] = -RLcp1_27;
    sens->J[1][5] = JTcp1_17_5;
    sens->J[1][6] = JTcp1_17_6;
    sens->J[1][7] = ROcp1_76;
    sens->J[2][2] = (1.0);
    sens->J[2][4] = RLcp1_17;
    sens->J[2][5] = JTcp1_27_5;
    sens->J[2][6] = JTcp1_27_6;
    sens->J[2][7] = ROcp1_86;
    sens->J[3][3] = (1.0);
    sens->J[3][5] = JTcp1_37_5;
    sens->J[3][6] = JTcp1_37_6;
    sens->J[3][7] = ROcp1_96;
    sens->J[4][5] = -S4;
    sens->J[4][6] = ROcp1_15;
    sens->J[5][5] = C4;
    sens->J[5][6] = ROcp1_25;
    sens->J[6][4] = (1.0);
    sens->J[6][6] = -S5;
    sens->A[1] = ACcp1_17;
    sens->A[2] = ACcp1_27;
    sens->A[3] = ACcp1_37;
    sens->OMP[1] = OPcp1_16;
    sens->OMP[2] = OPcp1_26;
    sens->OMP[3] = OPcp1_36;
 
// 
break;
case 3:
 


// = = Block_1_0_0_3_0_1 = = 
 
// Sensor Kinematics 


    ROcp2_15 = C4*C5;
    ROcp2_25 = S4*C5;
    ROcp2_75 = C4*S5;
    ROcp2_85 = S4*S5;
    ROcp2_46 = ROcp2_75*S6-S4*C6;
    ROcp2_56 = ROcp2_85*S6+C4*C6;
    ROcp2_66 = C5*S6;
    ROcp2_76 = ROcp2_75*C6+S4*S6;
    ROcp2_86 = ROcp2_85*C6-C4*S6;
    ROcp2_96 = C5*C6;
    OMcp2_15 = -qd[5]*S4;
    OMcp2_25 = qd[5]*C4;
    OMcp2_16 = OMcp2_15+ROcp2_15*qd[6];
    OMcp2_26 = OMcp2_25+ROcp2_25*qd[6];
    OMcp2_36 = qd[4]-qd[6]*S5;
    OPcp2_16 = ROcp2_15*qdd[6]-qdd[5]*S4-qd[4]*qd[5]*C4-qd[6]*(OMcp2_25*S5+ROcp2_25*qd[4]);
    OPcp2_26 = ROcp2_25*qdd[6]+qdd[5]*C4-qd[4]*qd[5]*S4+qd[6]*(OMcp2_15*S5+ROcp2_15*qd[4]);
    OPcp2_36 = qdd[4]-qdd[6]*S5-qd[5]*qd[6]*C5;

// = = Block_1_0_0_3_0_2 = = 
 
// Sensor Kinematics 


    RLcp2_17 = Dz73*ROcp2_76;
    RLcp2_27 = Dz73*ROcp2_86;
    RLcp2_37 = Dz73*ROcp2_96;
    ORcp2_17 = OMcp2_26*RLcp2_37-OMcp2_36*RLcp2_27;
    ORcp2_27 = -(OMcp2_16*RLcp2_37-OMcp2_36*RLcp2_17);
    ORcp2_37 = OMcp2_16*RLcp2_27-OMcp2_26*RLcp2_17;

// = = Block_1_0_0_3_0_3 = = 
 
// Sensor Kinematics 


    ROcp2_18 = ROcp2_15*C8-ROcp2_76*S8;
    ROcp2_28 = ROcp2_25*C8-ROcp2_86*S8;
    ROcp2_38 = -(ROcp2_96*S8+S5*C8);
    ROcp2_78 = ROcp2_15*S8+ROcp2_76*C8;
    ROcp2_88 = ROcp2_25*S8+ROcp2_86*C8;
    ROcp2_98 = ROcp2_96*C8-S5*S8;
    RLcp2_18 = ROcp2_46*s->dpt[2][6];
    RLcp2_28 = ROcp2_56*s->dpt[2][6];
    RLcp2_38 = ROcp2_66*s->dpt[2][6];
    POcp2_18 = RLcp2_17+RLcp2_18+q[1];
    POcp2_28 = RLcp2_27+RLcp2_28+q[2];
    POcp2_38 = RLcp2_37+RLcp2_38+q[3];
    JTcp2_18_4 = -(RLcp2_27+RLcp2_28);
    JTcp2_28_4 = RLcp2_17+RLcp2_18;
    JTcp2_18_5 = C4*(RLcp2_37+RLcp2_38);
    JTcp2_28_5 = S4*(RLcp2_37+RLcp2_38);
    JTcp2_38_5 = -(C4*(RLcp2_17+RLcp2_18)+S4*(RLcp2_27+RLcp2_28));
    JTcp2_18_6 = ROcp2_25*(RLcp2_37+RLcp2_38)+S5*(RLcp2_27+RLcp2_28);
    JTcp2_28_6 = -(ROcp2_15*(RLcp2_37+RLcp2_38)+S5*(RLcp2_17+RLcp2_18));
    JTcp2_38_6 = ROcp2_15*(RLcp2_27+RLcp2_28)-ROcp2_25*(RLcp2_17+RLcp2_18);
    OMcp2_18 = OMcp2_16+ROcp2_46*qd[8];
    OMcp2_28 = OMcp2_26+ROcp2_56*qd[8];
    OMcp2_38 = OMcp2_36+ROcp2_66*qd[8];
    ORcp2_18 = OMcp2_26*RLcp2_38-OMcp2_36*RLcp2_28;
    ORcp2_28 = -(OMcp2_16*RLcp2_38-OMcp2_36*RLcp2_18);
    ORcp2_38 = OMcp2_16*RLcp2_28-OMcp2_26*RLcp2_18;
    VIcp2_18 = ORcp2_17+ORcp2_18+qd[1]+ROcp2_76*qd[7];
    VIcp2_28 = ORcp2_27+ORcp2_28+qd[2]+ROcp2_86*qd[7];
    VIcp2_38 = ORcp2_37+ORcp2_38+qd[3]+ROcp2_96*qd[7];
    OPcp2_18 = OPcp2_16+ROcp2_46*qdd[8]+qd[8]*(OMcp2_26*ROcp2_66-OMcp2_36*ROcp2_56);
    OPcp2_28 = OPcp2_26+ROcp2_56*qdd[8]-qd[8]*(OMcp2_16*ROcp2_66-OMcp2_36*ROcp2_46);
    OPcp2_38 = OPcp2_36+ROcp2_66*qdd[8]+qd[8]*(OMcp2_16*ROcp2_56-OMcp2_26*ROcp2_46);
    ACcp2_18 = qdd[1]+OMcp2_26*ORcp2_37+OMcp2_26*ORcp2_38-OMcp2_36*ORcp2_27-OMcp2_36*ORcp2_28+OPcp2_26*RLcp2_37+OPcp2_26*
 RLcp2_38-OPcp2_36*RLcp2_27-OPcp2_36*RLcp2_28+ROcp2_76*qdd[7]+(2.0)*qd[7]*(OMcp2_26*ROcp2_96-OMcp2_36*ROcp2_86);
    ACcp2_28 = qdd[2]-OMcp2_16*ORcp2_37-OMcp2_16*ORcp2_38+OMcp2_36*ORcp2_17+OMcp2_36*ORcp2_18-OPcp2_16*RLcp2_37-OPcp2_16*
 RLcp2_38+OPcp2_36*RLcp2_17+OPcp2_36*RLcp2_18+ROcp2_86*qdd[7]-(2.0)*qd[7]*(OMcp2_16*ROcp2_96-OMcp2_36*ROcp2_76);
    ACcp2_38 = qdd[3]+OMcp2_16*ORcp2_27+OMcp2_16*ORcp2_28-OMcp2_26*ORcp2_17-OMcp2_26*ORcp2_18+OPcp2_16*RLcp2_27+OPcp2_16*
 RLcp2_28-OPcp2_26*RLcp2_17-OPcp2_26*RLcp2_18+ROcp2_96*qdd[7]+(2.0)*qd[7]*(OMcp2_16*ROcp2_86-OMcp2_26*ROcp2_76);

// = = Block_1_0_0_3_1_0 = = 
 
// Symbolic Outputs  

    sens->P[1] = POcp2_18;
    sens->P[2] = POcp2_28;
    sens->P[3] = POcp2_38;
    sens->R[1][1] = ROcp2_18;
    sens->R[1][2] = ROcp2_28;
    sens->R[1][3] = ROcp2_38;
    sens->R[2][1] = ROcp2_46;
    sens->R[2][2] = ROcp2_56;
    sens->R[2][3] = ROcp2_66;
    sens->R[3][1] = ROcp2_78;
    sens->R[3][2] = ROcp2_88;
    sens->R[3][3] = ROcp2_98;
    sens->V[1] = VIcp2_18;
    sens->V[2] = VIcp2_28;
    sens->V[3] = VIcp2_38;
    sens->OM[1] = OMcp2_18;
    sens->OM[2] = OMcp2_28;
    sens->OM[3] = OMcp2_38;
    sens->J[1][1] = (1.0);
    sens->J[1][4] = JTcp2_18_4;
    sens->J[1][5] = JTcp2_18_5;
    sens->J[1][6] = JTcp2_18_6;
    sens->J[1][7] = ROcp2_76;
    sens->J[2][2] = (1.0);
    sens->J[2][4] = JTcp2_28_4;
    sens->J[2][5] = JTcp2_28_5;
    sens->J[2][6] = JTcp2_28_6;
    sens->J[2][7] = ROcp2_86;
    sens->J[3][3] = (1.0);
    sens->J[3][5] = JTcp2_38_5;
    sens->J[3][6] = JTcp2_38_6;
    sens->J[3][7] = ROcp2_96;
    sens->J[4][5] = -S4;
    sens->J[4][6] = ROcp2_15;
    sens->J[4][8] = ROcp2_46;
    sens->J[5][5] = C4;
    sens->J[5][6] = ROcp2_25;
    sens->J[5][8] = ROcp2_56;
    sens->J[6][4] = (1.0);
    sens->J[6][6] = -S5;
    sens->J[6][8] = ROcp2_66;
    sens->A[1] = ACcp2_18;
    sens->A[2] = ACcp2_28;
    sens->A[3] = ACcp2_38;
    sens->OMP[1] = OPcp2_18;
    sens->OMP[2] = OPcp2_28;
    sens->OMP[3] = OPcp2_38;
 
// 
break;
case 4:
 


// = = Block_1_0_0_4_0_1 = = 
 
// Sensor Kinematics 


    ROcp3_15 = C4*C5;
    ROcp3_25 = S4*C5;
    ROcp3_75 = C4*S5;
    ROcp3_85 = S4*S5;
    ROcp3_46 = ROcp3_75*S6-S4*C6;
    ROcp3_56 = ROcp3_85*S6+C4*C6;
    ROcp3_66 = C5*S6;
    ROcp3_76 = ROcp3_75*C6+S4*S6;
    ROcp3_86 = ROcp3_85*C6-C4*S6;
    ROcp3_96 = C5*C6;
    OMcp3_15 = -qd[5]*S4;
    OMcp3_25 = qd[5]*C4;
    OMcp3_16 = OMcp3_15+ROcp3_15*qd[6];
    OMcp3_26 = OMcp3_25+ROcp3_25*qd[6];
    OMcp3_36 = qd[4]-qd[6]*S5;
    OPcp3_16 = ROcp3_15*qdd[6]-qdd[5]*S4-qd[4]*qd[5]*C4-qd[6]*(OMcp3_25*S5+ROcp3_25*qd[4]);
    OPcp3_26 = ROcp3_25*qdd[6]+qdd[5]*C4-qd[4]*qd[5]*S4+qd[6]*(OMcp3_15*S5+ROcp3_15*qd[4]);
    OPcp3_36 = qdd[4]-qdd[6]*S5-qd[5]*qd[6]*C5;

// = = Block_1_0_0_4_0_2 = = 
 
// Sensor Kinematics 


    RLcp3_17 = Dz73*ROcp3_76;
    RLcp3_27 = Dz73*ROcp3_86;
    RLcp3_37 = Dz73*ROcp3_96;
    ORcp3_17 = OMcp3_26*RLcp3_37-OMcp3_36*RLcp3_27;
    ORcp3_27 = -(OMcp3_16*RLcp3_37-OMcp3_36*RLcp3_17);
    ORcp3_37 = OMcp3_16*RLcp3_27-OMcp3_26*RLcp3_17;

// = = Block_1_0_0_4_0_3 = = 
 
// Sensor Kinematics 


    ROcp3_18 = ROcp3_15*C8-ROcp3_76*S8;
    ROcp3_28 = ROcp3_25*C8-ROcp3_86*S8;
    ROcp3_38 = -(ROcp3_96*S8+S5*C8);
    ROcp3_78 = ROcp3_15*S8+ROcp3_76*C8;
    ROcp3_88 = ROcp3_25*S8+ROcp3_86*C8;
    ROcp3_98 = ROcp3_96*C8-S5*S8;
    ROcp3_49 = ROcp3_46*C9+ROcp3_78*S9;
    ROcp3_59 = ROcp3_56*C9+ROcp3_88*S9;
    ROcp3_69 = ROcp3_66*C9+ROcp3_98*S9;
    ROcp3_79 = -(ROcp3_46*S9-ROcp3_78*C9);
    ROcp3_89 = -(ROcp3_56*S9-ROcp3_88*C9);
    ROcp3_99 = -(ROcp3_66*S9-ROcp3_98*C9);
    RLcp3_18 = ROcp3_46*s->dpt[2][6];
    RLcp3_28 = ROcp3_56*s->dpt[2][6];
    RLcp3_38 = ROcp3_66*s->dpt[2][6];
    POcp3_18 = RLcp3_17+RLcp3_18+q[1];
    POcp3_28 = RLcp3_27+RLcp3_28+q[2];
    POcp3_38 = RLcp3_37+RLcp3_38+q[3];
    JTcp3_18_4 = -(RLcp3_27+RLcp3_28);
    JTcp3_28_4 = RLcp3_17+RLcp3_18;
    JTcp3_18_5 = C4*(RLcp3_37+RLcp3_38);
    JTcp3_28_5 = S4*(RLcp3_37+RLcp3_38);
    JTcp3_38_5 = -(C4*(RLcp3_17+RLcp3_18)+S4*(RLcp3_27+RLcp3_28));
    JTcp3_18_6 = ROcp3_25*(RLcp3_37+RLcp3_38)+S5*(RLcp3_27+RLcp3_28);
    JTcp3_28_6 = -(ROcp3_15*(RLcp3_37+RLcp3_38)+S5*(RLcp3_17+RLcp3_18));
    JTcp3_38_6 = ROcp3_15*(RLcp3_27+RLcp3_28)-ROcp3_25*(RLcp3_17+RLcp3_18);
    OMcp3_18 = OMcp3_16+ROcp3_46*qd[8];
    OMcp3_28 = OMcp3_26+ROcp3_56*qd[8];
    OMcp3_38 = OMcp3_36+ROcp3_66*qd[8];
    ORcp3_18 = OMcp3_26*RLcp3_38-OMcp3_36*RLcp3_28;
    ORcp3_28 = -(OMcp3_16*RLcp3_38-OMcp3_36*RLcp3_18);
    ORcp3_38 = OMcp3_16*RLcp3_28-OMcp3_26*RLcp3_18;
    VIcp3_18 = ORcp3_17+ORcp3_18+qd[1]+ROcp3_76*qd[7];
    VIcp3_28 = ORcp3_27+ORcp3_28+qd[2]+ROcp3_86*qd[7];
    VIcp3_38 = ORcp3_37+ORcp3_38+qd[3]+ROcp3_96*qd[7];
    ACcp3_18 = qdd[1]+OMcp3_26*ORcp3_37+OMcp3_26*ORcp3_38-OMcp3_36*ORcp3_27-OMcp3_36*ORcp3_28+OPcp3_26*RLcp3_37+OPcp3_26*
 RLcp3_38-OPcp3_36*RLcp3_27-OPcp3_36*RLcp3_28+ROcp3_76*qdd[7]+(2.0)*qd[7]*(OMcp3_26*ROcp3_96-OMcp3_36*ROcp3_86);
    ACcp3_28 = qdd[2]-OMcp3_16*ORcp3_37-OMcp3_16*ORcp3_38+OMcp3_36*ORcp3_17+OMcp3_36*ORcp3_18-OPcp3_16*RLcp3_37-OPcp3_16*
 RLcp3_38+OPcp3_36*RLcp3_17+OPcp3_36*RLcp3_18+ROcp3_86*qdd[7]-(2.0)*qd[7]*(OMcp3_16*ROcp3_96-OMcp3_36*ROcp3_76);
    ACcp3_38 = qdd[3]+OMcp3_16*ORcp3_27+OMcp3_16*ORcp3_28-OMcp3_26*ORcp3_17-OMcp3_26*ORcp3_18+OPcp3_16*RLcp3_27+OPcp3_16*
 RLcp3_28-OPcp3_26*RLcp3_17-OPcp3_26*RLcp3_18+ROcp3_96*qdd[7]+(2.0)*qd[7]*(OMcp3_16*ROcp3_86-OMcp3_26*ROcp3_76);
    OMcp3_19 = OMcp3_18+ROcp3_18*qd[9];
    OMcp3_29 = OMcp3_28+ROcp3_28*qd[9];
    OMcp3_39 = OMcp3_38+ROcp3_38*qd[9];
    OPcp3_19 = OPcp3_16+ROcp3_18*qdd[9]+ROcp3_46*qdd[8]+qd[8]*(OMcp3_26*ROcp3_66-OMcp3_36*ROcp3_56)+qd[9]*(OMcp3_28*
 ROcp3_38-OMcp3_38*ROcp3_28);
    OPcp3_29 = OPcp3_26+ROcp3_28*qdd[9]+ROcp3_56*qdd[8]-qd[8]*(OMcp3_16*ROcp3_66-OMcp3_36*ROcp3_46)-qd[9]*(OMcp3_18*
 ROcp3_38-OMcp3_38*ROcp3_18);
    OPcp3_39 = OPcp3_36+ROcp3_38*qdd[9]+ROcp3_66*qdd[8]+qd[8]*(OMcp3_16*ROcp3_56-OMcp3_26*ROcp3_46)+qd[9]*(OMcp3_18*
 ROcp3_28-OMcp3_28*ROcp3_18);

// = = Block_1_0_0_4_1_0 = = 
 
// Symbolic Outputs  

    sens->P[1] = POcp3_18;
    sens->P[2] = POcp3_28;
    sens->P[3] = POcp3_38;
    sens->R[1][1] = ROcp3_18;
    sens->R[1][2] = ROcp3_28;
    sens->R[1][3] = ROcp3_38;
    sens->R[2][1] = ROcp3_49;
    sens->R[2][2] = ROcp3_59;
    sens->R[2][3] = ROcp3_69;
    sens->R[3][1] = ROcp3_79;
    sens->R[3][2] = ROcp3_89;
    sens->R[3][3] = ROcp3_99;
    sens->V[1] = VIcp3_18;
    sens->V[2] = VIcp3_28;
    sens->V[3] = VIcp3_38;
    sens->OM[1] = OMcp3_19;
    sens->OM[2] = OMcp3_29;
    sens->OM[3] = OMcp3_39;
    sens->J[1][1] = (1.0);
    sens->J[1][4] = JTcp3_18_4;
    sens->J[1][5] = JTcp3_18_5;
    sens->J[1][6] = JTcp3_18_6;
    sens->J[1][7] = ROcp3_76;
    sens->J[2][2] = (1.0);
    sens->J[2][4] = JTcp3_28_4;
    sens->J[2][5] = JTcp3_28_5;
    sens->J[2][6] = JTcp3_28_6;
    sens->J[2][7] = ROcp3_86;
    sens->J[3][3] = (1.0);
    sens->J[3][5] = JTcp3_38_5;
    sens->J[3][6] = JTcp3_38_6;
    sens->J[3][7] = ROcp3_96;
    sens->J[4][5] = -S4;
    sens->J[4][6] = ROcp3_15;
    sens->J[4][8] = ROcp3_46;
    sens->J[4][9] = ROcp3_18;
    sens->J[5][5] = C4;
    sens->J[5][6] = ROcp3_25;
    sens->J[5][8] = ROcp3_56;
    sens->J[5][9] = ROcp3_28;
    sens->J[6][4] = (1.0);
    sens->J[6][6] = -S5;
    sens->J[6][8] = ROcp3_66;
    sens->J[6][9] = ROcp3_38;
    sens->A[1] = ACcp3_18;
    sens->A[2] = ACcp3_28;
    sens->A[3] = ACcp3_38;
    sens->OMP[1] = OPcp3_19;
    sens->OMP[2] = OPcp3_29;
    sens->OMP[3] = OPcp3_39;
 
// 
break;
case 5:
 


// = = Block_1_0_0_5_0_1 = = 
 
// Sensor Kinematics 


    ROcp4_15 = C4*C5;
    ROcp4_25 = S4*C5;
    ROcp4_75 = C4*S5;
    ROcp4_85 = S4*S5;
    ROcp4_46 = ROcp4_75*S6-S4*C6;
    ROcp4_56 = ROcp4_85*S6+C4*C6;
    ROcp4_66 = C5*S6;
    ROcp4_76 = ROcp4_75*C6+S4*S6;
    ROcp4_86 = ROcp4_85*C6-C4*S6;
    ROcp4_96 = C5*C6;
    OMcp4_15 = -qd[5]*S4;
    OMcp4_25 = qd[5]*C4;
    OMcp4_16 = OMcp4_15+ROcp4_15*qd[6];
    OMcp4_26 = OMcp4_25+ROcp4_25*qd[6];
    OMcp4_36 = qd[4]-qd[6]*S5;
    OPcp4_16 = ROcp4_15*qdd[6]-qdd[5]*S4-qd[4]*qd[5]*C4-qd[6]*(OMcp4_25*S5+ROcp4_25*qd[4]);
    OPcp4_26 = ROcp4_25*qdd[6]+qdd[5]*C4-qd[4]*qd[5]*S4+qd[6]*(OMcp4_15*S5+ROcp4_15*qd[4]);
    OPcp4_36 = qdd[4]-qdd[6]*S5-qd[5]*qd[6]*C5;

// = = Block_1_0_0_5_0_2 = = 
 
// Sensor Kinematics 


    RLcp4_17 = Dz73*ROcp4_76;
    RLcp4_27 = Dz73*ROcp4_86;
    RLcp4_37 = Dz73*ROcp4_96;
    ORcp4_17 = OMcp4_26*RLcp4_37-OMcp4_36*RLcp4_27;
    ORcp4_27 = -(OMcp4_16*RLcp4_37-OMcp4_36*RLcp4_17);
    ORcp4_37 = OMcp4_16*RLcp4_27-OMcp4_26*RLcp4_17;

// = = Block_1_0_0_5_0_3 = = 
 
// Sensor Kinematics 


    ROcp4_18 = ROcp4_15*C8-ROcp4_76*S8;
    ROcp4_28 = ROcp4_25*C8-ROcp4_86*S8;
    ROcp4_38 = -(ROcp4_96*S8+S5*C8);
    ROcp4_78 = ROcp4_15*S8+ROcp4_76*C8;
    ROcp4_88 = ROcp4_25*S8+ROcp4_86*C8;
    ROcp4_98 = ROcp4_96*C8-S5*S8;
    ROcp4_49 = ROcp4_46*C9+ROcp4_78*S9;
    ROcp4_59 = ROcp4_56*C9+ROcp4_88*S9;
    ROcp4_69 = ROcp4_66*C9+ROcp4_98*S9;
    ROcp4_79 = -(ROcp4_46*S9-ROcp4_78*C9);
    ROcp4_89 = -(ROcp4_56*S9-ROcp4_88*C9);
    ROcp4_99 = -(ROcp4_66*S9-ROcp4_98*C9);
    ROcp4_110 = ROcp4_18*C10+ROcp4_49*S10;
    ROcp4_210 = ROcp4_28*C10+ROcp4_59*S10;
    ROcp4_310 = ROcp4_38*C10+ROcp4_69*S10;
    ROcp4_410 = -(ROcp4_18*S10-ROcp4_49*C10);
    ROcp4_510 = -(ROcp4_28*S10-ROcp4_59*C10);
    ROcp4_610 = -(ROcp4_38*S10-ROcp4_69*C10);
    RLcp4_18 = ROcp4_46*s->dpt[2][6];
    RLcp4_28 = ROcp4_56*s->dpt[2][6];
    RLcp4_38 = ROcp4_66*s->dpt[2][6];
    OMcp4_18 = OMcp4_16+ROcp4_46*qd[8];
    OMcp4_28 = OMcp4_26+ROcp4_56*qd[8];
    OMcp4_38 = OMcp4_36+ROcp4_66*qd[8];
    ORcp4_18 = OMcp4_26*RLcp4_38-OMcp4_36*RLcp4_28;
    ORcp4_28 = -(OMcp4_16*RLcp4_38-OMcp4_36*RLcp4_18);
    ORcp4_38 = OMcp4_16*RLcp4_28-OMcp4_26*RLcp4_18;
    OMcp4_19 = OMcp4_18+ROcp4_18*qd[9];
    OMcp4_29 = OMcp4_28+ROcp4_28*qd[9];
    OMcp4_39 = OMcp4_38+ROcp4_38*qd[9];
    OPcp4_19 = OPcp4_16+ROcp4_18*qdd[9]+ROcp4_46*qdd[8]+qd[8]*(OMcp4_26*ROcp4_66-OMcp4_36*ROcp4_56)+qd[9]*(OMcp4_28*
 ROcp4_38-OMcp4_38*ROcp4_28);
    OPcp4_29 = OPcp4_26+ROcp4_28*qdd[9]+ROcp4_56*qdd[8]-qd[8]*(OMcp4_16*ROcp4_66-OMcp4_36*ROcp4_46)-qd[9]*(OMcp4_18*
 ROcp4_38-OMcp4_38*ROcp4_18);
    OPcp4_39 = OPcp4_36+ROcp4_38*qdd[9]+ROcp4_66*qdd[8]+qd[8]*(OMcp4_16*ROcp4_56-OMcp4_26*ROcp4_46)+qd[9]*(OMcp4_18*
 ROcp4_28-OMcp4_28*ROcp4_18);
    RLcp4_110 = ROcp4_79*s->dpt[3][11];
    RLcp4_210 = ROcp4_89*s->dpt[3][11];
    RLcp4_310 = ROcp4_99*s->dpt[3][11];
    POcp4_110 = RLcp4_110+RLcp4_17+RLcp4_18+q[1];
    POcp4_210 = RLcp4_210+RLcp4_27+RLcp4_28+q[2];
    POcp4_310 = RLcp4_310+RLcp4_37+RLcp4_38+q[3];
    JTcp4_110_4 = -(RLcp4_210+RLcp4_27+RLcp4_28);
    JTcp4_210_4 = RLcp4_110+RLcp4_17+RLcp4_18;
    JTcp4_110_5 = C4*(RLcp4_310+RLcp4_37+RLcp4_38);
    JTcp4_210_5 = S4*(RLcp4_310+RLcp4_37+RLcp4_38);
    JTcp4_310_5 = -(RLcp4_210*S4+C4*(RLcp4_110+RLcp4_17+RLcp4_18)+S4*(RLcp4_27+RLcp4_28));
    JTcp4_110_6 = RLcp4_210*S5+RLcp4_310*ROcp4_25+ROcp4_25*(RLcp4_37+RLcp4_38)+S5*(RLcp4_27+RLcp4_28);
    JTcp4_210_6 = -(RLcp4_110*S5+RLcp4_310*ROcp4_15+ROcp4_15*(RLcp4_37+RLcp4_38)+S5*(RLcp4_17+RLcp4_18));
    JTcp4_310_6 = ROcp4_15*(RLcp4_27+RLcp4_28)-ROcp4_25*(RLcp4_17+RLcp4_18)-RLcp4_110*ROcp4_25+RLcp4_210*ROcp4_15;
    JTcp4_110_8 = -(RLcp4_210*ROcp4_66-RLcp4_310*ROcp4_56);
    JTcp4_210_8 = RLcp4_110*ROcp4_66-RLcp4_310*ROcp4_46;
    JTcp4_310_8 = -(RLcp4_110*ROcp4_56-RLcp4_210*ROcp4_46);
    JTcp4_110_9 = -(RLcp4_210*ROcp4_38-RLcp4_310*ROcp4_28);
    JTcp4_210_9 = RLcp4_110*ROcp4_38-RLcp4_310*ROcp4_18;
    JTcp4_310_9 = -(RLcp4_110*ROcp4_28-RLcp4_210*ROcp4_18);
    OMcp4_110 = OMcp4_19+ROcp4_79*qd[10];
    OMcp4_210 = OMcp4_29+ROcp4_89*qd[10];
    OMcp4_310 = OMcp4_39+ROcp4_99*qd[10];
    ORcp4_110 = OMcp4_29*RLcp4_310-OMcp4_39*RLcp4_210;
    ORcp4_210 = -(OMcp4_19*RLcp4_310-OMcp4_39*RLcp4_110);
    ORcp4_310 = OMcp4_19*RLcp4_210-OMcp4_29*RLcp4_110;
    VIcp4_110 = ORcp4_110+ORcp4_17+ORcp4_18+qd[1]+ROcp4_76*qd[7];
    VIcp4_210 = ORcp4_210+ORcp4_27+ORcp4_28+qd[2]+ROcp4_86*qd[7];
    VIcp4_310 = ORcp4_310+ORcp4_37+ORcp4_38+qd[3]+ROcp4_96*qd[7];
    OPcp4_110 = OPcp4_19+ROcp4_79*qdd[10]+qd[10]*(OMcp4_29*ROcp4_99-OMcp4_39*ROcp4_89);
    OPcp4_210 = OPcp4_29+ROcp4_89*qdd[10]-qd[10]*(OMcp4_19*ROcp4_99-OMcp4_39*ROcp4_79);
    OPcp4_310 = OPcp4_39+ROcp4_99*qdd[10]+qd[10]*(OMcp4_19*ROcp4_89-OMcp4_29*ROcp4_79);
    ACcp4_110 = qdd[1]+OMcp4_26*ORcp4_37+OMcp4_26*ORcp4_38+OMcp4_29*ORcp4_310-OMcp4_36*ORcp4_27-OMcp4_36*ORcp4_28-OMcp4_39
 *ORcp4_210+OPcp4_26*RLcp4_37+OPcp4_26*RLcp4_38+OPcp4_29*RLcp4_310-OPcp4_36*RLcp4_27-OPcp4_36*RLcp4_28-OPcp4_39*RLcp4_210+
 ROcp4_76*qdd[7]+(2.0)*qd[7]*(OMcp4_26*ROcp4_96-OMcp4_36*ROcp4_86);
    ACcp4_210 = qdd[2]-OMcp4_16*ORcp4_37-OMcp4_16*ORcp4_38-OMcp4_19*ORcp4_310+OMcp4_36*ORcp4_17+OMcp4_36*ORcp4_18+OMcp4_39
 *ORcp4_110-OPcp4_16*RLcp4_37-OPcp4_16*RLcp4_38-OPcp4_19*RLcp4_310+OPcp4_36*RLcp4_17+OPcp4_36*RLcp4_18+OPcp4_39*RLcp4_110+
 ROcp4_86*qdd[7]-(2.0)*qd[7]*(OMcp4_16*ROcp4_96-OMcp4_36*ROcp4_76);
    ACcp4_310 = qdd[3]+OMcp4_16*ORcp4_27+OMcp4_16*ORcp4_28+OMcp4_19*ORcp4_210-OMcp4_26*ORcp4_17-OMcp4_26*ORcp4_18-OMcp4_29
 *ORcp4_110+OPcp4_16*RLcp4_27+OPcp4_16*RLcp4_28+OPcp4_19*RLcp4_210-OPcp4_26*RLcp4_17-OPcp4_26*RLcp4_18-OPcp4_29*RLcp4_110+
 ROcp4_96*qdd[7]+(2.0)*qd[7]*(OMcp4_16*ROcp4_86-OMcp4_26*ROcp4_76);

// = = Block_1_0_0_5_1_0 = = 
 
// Symbolic Outputs  

    sens->P[1] = POcp4_110;
    sens->P[2] = POcp4_210;
    sens->P[3] = POcp4_310;
    sens->R[1][1] = ROcp4_110;
    sens->R[1][2] = ROcp4_210;
    sens->R[1][3] = ROcp4_310;
    sens->R[2][1] = ROcp4_410;
    sens->R[2][2] = ROcp4_510;
    sens->R[2][3] = ROcp4_610;
    sens->R[3][1] = ROcp4_79;
    sens->R[3][2] = ROcp4_89;
    sens->R[3][3] = ROcp4_99;
    sens->V[1] = VIcp4_110;
    sens->V[2] = VIcp4_210;
    sens->V[3] = VIcp4_310;
    sens->OM[1] = OMcp4_110;
    sens->OM[2] = OMcp4_210;
    sens->OM[3] = OMcp4_310;
    sens->J[1][1] = (1.0);
    sens->J[1][4] = JTcp4_110_4;
    sens->J[1][5] = JTcp4_110_5;
    sens->J[1][6] = JTcp4_110_6;
    sens->J[1][7] = ROcp4_76;
    sens->J[1][8] = JTcp4_110_8;
    sens->J[1][9] = JTcp4_110_9;
    sens->J[2][2] = (1.0);
    sens->J[2][4] = JTcp4_210_4;
    sens->J[2][5] = JTcp4_210_5;
    sens->J[2][6] = JTcp4_210_6;
    sens->J[2][7] = ROcp4_86;
    sens->J[2][8] = JTcp4_210_8;
    sens->J[2][9] = JTcp4_210_9;
    sens->J[3][3] = (1.0);
    sens->J[3][5] = JTcp4_310_5;
    sens->J[3][6] = JTcp4_310_6;
    sens->J[3][7] = ROcp4_96;
    sens->J[3][8] = JTcp4_310_8;
    sens->J[3][9] = JTcp4_310_9;
    sens->J[4][5] = -S4;
    sens->J[4][6] = ROcp4_15;
    sens->J[4][8] = ROcp4_46;
    sens->J[4][9] = ROcp4_18;
    sens->J[4][10] = ROcp4_79;
    sens->J[5][5] = C4;
    sens->J[5][6] = ROcp4_25;
    sens->J[5][8] = ROcp4_56;
    sens->J[5][9] = ROcp4_28;
    sens->J[5][10] = ROcp4_89;
    sens->J[6][4] = (1.0);
    sens->J[6][6] = -S5;
    sens->J[6][8] = ROcp4_66;
    sens->J[6][9] = ROcp4_38;
    sens->J[6][10] = ROcp4_99;
    sens->A[1] = ACcp4_110;
    sens->A[2] = ACcp4_210;
    sens->A[3] = ACcp4_310;
    sens->OMP[1] = OPcp4_110;
    sens->OMP[2] = OPcp4_210;
    sens->OMP[3] = OPcp4_310;
 
// 
break;
case 6:
 


// = = Block_1_0_0_6_0_1 = = 
 
// Sensor Kinematics 


    ROcp5_15 = C4*C5;
    ROcp5_25 = S4*C5;
    ROcp5_75 = C4*S5;
    ROcp5_85 = S4*S5;
    ROcp5_46 = ROcp5_75*S6-S4*C6;
    ROcp5_56 = ROcp5_85*S6+C4*C6;
    ROcp5_66 = C5*S6;
    ROcp5_76 = ROcp5_75*C6+S4*S6;
    ROcp5_86 = ROcp5_85*C6-C4*S6;
    ROcp5_96 = C5*C6;
    OMcp5_15 = -qd[5]*S4;
    OMcp5_25 = qd[5]*C4;
    OMcp5_16 = OMcp5_15+ROcp5_15*qd[6];
    OMcp5_26 = OMcp5_25+ROcp5_25*qd[6];
    OMcp5_36 = qd[4]-qd[6]*S5;
    OPcp5_16 = ROcp5_15*qdd[6]-qdd[5]*S4-qd[4]*qd[5]*C4-qd[6]*(OMcp5_25*S5+ROcp5_25*qd[4]);
    OPcp5_26 = ROcp5_25*qdd[6]+qdd[5]*C4-qd[4]*qd[5]*S4+qd[6]*(OMcp5_15*S5+ROcp5_15*qd[4]);
    OPcp5_36 = qdd[4]-qdd[6]*S5-qd[5]*qd[6]*C5;

// = = Block_1_0_0_6_0_2 = = 
 
// Sensor Kinematics 


    RLcp5_17 = Dz73*ROcp5_76;
    RLcp5_27 = Dz73*ROcp5_86;
    RLcp5_37 = Dz73*ROcp5_96;
    ORcp5_17 = OMcp5_26*RLcp5_37-OMcp5_36*RLcp5_27;
    ORcp5_27 = -(OMcp5_16*RLcp5_37-OMcp5_36*RLcp5_17);
    ORcp5_37 = OMcp5_16*RLcp5_27-OMcp5_26*RLcp5_17;

// = = Block_1_0_0_6_0_4 = = 
 
// Sensor Kinematics 


    ROcp5_111 = ROcp5_15*C11-ROcp5_76*S11;
    ROcp5_211 = ROcp5_25*C11-ROcp5_86*S11;
    ROcp5_311 = -(ROcp5_96*S11+C11*S5);
    ROcp5_711 = ROcp5_15*S11+ROcp5_76*C11;
    ROcp5_811 = ROcp5_25*S11+ROcp5_86*C11;
    ROcp5_911 = ROcp5_96*C11-S11*S5;
    RLcp5_111 = ROcp5_46*s->dpt[2][7];
    RLcp5_211 = ROcp5_56*s->dpt[2][7];
    RLcp5_311 = ROcp5_66*s->dpt[2][7];
    POcp5_111 = RLcp5_111+RLcp5_17+q[1];
    POcp5_211 = RLcp5_211+RLcp5_27+q[2];
    POcp5_311 = RLcp5_311+RLcp5_37+q[3];
    JTcp5_111_4 = -(RLcp5_211+RLcp5_27);
    JTcp5_211_4 = RLcp5_111+RLcp5_17;
    JTcp5_111_5 = C4*(RLcp5_311+RLcp5_37);
    JTcp5_211_5 = S4*(RLcp5_311+RLcp5_37);
    JTcp5_311_5 = -(C4*(RLcp5_111+RLcp5_17)+S4*(RLcp5_211+RLcp5_27));
    JTcp5_111_6 = ROcp5_25*(RLcp5_311+RLcp5_37)+S5*(RLcp5_211+RLcp5_27);
    JTcp5_211_6 = -(ROcp5_15*(RLcp5_311+RLcp5_37)+S5*(RLcp5_111+RLcp5_17));
    JTcp5_311_6 = ROcp5_15*(RLcp5_211+RLcp5_27)-ROcp5_25*(RLcp5_111+RLcp5_17);
    OMcp5_111 = OMcp5_16+ROcp5_46*qd[11];
    OMcp5_211 = OMcp5_26+ROcp5_56*qd[11];
    OMcp5_311 = OMcp5_36+ROcp5_66*qd[11];
    ORcp5_111 = OMcp5_26*RLcp5_311-OMcp5_36*RLcp5_211;
    ORcp5_211 = -(OMcp5_16*RLcp5_311-OMcp5_36*RLcp5_111);
    ORcp5_311 = OMcp5_16*RLcp5_211-OMcp5_26*RLcp5_111;
    VIcp5_111 = ORcp5_111+ORcp5_17+qd[1]+ROcp5_76*qd[7];
    VIcp5_211 = ORcp5_211+ORcp5_27+qd[2]+ROcp5_86*qd[7];
    VIcp5_311 = ORcp5_311+ORcp5_37+qd[3]+ROcp5_96*qd[7];
    OPcp5_111 = OPcp5_16+ROcp5_46*qdd[11]+qd[11]*(OMcp5_26*ROcp5_66-OMcp5_36*ROcp5_56);
    OPcp5_211 = OPcp5_26+ROcp5_56*qdd[11]-qd[11]*(OMcp5_16*ROcp5_66-OMcp5_36*ROcp5_46);
    OPcp5_311 = OPcp5_36+ROcp5_66*qdd[11]+qd[11]*(OMcp5_16*ROcp5_56-OMcp5_26*ROcp5_46);
    ACcp5_111 = qdd[1]+OMcp5_26*ORcp5_311+OMcp5_26*ORcp5_37-OMcp5_36*ORcp5_211-OMcp5_36*ORcp5_27+OPcp5_26*RLcp5_311+
 OPcp5_26*RLcp5_37-OPcp5_36*RLcp5_211-OPcp5_36*RLcp5_27+ROcp5_76*qdd[7]+(2.0)*qd[7]*(OMcp5_26*ROcp5_96-OMcp5_36*ROcp5_86);
    ACcp5_211 = qdd[2]-OMcp5_16*ORcp5_311-OMcp5_16*ORcp5_37+OMcp5_36*ORcp5_111+OMcp5_36*ORcp5_17-OPcp5_16*RLcp5_311-
 OPcp5_16*RLcp5_37+OPcp5_36*RLcp5_111+OPcp5_36*RLcp5_17+ROcp5_86*qdd[7]-(2.0)*qd[7]*(OMcp5_16*ROcp5_96-OMcp5_36*ROcp5_76);
    ACcp5_311 = qdd[3]+OMcp5_16*ORcp5_211+OMcp5_16*ORcp5_27-OMcp5_26*ORcp5_111-OMcp5_26*ORcp5_17+OPcp5_16*RLcp5_211+
 OPcp5_16*RLcp5_27-OPcp5_26*RLcp5_111-OPcp5_26*RLcp5_17+ROcp5_96*qdd[7]+(2.0)*qd[7]*(OMcp5_16*ROcp5_86-OMcp5_26*ROcp5_76);

// = = Block_1_0_0_6_1_0 = = 
 
// Symbolic Outputs  

    sens->P[1] = POcp5_111;
    sens->P[2] = POcp5_211;
    sens->P[3] = POcp5_311;
    sens->R[1][1] = ROcp5_111;
    sens->R[1][2] = ROcp5_211;
    sens->R[1][3] = ROcp5_311;
    sens->R[2][1] = ROcp5_46;
    sens->R[2][2] = ROcp5_56;
    sens->R[2][3] = ROcp5_66;
    sens->R[3][1] = ROcp5_711;
    sens->R[3][2] = ROcp5_811;
    sens->R[3][3] = ROcp5_911;
    sens->V[1] = VIcp5_111;
    sens->V[2] = VIcp5_211;
    sens->V[3] = VIcp5_311;
    sens->OM[1] = OMcp5_111;
    sens->OM[2] = OMcp5_211;
    sens->OM[3] = OMcp5_311;
    sens->J[1][1] = (1.0);
    sens->J[1][4] = JTcp5_111_4;
    sens->J[1][5] = JTcp5_111_5;
    sens->J[1][6] = JTcp5_111_6;
    sens->J[1][7] = ROcp5_76;
    sens->J[2][2] = (1.0);
    sens->J[2][4] = JTcp5_211_4;
    sens->J[2][5] = JTcp5_211_5;
    sens->J[2][6] = JTcp5_211_6;
    sens->J[2][7] = ROcp5_86;
    sens->J[3][3] = (1.0);
    sens->J[3][5] = JTcp5_311_5;
    sens->J[3][6] = JTcp5_311_6;
    sens->J[3][7] = ROcp5_96;
    sens->J[4][5] = -S4;
    sens->J[4][6] = ROcp5_15;
    sens->J[4][11] = ROcp5_46;
    sens->J[5][5] = C4;
    sens->J[5][6] = ROcp5_25;
    sens->J[5][11] = ROcp5_56;
    sens->J[6][4] = (1.0);
    sens->J[6][6] = -S5;
    sens->J[6][11] = ROcp5_66;
    sens->A[1] = ACcp5_111;
    sens->A[2] = ACcp5_211;
    sens->A[3] = ACcp5_311;
    sens->OMP[1] = OPcp5_111;
    sens->OMP[2] = OPcp5_211;
    sens->OMP[3] = OPcp5_311;
 
// 
break;
case 7:
 


// = = Block_1_0_0_7_0_1 = = 
 
// Sensor Kinematics 


    ROcp6_15 = C4*C5;
    ROcp6_25 = S4*C5;
    ROcp6_75 = C4*S5;
    ROcp6_85 = S4*S5;
    ROcp6_46 = ROcp6_75*S6-S4*C6;
    ROcp6_56 = ROcp6_85*S6+C4*C6;
    ROcp6_66 = C5*S6;
    ROcp6_76 = ROcp6_75*C6+S4*S6;
    ROcp6_86 = ROcp6_85*C6-C4*S6;
    ROcp6_96 = C5*C6;
    OMcp6_15 = -qd[5]*S4;
    OMcp6_25 = qd[5]*C4;
    OMcp6_16 = OMcp6_15+ROcp6_15*qd[6];
    OMcp6_26 = OMcp6_25+ROcp6_25*qd[6];
    OMcp6_36 = qd[4]-qd[6]*S5;
    OPcp6_16 = ROcp6_15*qdd[6]-qdd[5]*S4-qd[4]*qd[5]*C4-qd[6]*(OMcp6_25*S5+ROcp6_25*qd[4]);
    OPcp6_26 = ROcp6_25*qdd[6]+qdd[5]*C4-qd[4]*qd[5]*S4+qd[6]*(OMcp6_15*S5+ROcp6_15*qd[4]);
    OPcp6_36 = qdd[4]-qdd[6]*S5-qd[5]*qd[6]*C5;

// = = Block_1_0_0_7_0_2 = = 
 
// Sensor Kinematics 


    RLcp6_17 = Dz73*ROcp6_76;
    RLcp6_27 = Dz73*ROcp6_86;
    RLcp6_37 = Dz73*ROcp6_96;
    ORcp6_17 = OMcp6_26*RLcp6_37-OMcp6_36*RLcp6_27;
    ORcp6_27 = -(OMcp6_16*RLcp6_37-OMcp6_36*RLcp6_17);
    ORcp6_37 = OMcp6_16*RLcp6_27-OMcp6_26*RLcp6_17;

// = = Block_1_0_0_7_0_4 = = 
 
// Sensor Kinematics 


    ROcp6_111 = ROcp6_15*C11-ROcp6_76*S11;
    ROcp6_211 = ROcp6_25*C11-ROcp6_86*S11;
    ROcp6_311 = -(ROcp6_96*S11+C11*S5);
    ROcp6_711 = ROcp6_15*S11+ROcp6_76*C11;
    ROcp6_811 = ROcp6_25*S11+ROcp6_86*C11;
    ROcp6_911 = ROcp6_96*C11-S11*S5;
    ROcp6_412 = ROcp6_46*C12+ROcp6_711*S12;
    ROcp6_512 = ROcp6_56*C12+ROcp6_811*S12;
    ROcp6_612 = ROcp6_66*C12+ROcp6_911*S12;
    ROcp6_712 = -(ROcp6_46*S12-ROcp6_711*C12);
    ROcp6_812 = -(ROcp6_56*S12-ROcp6_811*C12);
    ROcp6_912 = -(ROcp6_66*S12-ROcp6_911*C12);
    RLcp6_111 = ROcp6_46*s->dpt[2][7];
    RLcp6_211 = ROcp6_56*s->dpt[2][7];
    RLcp6_311 = ROcp6_66*s->dpt[2][7];
    POcp6_111 = RLcp6_111+RLcp6_17+q[1];
    POcp6_211 = RLcp6_211+RLcp6_27+q[2];
    POcp6_311 = RLcp6_311+RLcp6_37+q[3];
    JTcp6_111_4 = -(RLcp6_211+RLcp6_27);
    JTcp6_211_4 = RLcp6_111+RLcp6_17;
    JTcp6_111_5 = C4*(RLcp6_311+RLcp6_37);
    JTcp6_211_5 = S4*(RLcp6_311+RLcp6_37);
    JTcp6_311_5 = -(C4*(RLcp6_111+RLcp6_17)+S4*(RLcp6_211+RLcp6_27));
    JTcp6_111_6 = ROcp6_25*(RLcp6_311+RLcp6_37)+S5*(RLcp6_211+RLcp6_27);
    JTcp6_211_6 = -(ROcp6_15*(RLcp6_311+RLcp6_37)+S5*(RLcp6_111+RLcp6_17));
    JTcp6_311_6 = ROcp6_15*(RLcp6_211+RLcp6_27)-ROcp6_25*(RLcp6_111+RLcp6_17);
    OMcp6_111 = OMcp6_16+ROcp6_46*qd[11];
    OMcp6_211 = OMcp6_26+ROcp6_56*qd[11];
    OMcp6_311 = OMcp6_36+ROcp6_66*qd[11];
    ORcp6_111 = OMcp6_26*RLcp6_311-OMcp6_36*RLcp6_211;
    ORcp6_211 = -(OMcp6_16*RLcp6_311-OMcp6_36*RLcp6_111);
    ORcp6_311 = OMcp6_16*RLcp6_211-OMcp6_26*RLcp6_111;
    VIcp6_111 = ORcp6_111+ORcp6_17+qd[1]+ROcp6_76*qd[7];
    VIcp6_211 = ORcp6_211+ORcp6_27+qd[2]+ROcp6_86*qd[7];
    VIcp6_311 = ORcp6_311+ORcp6_37+qd[3]+ROcp6_96*qd[7];
    ACcp6_111 = qdd[1]+OMcp6_26*ORcp6_311+OMcp6_26*ORcp6_37-OMcp6_36*ORcp6_211-OMcp6_36*ORcp6_27+OPcp6_26*RLcp6_311+
 OPcp6_26*RLcp6_37-OPcp6_36*RLcp6_211-OPcp6_36*RLcp6_27+ROcp6_76*qdd[7]+(2.0)*qd[7]*(OMcp6_26*ROcp6_96-OMcp6_36*ROcp6_86);
    ACcp6_211 = qdd[2]-OMcp6_16*ORcp6_311-OMcp6_16*ORcp6_37+OMcp6_36*ORcp6_111+OMcp6_36*ORcp6_17-OPcp6_16*RLcp6_311-
 OPcp6_16*RLcp6_37+OPcp6_36*RLcp6_111+OPcp6_36*RLcp6_17+ROcp6_86*qdd[7]-(2.0)*qd[7]*(OMcp6_16*ROcp6_96-OMcp6_36*ROcp6_76);
    ACcp6_311 = qdd[3]+OMcp6_16*ORcp6_211+OMcp6_16*ORcp6_27-OMcp6_26*ORcp6_111-OMcp6_26*ORcp6_17+OPcp6_16*RLcp6_211+
 OPcp6_16*RLcp6_27-OPcp6_26*RLcp6_111-OPcp6_26*RLcp6_17+ROcp6_96*qdd[7]+(2.0)*qd[7]*(OMcp6_16*ROcp6_86-OMcp6_26*ROcp6_76);
    OMcp6_112 = OMcp6_111+ROcp6_111*qd[12];
    OMcp6_212 = OMcp6_211+ROcp6_211*qd[12];
    OMcp6_312 = OMcp6_311+ROcp6_311*qd[12];
    OPcp6_112 = OPcp6_16+ROcp6_111*qdd[12]+ROcp6_46*qdd[11]+qd[11]*(OMcp6_26*ROcp6_66-OMcp6_36*ROcp6_56)+qd[12]*(OMcp6_211
 *ROcp6_311-OMcp6_311*ROcp6_211);
    OPcp6_212 = OPcp6_26+ROcp6_211*qdd[12]+ROcp6_56*qdd[11]-qd[11]*(OMcp6_16*ROcp6_66-OMcp6_36*ROcp6_46)-qd[12]*(OMcp6_111
 *ROcp6_311-OMcp6_311*ROcp6_111);
    OPcp6_312 = OPcp6_36+ROcp6_311*qdd[12]+ROcp6_66*qdd[11]+qd[11]*(OMcp6_16*ROcp6_56-OMcp6_26*ROcp6_46)+qd[12]*(OMcp6_111
 *ROcp6_211-OMcp6_211*ROcp6_111);

// = = Block_1_0_0_7_1_0 = = 
 
// Symbolic Outputs  

    sens->P[1] = POcp6_111;
    sens->P[2] = POcp6_211;
    sens->P[3] = POcp6_311;
    sens->R[1][1] = ROcp6_111;
    sens->R[1][2] = ROcp6_211;
    sens->R[1][3] = ROcp6_311;
    sens->R[2][1] = ROcp6_412;
    sens->R[2][2] = ROcp6_512;
    sens->R[2][3] = ROcp6_612;
    sens->R[3][1] = ROcp6_712;
    sens->R[3][2] = ROcp6_812;
    sens->R[3][3] = ROcp6_912;
    sens->V[1] = VIcp6_111;
    sens->V[2] = VIcp6_211;
    sens->V[3] = VIcp6_311;
    sens->OM[1] = OMcp6_112;
    sens->OM[2] = OMcp6_212;
    sens->OM[3] = OMcp6_312;
    sens->J[1][1] = (1.0);
    sens->J[1][4] = JTcp6_111_4;
    sens->J[1][5] = JTcp6_111_5;
    sens->J[1][6] = JTcp6_111_6;
    sens->J[1][7] = ROcp6_76;
    sens->J[2][2] = (1.0);
    sens->J[2][4] = JTcp6_211_4;
    sens->J[2][5] = JTcp6_211_5;
    sens->J[2][6] = JTcp6_211_6;
    sens->J[2][7] = ROcp6_86;
    sens->J[3][3] = (1.0);
    sens->J[3][5] = JTcp6_311_5;
    sens->J[3][6] = JTcp6_311_6;
    sens->J[3][7] = ROcp6_96;
    sens->J[4][5] = -S4;
    sens->J[4][6] = ROcp6_15;
    sens->J[4][11] = ROcp6_46;
    sens->J[4][12] = ROcp6_111;
    sens->J[5][5] = C4;
    sens->J[5][6] = ROcp6_25;
    sens->J[5][11] = ROcp6_56;
    sens->J[5][12] = ROcp6_211;
    sens->J[6][4] = (1.0);
    sens->J[6][6] = -S5;
    sens->J[6][11] = ROcp6_66;
    sens->J[6][12] = ROcp6_311;
    sens->A[1] = ACcp6_111;
    sens->A[2] = ACcp6_211;
    sens->A[3] = ACcp6_311;
    sens->OMP[1] = OPcp6_112;
    sens->OMP[2] = OPcp6_212;
    sens->OMP[3] = OPcp6_312;
 
// 
break;
case 8:
 


// = = Block_1_0_0_8_0_1 = = 
 
// Sensor Kinematics 


    ROcp7_15 = C4*C5;
    ROcp7_25 = S4*C5;
    ROcp7_75 = C4*S5;
    ROcp7_85 = S4*S5;
    ROcp7_46 = ROcp7_75*S6-S4*C6;
    ROcp7_56 = ROcp7_85*S6+C4*C6;
    ROcp7_66 = C5*S6;
    ROcp7_76 = ROcp7_75*C6+S4*S6;
    ROcp7_86 = ROcp7_85*C6-C4*S6;
    ROcp7_96 = C5*C6;
    OMcp7_15 = -qd[5]*S4;
    OMcp7_25 = qd[5]*C4;
    OMcp7_16 = OMcp7_15+ROcp7_15*qd[6];
    OMcp7_26 = OMcp7_25+ROcp7_25*qd[6];
    OMcp7_36 = qd[4]-qd[6]*S5;
    OPcp7_16 = ROcp7_15*qdd[6]-qdd[5]*S4-qd[4]*qd[5]*C4-qd[6]*(OMcp7_25*S5+ROcp7_25*qd[4]);
    OPcp7_26 = ROcp7_25*qdd[6]+qdd[5]*C4-qd[4]*qd[5]*S4+qd[6]*(OMcp7_15*S5+ROcp7_15*qd[4]);
    OPcp7_36 = qdd[4]-qdd[6]*S5-qd[5]*qd[6]*C5;

// = = Block_1_0_0_8_0_2 = = 
 
// Sensor Kinematics 


    RLcp7_17 = Dz73*ROcp7_76;
    RLcp7_27 = Dz73*ROcp7_86;
    RLcp7_37 = Dz73*ROcp7_96;
    ORcp7_17 = OMcp7_26*RLcp7_37-OMcp7_36*RLcp7_27;
    ORcp7_27 = -(OMcp7_16*RLcp7_37-OMcp7_36*RLcp7_17);
    ORcp7_37 = OMcp7_16*RLcp7_27-OMcp7_26*RLcp7_17;

// = = Block_1_0_0_8_0_4 = = 
 
// Sensor Kinematics 


    ROcp7_111 = ROcp7_15*C11-ROcp7_76*S11;
    ROcp7_211 = ROcp7_25*C11-ROcp7_86*S11;
    ROcp7_311 = -(ROcp7_96*S11+C11*S5);
    ROcp7_711 = ROcp7_15*S11+ROcp7_76*C11;
    ROcp7_811 = ROcp7_25*S11+ROcp7_86*C11;
    ROcp7_911 = ROcp7_96*C11-S11*S5;
    ROcp7_412 = ROcp7_46*C12+ROcp7_711*S12;
    ROcp7_512 = ROcp7_56*C12+ROcp7_811*S12;
    ROcp7_612 = ROcp7_66*C12+ROcp7_911*S12;
    ROcp7_712 = -(ROcp7_46*S12-ROcp7_711*C12);
    ROcp7_812 = -(ROcp7_56*S12-ROcp7_811*C12);
    ROcp7_912 = -(ROcp7_66*S12-ROcp7_911*C12);
    ROcp7_113 = ROcp7_111*C13+ROcp7_412*S13;
    ROcp7_213 = ROcp7_211*C13+ROcp7_512*S13;
    ROcp7_313 = ROcp7_311*C13+ROcp7_612*S13;
    ROcp7_413 = -(ROcp7_111*S13-ROcp7_412*C13);
    ROcp7_513 = -(ROcp7_211*S13-ROcp7_512*C13);
    ROcp7_613 = -(ROcp7_311*S13-ROcp7_612*C13);
    RLcp7_111 = ROcp7_46*s->dpt[2][7];
    RLcp7_211 = ROcp7_56*s->dpt[2][7];
    RLcp7_311 = ROcp7_66*s->dpt[2][7];
    OMcp7_111 = OMcp7_16+ROcp7_46*qd[11];
    OMcp7_211 = OMcp7_26+ROcp7_56*qd[11];
    OMcp7_311 = OMcp7_36+ROcp7_66*qd[11];
    ORcp7_111 = OMcp7_26*RLcp7_311-OMcp7_36*RLcp7_211;
    ORcp7_211 = -(OMcp7_16*RLcp7_311-OMcp7_36*RLcp7_111);
    ORcp7_311 = OMcp7_16*RLcp7_211-OMcp7_26*RLcp7_111;
    OMcp7_112 = OMcp7_111+ROcp7_111*qd[12];
    OMcp7_212 = OMcp7_211+ROcp7_211*qd[12];
    OMcp7_312 = OMcp7_311+ROcp7_311*qd[12];
    OPcp7_112 = OPcp7_16+ROcp7_111*qdd[12]+ROcp7_46*qdd[11]+qd[11]*(OMcp7_26*ROcp7_66-OMcp7_36*ROcp7_56)+qd[12]*(OMcp7_211
 *ROcp7_311-OMcp7_311*ROcp7_211);
    OPcp7_212 = OPcp7_26+ROcp7_211*qdd[12]+ROcp7_56*qdd[11]-qd[11]*(OMcp7_16*ROcp7_66-OMcp7_36*ROcp7_46)-qd[12]*(OMcp7_111
 *ROcp7_311-OMcp7_311*ROcp7_111);
    OPcp7_312 = OPcp7_36+ROcp7_311*qdd[12]+ROcp7_66*qdd[11]+qd[11]*(OMcp7_16*ROcp7_56-OMcp7_26*ROcp7_46)+qd[12]*(OMcp7_111
 *ROcp7_211-OMcp7_211*ROcp7_111);
    RLcp7_113 = ROcp7_712*s->dpt[3][15];
    RLcp7_213 = ROcp7_812*s->dpt[3][15];
    RLcp7_313 = ROcp7_912*s->dpt[3][15];
    POcp7_113 = RLcp7_111+RLcp7_113+RLcp7_17+q[1];
    POcp7_213 = RLcp7_211+RLcp7_213+RLcp7_27+q[2];
    POcp7_313 = RLcp7_311+RLcp7_313+RLcp7_37+q[3];
    JTcp7_113_4 = -(RLcp7_211+RLcp7_213+RLcp7_27);
    JTcp7_213_4 = RLcp7_111+RLcp7_113+RLcp7_17;
    JTcp7_113_5 = C4*(RLcp7_311+RLcp7_313+RLcp7_37);
    JTcp7_213_5 = S4*(RLcp7_311+RLcp7_313+RLcp7_37);
    JTcp7_313_5 = -(RLcp7_213*S4+C4*(RLcp7_111+RLcp7_113+RLcp7_17)+S4*(RLcp7_211+RLcp7_27));
    JTcp7_113_6 = RLcp7_213*S5+RLcp7_313*ROcp7_25+ROcp7_25*(RLcp7_311+RLcp7_37)+S5*(RLcp7_211+RLcp7_27);
    JTcp7_213_6 = -(RLcp7_113*S5+RLcp7_313*ROcp7_15+ROcp7_15*(RLcp7_311+RLcp7_37)+S5*(RLcp7_111+RLcp7_17));
    JTcp7_313_6 = ROcp7_15*(RLcp7_211+RLcp7_27)-ROcp7_25*(RLcp7_111+RLcp7_17)-RLcp7_113*ROcp7_25+RLcp7_213*ROcp7_15;
    JTcp7_113_8 = -(RLcp7_213*ROcp7_66-RLcp7_313*ROcp7_56);
    JTcp7_213_8 = RLcp7_113*ROcp7_66-RLcp7_313*ROcp7_46;
    JTcp7_313_8 = -(RLcp7_113*ROcp7_56-RLcp7_213*ROcp7_46);
    JTcp7_113_9 = -(RLcp7_213*ROcp7_311-RLcp7_313*ROcp7_211);
    JTcp7_213_9 = RLcp7_113*ROcp7_311-RLcp7_313*ROcp7_111;
    JTcp7_313_9 = -(RLcp7_113*ROcp7_211-RLcp7_213*ROcp7_111);
    OMcp7_113 = OMcp7_112+ROcp7_712*qd[13];
    OMcp7_213 = OMcp7_212+ROcp7_812*qd[13];
    OMcp7_313 = OMcp7_312+ROcp7_912*qd[13];
    ORcp7_113 = OMcp7_212*RLcp7_313-OMcp7_312*RLcp7_213;
    ORcp7_213 = -(OMcp7_112*RLcp7_313-OMcp7_312*RLcp7_113);
    ORcp7_313 = OMcp7_112*RLcp7_213-OMcp7_212*RLcp7_113;
    VIcp7_113 = ORcp7_111+ORcp7_113+ORcp7_17+qd[1]+ROcp7_76*qd[7];
    VIcp7_213 = ORcp7_211+ORcp7_213+ORcp7_27+qd[2]+ROcp7_86*qd[7];
    VIcp7_313 = ORcp7_311+ORcp7_313+ORcp7_37+qd[3]+ROcp7_96*qd[7];
    OPcp7_113 = OPcp7_112+ROcp7_712*qdd[13]+qd[13]*(OMcp7_212*ROcp7_912-OMcp7_312*ROcp7_812);
    OPcp7_213 = OPcp7_212+ROcp7_812*qdd[13]-qd[13]*(OMcp7_112*ROcp7_912-OMcp7_312*ROcp7_712);
    OPcp7_313 = OPcp7_312+ROcp7_912*qdd[13]+qd[13]*(OMcp7_112*ROcp7_812-OMcp7_212*ROcp7_712);
    ACcp7_113 = qdd[1]+OMcp7_212*ORcp7_313+OMcp7_26*ORcp7_311+OMcp7_26*ORcp7_37-OMcp7_312*ORcp7_213-OMcp7_36*ORcp7_211-
 OMcp7_36*ORcp7_27+OPcp7_212*RLcp7_313+OPcp7_26*RLcp7_311+OPcp7_26*RLcp7_37-OPcp7_312*RLcp7_213-OPcp7_36*RLcp7_211-OPcp7_36*
 RLcp7_27+ROcp7_76*qdd[7]+(2.0)*qd[7]*(OMcp7_26*ROcp7_96-OMcp7_36*ROcp7_86);
    ACcp7_213 = qdd[2]-OMcp7_112*ORcp7_313-OMcp7_16*ORcp7_311-OMcp7_16*ORcp7_37+OMcp7_312*ORcp7_113+OMcp7_36*ORcp7_111+
 OMcp7_36*ORcp7_17-OPcp7_112*RLcp7_313-OPcp7_16*RLcp7_311-OPcp7_16*RLcp7_37+OPcp7_312*RLcp7_113+OPcp7_36*RLcp7_111+OPcp7_36*
 RLcp7_17+ROcp7_86*qdd[7]-(2.0)*qd[7]*(OMcp7_16*ROcp7_96-OMcp7_36*ROcp7_76);
    ACcp7_313 = qdd[3]+OMcp7_112*ORcp7_213+OMcp7_16*ORcp7_211+OMcp7_16*ORcp7_27-OMcp7_212*ORcp7_113-OMcp7_26*ORcp7_111-
 OMcp7_26*ORcp7_17+OPcp7_112*RLcp7_213+OPcp7_16*RLcp7_211+OPcp7_16*RLcp7_27-OPcp7_212*RLcp7_113-OPcp7_26*RLcp7_111-OPcp7_26*
 RLcp7_17+ROcp7_96*qdd[7]+(2.0)*qd[7]*(OMcp7_16*ROcp7_86-OMcp7_26*ROcp7_76);

// = = Block_1_0_0_8_1_0 = = 
 
// Symbolic Outputs  

    sens->P[1] = POcp7_113;
    sens->P[2] = POcp7_213;
    sens->P[3] = POcp7_313;
    sens->R[1][1] = ROcp7_113;
    sens->R[1][2] = ROcp7_213;
    sens->R[1][3] = ROcp7_313;
    sens->R[2][1] = ROcp7_413;
    sens->R[2][2] = ROcp7_513;
    sens->R[2][3] = ROcp7_613;
    sens->R[3][1] = ROcp7_712;
    sens->R[3][2] = ROcp7_812;
    sens->R[3][3] = ROcp7_912;
    sens->V[1] = VIcp7_113;
    sens->V[2] = VIcp7_213;
    sens->V[3] = VIcp7_313;
    sens->OM[1] = OMcp7_113;
    sens->OM[2] = OMcp7_213;
    sens->OM[3] = OMcp7_313;
    sens->J[1][1] = (1.0);
    sens->J[1][4] = JTcp7_113_4;
    sens->J[1][5] = JTcp7_113_5;
    sens->J[1][6] = JTcp7_113_6;
    sens->J[1][7] = ROcp7_76;
    sens->J[1][11] = JTcp7_113_8;
    sens->J[1][12] = JTcp7_113_9;
    sens->J[2][2] = (1.0);
    sens->J[2][4] = JTcp7_213_4;
    sens->J[2][5] = JTcp7_213_5;
    sens->J[2][6] = JTcp7_213_6;
    sens->J[2][7] = ROcp7_86;
    sens->J[2][11] = JTcp7_213_8;
    sens->J[2][12] = JTcp7_213_9;
    sens->J[3][3] = (1.0);
    sens->J[3][5] = JTcp7_313_5;
    sens->J[3][6] = JTcp7_313_6;
    sens->J[3][7] = ROcp7_96;
    sens->J[3][11] = JTcp7_313_8;
    sens->J[3][12] = JTcp7_313_9;
    sens->J[4][5] = -S4;
    sens->J[4][6] = ROcp7_15;
    sens->J[4][11] = ROcp7_46;
    sens->J[4][12] = ROcp7_111;
    sens->J[4][13] = ROcp7_712;
    sens->J[5][5] = C4;
    sens->J[5][6] = ROcp7_25;
    sens->J[5][11] = ROcp7_56;
    sens->J[5][12] = ROcp7_211;
    sens->J[5][13] = ROcp7_812;
    sens->J[6][4] = (1.0);
    sens->J[6][6] = -S5;
    sens->J[6][11] = ROcp7_66;
    sens->J[6][12] = ROcp7_311;
    sens->J[6][13] = ROcp7_912;
    sens->A[1] = ACcp7_113;
    sens->A[2] = ACcp7_213;
    sens->A[3] = ACcp7_313;
    sens->OMP[1] = OPcp7_113;
    sens->OMP[2] = OPcp7_213;
    sens->OMP[3] = OPcp7_313;
 
// 
break;
case 9:
 


// = = Block_1_0_0_9_0_1 = = 
 
// Sensor Kinematics 


    ROcp8_15 = C4*C5;
    ROcp8_25 = S4*C5;
    ROcp8_75 = C4*S5;
    ROcp8_85 = S4*S5;
    ROcp8_46 = ROcp8_75*S6-S4*C6;
    ROcp8_56 = ROcp8_85*S6+C4*C6;
    ROcp8_66 = C5*S6;
    ROcp8_76 = ROcp8_75*C6+S4*S6;
    ROcp8_86 = ROcp8_85*C6-C4*S6;
    ROcp8_96 = C5*C6;
    OMcp8_15 = -qd[5]*S4;
    OMcp8_25 = qd[5]*C4;
    OMcp8_16 = OMcp8_15+ROcp8_15*qd[6];
    OMcp8_26 = OMcp8_25+ROcp8_25*qd[6];
    OMcp8_36 = qd[4]-qd[6]*S5;
    OPcp8_16 = ROcp8_15*qdd[6]-qdd[5]*S4-qd[4]*qd[5]*C4-qd[6]*(OMcp8_25*S5+ROcp8_25*qd[4]);
    OPcp8_26 = ROcp8_25*qdd[6]+qdd[5]*C4-qd[4]*qd[5]*S4+qd[6]*(OMcp8_15*S5+ROcp8_15*qd[4]);
    OPcp8_36 = qdd[4]-qdd[6]*S5-qd[5]*qd[6]*C5;

// = = Block_1_0_0_9_0_2 = = 
 
// Sensor Kinematics 


    RLcp8_17 = Dz73*ROcp8_76;
    RLcp8_27 = Dz73*ROcp8_86;
    RLcp8_37 = Dz73*ROcp8_96;
    ORcp8_17 = OMcp8_26*RLcp8_37-OMcp8_36*RLcp8_27;
    ORcp8_27 = -(OMcp8_16*RLcp8_37-OMcp8_36*RLcp8_17);
    ORcp8_37 = OMcp8_16*RLcp8_27-OMcp8_26*RLcp8_17;

// = = Block_1_0_0_9_0_5 = = 
 
// Sensor Kinematics 


    ROcp8_114 = ROcp8_15*C14+ROcp8_46*S14;
    ROcp8_214 = ROcp8_25*C14+ROcp8_56*S14;
    ROcp8_314 = ROcp8_66*S14-C14*S5;
    ROcp8_414 = -(ROcp8_15*S14-ROcp8_46*C14);
    ROcp8_514 = -(ROcp8_25*S14-ROcp8_56*C14);
    ROcp8_614 = ROcp8_66*C14+S14*S5;
    RLcp8_114 = ROcp8_76*s->dpt[3][8];
    RLcp8_214 = ROcp8_86*s->dpt[3][8];
    RLcp8_314 = ROcp8_96*s->dpt[3][8];
    POcp8_114 = RLcp8_114+RLcp8_17+q[1];
    POcp8_214 = RLcp8_214+RLcp8_27+q[2];
    POcp8_314 = RLcp8_314+RLcp8_37+q[3];
    JTcp8_114_4 = -(RLcp8_214+RLcp8_27);
    JTcp8_214_4 = RLcp8_114+RLcp8_17;
    JTcp8_114_5 = C4*(RLcp8_314+RLcp8_37);
    JTcp8_214_5 = S4*(RLcp8_314+RLcp8_37);
    JTcp8_314_5 = -(C4*(RLcp8_114+RLcp8_17)+S4*(RLcp8_214+RLcp8_27));
    JTcp8_114_6 = ROcp8_25*(RLcp8_314+RLcp8_37)+S5*(RLcp8_214+RLcp8_27);
    JTcp8_214_6 = -(ROcp8_15*(RLcp8_314+RLcp8_37)+S5*(RLcp8_114+RLcp8_17));
    JTcp8_314_6 = ROcp8_15*(RLcp8_214+RLcp8_27)-ROcp8_25*(RLcp8_114+RLcp8_17);
    OMcp8_114 = OMcp8_16+ROcp8_76*qd[14];
    OMcp8_214 = OMcp8_26+ROcp8_86*qd[14];
    OMcp8_314 = OMcp8_36+ROcp8_96*qd[14];
    ORcp8_114 = OMcp8_26*RLcp8_314-OMcp8_36*RLcp8_214;
    ORcp8_214 = -(OMcp8_16*RLcp8_314-OMcp8_36*RLcp8_114);
    ORcp8_314 = OMcp8_16*RLcp8_214-OMcp8_26*RLcp8_114;
    VIcp8_114 = ORcp8_114+ORcp8_17+qd[1]+ROcp8_76*qd[7];
    VIcp8_214 = ORcp8_214+ORcp8_27+qd[2]+ROcp8_86*qd[7];
    VIcp8_314 = ORcp8_314+ORcp8_37+qd[3]+ROcp8_96*qd[7];
    OPcp8_114 = OPcp8_16+ROcp8_76*qdd[14]+qd[14]*(OMcp8_26*ROcp8_96-OMcp8_36*ROcp8_86);
    OPcp8_214 = OPcp8_26+ROcp8_86*qdd[14]-qd[14]*(OMcp8_16*ROcp8_96-OMcp8_36*ROcp8_76);
    OPcp8_314 = OPcp8_36+ROcp8_96*qdd[14]+qd[14]*(OMcp8_16*ROcp8_86-OMcp8_26*ROcp8_76);
    ACcp8_114 = qdd[1]+OMcp8_26*ORcp8_314+OMcp8_26*ORcp8_37-OMcp8_36*ORcp8_214-OMcp8_36*ORcp8_27+OPcp8_26*RLcp8_314+
 OPcp8_26*RLcp8_37-OPcp8_36*RLcp8_214-OPcp8_36*RLcp8_27+ROcp8_76*qdd[7]+(2.0)*qd[7]*(OMcp8_26*ROcp8_96-OMcp8_36*ROcp8_86);
    ACcp8_214 = qdd[2]-OMcp8_16*ORcp8_314-OMcp8_16*ORcp8_37+OMcp8_36*ORcp8_114+OMcp8_36*ORcp8_17-OPcp8_16*RLcp8_314-
 OPcp8_16*RLcp8_37+OPcp8_36*RLcp8_114+OPcp8_36*RLcp8_17+ROcp8_86*qdd[7]-(2.0)*qd[7]*(OMcp8_16*ROcp8_96-OMcp8_36*ROcp8_76);
    ACcp8_314 = qdd[3]+OMcp8_16*ORcp8_214+OMcp8_16*ORcp8_27-OMcp8_26*ORcp8_114-OMcp8_26*ORcp8_17+OPcp8_16*RLcp8_214+
 OPcp8_16*RLcp8_27-OPcp8_26*RLcp8_114-OPcp8_26*RLcp8_17+ROcp8_96*qdd[7]+(2.0)*qd[7]*(OMcp8_16*ROcp8_86-OMcp8_26*ROcp8_76);

// = = Block_1_0_0_9_1_0 = = 
 
// Symbolic Outputs  

    sens->P[1] = POcp8_114;
    sens->P[2] = POcp8_214;
    sens->P[3] = POcp8_314;
    sens->R[1][1] = ROcp8_114;
    sens->R[1][2] = ROcp8_214;
    sens->R[1][3] = ROcp8_314;
    sens->R[2][1] = ROcp8_414;
    sens->R[2][2] = ROcp8_514;
    sens->R[2][3] = ROcp8_614;
    sens->R[3][1] = ROcp8_76;
    sens->R[3][2] = ROcp8_86;
    sens->R[3][3] = ROcp8_96;
    sens->V[1] = VIcp8_114;
    sens->V[2] = VIcp8_214;
    sens->V[3] = VIcp8_314;
    sens->OM[1] = OMcp8_114;
    sens->OM[2] = OMcp8_214;
    sens->OM[3] = OMcp8_314;
    sens->J[1][1] = (1.0);
    sens->J[1][4] = JTcp8_114_4;
    sens->J[1][5] = JTcp8_114_5;
    sens->J[1][6] = JTcp8_114_6;
    sens->J[1][7] = ROcp8_76;
    sens->J[2][2] = (1.0);
    sens->J[2][4] = JTcp8_214_4;
    sens->J[2][5] = JTcp8_214_5;
    sens->J[2][6] = JTcp8_214_6;
    sens->J[2][7] = ROcp8_86;
    sens->J[3][3] = (1.0);
    sens->J[3][5] = JTcp8_314_5;
    sens->J[3][6] = JTcp8_314_6;
    sens->J[3][7] = ROcp8_96;
    sens->J[4][5] = -S4;
    sens->J[4][6] = ROcp8_15;
    sens->J[4][14] = ROcp8_76;
    sens->J[5][5] = C4;
    sens->J[5][6] = ROcp8_25;
    sens->J[5][14] = ROcp8_86;
    sens->J[6][4] = (1.0);
    sens->J[6][6] = -S5;
    sens->J[6][14] = ROcp8_96;
    sens->A[1] = ACcp8_114;
    sens->A[2] = ACcp8_214;
    sens->A[3] = ACcp8_314;
    sens->OMP[1] = OPcp8_114;
    sens->OMP[2] = OPcp8_214;
    sens->OMP[3] = OPcp8_314;
 
// 
break;
case 10:
 


// = = Block_1_0_0_10_0_1 = = 
 
// Sensor Kinematics 


    ROcp9_15 = C4*C5;
    ROcp9_25 = S4*C5;
    ROcp9_75 = C4*S5;
    ROcp9_85 = S4*S5;
    ROcp9_46 = ROcp9_75*S6-S4*C6;
    ROcp9_56 = ROcp9_85*S6+C4*C6;
    ROcp9_66 = C5*S6;
    ROcp9_76 = ROcp9_75*C6+S4*S6;
    ROcp9_86 = ROcp9_85*C6-C4*S6;
    ROcp9_96 = C5*C6;
    OMcp9_15 = -qd[5]*S4;
    OMcp9_25 = qd[5]*C4;
    OMcp9_16 = OMcp9_15+ROcp9_15*qd[6];
    OMcp9_26 = OMcp9_25+ROcp9_25*qd[6];
    OMcp9_36 = qd[4]-qd[6]*S5;
    OPcp9_16 = ROcp9_15*qdd[6]-qdd[5]*S4-qd[4]*qd[5]*C4-qd[6]*(OMcp9_25*S5+ROcp9_25*qd[4]);
    OPcp9_26 = ROcp9_25*qdd[6]+qdd[5]*C4-qd[4]*qd[5]*S4+qd[6]*(OMcp9_15*S5+ROcp9_15*qd[4]);
    OPcp9_36 = qdd[4]-qdd[6]*S5-qd[5]*qd[6]*C5;

// = = Block_1_0_0_10_0_6 = = 
 
// Sensor Kinematics 


    ROcp9_116 = ROcp9_15*C16-ROcp9_76*S16;
    ROcp9_216 = ROcp9_25*C16-ROcp9_86*S16;
    ROcp9_316 = -(ROcp9_96*S16+C16*S5);
    ROcp9_716 = ROcp9_15*S16+ROcp9_76*C16;
    ROcp9_816 = ROcp9_25*S16+ROcp9_86*C16;
    ROcp9_916 = ROcp9_96*C16-S16*S5;
    RLcp9_116 = ROcp9_46*s->dpt[2][3];
    RLcp9_216 = ROcp9_56*s->dpt[2][3];
    RLcp9_316 = ROcp9_66*s->dpt[2][3];
    POcp9_116 = RLcp9_116+q[1];
    POcp9_216 = RLcp9_216+q[2];
    POcp9_316 = RLcp9_316+q[3];
    JTcp9_116_5 = RLcp9_316*C4;
    JTcp9_216_5 = RLcp9_316*S4;
    JTcp9_316_5 = -(RLcp9_116*C4+RLcp9_216*S4);
    JTcp9_116_6 = RLcp9_216*S5+RLcp9_316*ROcp9_25;
    JTcp9_216_6 = -(RLcp9_116*S5+RLcp9_316*ROcp9_15);
    JTcp9_316_6 = -(RLcp9_116*ROcp9_25-RLcp9_216*ROcp9_15);
    OMcp9_116 = OMcp9_16+ROcp9_46*qd[16];
    OMcp9_216 = OMcp9_26+ROcp9_56*qd[16];
    OMcp9_316 = OMcp9_36+ROcp9_66*qd[16];
    ORcp9_116 = OMcp9_26*RLcp9_316-OMcp9_36*RLcp9_216;
    ORcp9_216 = -(OMcp9_16*RLcp9_316-OMcp9_36*RLcp9_116);
    ORcp9_316 = OMcp9_16*RLcp9_216-OMcp9_26*RLcp9_116;
    VIcp9_116 = ORcp9_116+qd[1];
    VIcp9_216 = ORcp9_216+qd[2];
    VIcp9_316 = ORcp9_316+qd[3];
    OPcp9_116 = OPcp9_16+ROcp9_46*qdd[16]+qd[16]*(OMcp9_26*ROcp9_66-OMcp9_36*ROcp9_56);
    OPcp9_216 = OPcp9_26+ROcp9_56*qdd[16]-qd[16]*(OMcp9_16*ROcp9_66-OMcp9_36*ROcp9_46);
    OPcp9_316 = OPcp9_36+ROcp9_66*qdd[16]+qd[16]*(OMcp9_16*ROcp9_56-OMcp9_26*ROcp9_46);
    ACcp9_116 = qdd[1]+OMcp9_26*ORcp9_316-OMcp9_36*ORcp9_216+OPcp9_26*RLcp9_316-OPcp9_36*RLcp9_216;
    ACcp9_216 = qdd[2]-OMcp9_16*ORcp9_316+OMcp9_36*ORcp9_116-OPcp9_16*RLcp9_316+OPcp9_36*RLcp9_116;
    ACcp9_316 = qdd[3]+OMcp9_16*ORcp9_216-OMcp9_26*ORcp9_116+OPcp9_16*RLcp9_216-OPcp9_26*RLcp9_116;

// = = Block_1_0_0_10_1_0 = = 
 
// Symbolic Outputs  

    sens->P[1] = POcp9_116;
    sens->P[2] = POcp9_216;
    sens->P[3] = POcp9_316;
    sens->R[1][1] = ROcp9_116;
    sens->R[1][2] = ROcp9_216;
    sens->R[1][3] = ROcp9_316;
    sens->R[2][1] = ROcp9_46;
    sens->R[2][2] = ROcp9_56;
    sens->R[2][3] = ROcp9_66;
    sens->R[3][1] = ROcp9_716;
    sens->R[3][2] = ROcp9_816;
    sens->R[3][3] = ROcp9_916;
    sens->V[1] = VIcp9_116;
    sens->V[2] = VIcp9_216;
    sens->V[3] = VIcp9_316;
    sens->OM[1] = OMcp9_116;
    sens->OM[2] = OMcp9_216;
    sens->OM[3] = OMcp9_316;
    sens->J[1][1] = (1.0);
    sens->J[1][4] = -RLcp9_216;
    sens->J[1][5] = JTcp9_116_5;
    sens->J[1][6] = JTcp9_116_6;
    sens->J[2][2] = (1.0);
    sens->J[2][4] = RLcp9_116;
    sens->J[2][5] = JTcp9_216_5;
    sens->J[2][6] = JTcp9_216_6;
    sens->J[3][3] = (1.0);
    sens->J[3][5] = JTcp9_316_5;
    sens->J[3][6] = JTcp9_316_6;
    sens->J[4][5] = -S4;
    sens->J[4][6] = ROcp9_15;
    sens->J[4][16] = ROcp9_46;
    sens->J[5][5] = C4;
    sens->J[5][6] = ROcp9_25;
    sens->J[5][16] = ROcp9_56;
    sens->J[6][4] = (1.0);
    sens->J[6][6] = -S5;
    sens->J[6][16] = ROcp9_66;
    sens->A[1] = ACcp9_116;
    sens->A[2] = ACcp9_216;
    sens->A[3] = ACcp9_316;
    sens->OMP[1] = OPcp9_116;
    sens->OMP[2] = OPcp9_216;
    sens->OMP[3] = OPcp9_316;
 
// 
break;
case 11:
 


// = = Block_1_0_0_11_0_1 = = 
 
// Sensor Kinematics 


    ROcp10_15 = C4*C5;
    ROcp10_25 = S4*C5;
    ROcp10_75 = C4*S5;
    ROcp10_85 = S4*S5;
    ROcp10_46 = ROcp10_75*S6-S4*C6;
    ROcp10_56 = ROcp10_85*S6+C4*C6;
    ROcp10_66 = C5*S6;
    ROcp10_76 = ROcp10_75*C6+S4*S6;
    ROcp10_86 = ROcp10_85*C6-C4*S6;
    ROcp10_96 = C5*C6;
    OMcp10_15 = -qd[5]*S4;
    OMcp10_25 = qd[5]*C4;
    OMcp10_16 = OMcp10_15+ROcp10_15*qd[6];
    OMcp10_26 = OMcp10_25+ROcp10_25*qd[6];
    OMcp10_36 = qd[4]-qd[6]*S5;
    OPcp10_16 = ROcp10_15*qdd[6]-qdd[5]*S4-qd[4]*qd[5]*C4-qd[6]*(OMcp10_25*S5+ROcp10_25*qd[4]);
    OPcp10_26 = ROcp10_25*qdd[6]+qdd[5]*C4-qd[4]*qd[5]*S4+qd[6]*(OMcp10_15*S5+ROcp10_15*qd[4]);
    OPcp10_36 = qdd[4]-qdd[6]*S5-qd[5]*qd[6]*C5;

// = = Block_1_0_0_11_0_7 = = 
 
// Sensor Kinematics 


    ROcp10_117 = ROcp10_15*C17-ROcp10_76*S17;
    ROcp10_217 = ROcp10_25*C17-ROcp10_86*S17;
    ROcp10_317 = -(ROcp10_96*S17+C17*S5);
    ROcp10_717 = ROcp10_15*S17+ROcp10_76*C17;
    ROcp10_817 = ROcp10_25*S17+ROcp10_86*C17;
    ROcp10_917 = ROcp10_96*C17-S17*S5;
    RLcp10_117 = ROcp10_46*s->dpt[2][4];
    RLcp10_217 = ROcp10_56*s->dpt[2][4];
    RLcp10_317 = ROcp10_66*s->dpt[2][4];
    POcp10_117 = RLcp10_117+q[1];
    POcp10_217 = RLcp10_217+q[2];
    POcp10_317 = RLcp10_317+q[3];
    JTcp10_117_5 = RLcp10_317*C4;
    JTcp10_217_5 = RLcp10_317*S4;
    JTcp10_317_5 = -(RLcp10_117*C4+RLcp10_217*S4);
    JTcp10_117_6 = RLcp10_217*S5+RLcp10_317*ROcp10_25;
    JTcp10_217_6 = -(RLcp10_117*S5+RLcp10_317*ROcp10_15);
    JTcp10_317_6 = -(RLcp10_117*ROcp10_25-RLcp10_217*ROcp10_15);
    OMcp10_117 = OMcp10_16+ROcp10_46*qd[17];
    OMcp10_217 = OMcp10_26+ROcp10_56*qd[17];
    OMcp10_317 = OMcp10_36+ROcp10_66*qd[17];
    ORcp10_117 = OMcp10_26*RLcp10_317-OMcp10_36*RLcp10_217;
    ORcp10_217 = -(OMcp10_16*RLcp10_317-OMcp10_36*RLcp10_117);
    ORcp10_317 = OMcp10_16*RLcp10_217-OMcp10_26*RLcp10_117;
    VIcp10_117 = ORcp10_117+qd[1];
    VIcp10_217 = ORcp10_217+qd[2];
    VIcp10_317 = ORcp10_317+qd[3];
    OPcp10_117 = OPcp10_16+ROcp10_46*qdd[17]+qd[17]*(OMcp10_26*ROcp10_66-OMcp10_36*ROcp10_56);
    OPcp10_217 = OPcp10_26+ROcp10_56*qdd[17]-qd[17]*(OMcp10_16*ROcp10_66-OMcp10_36*ROcp10_46);
    OPcp10_317 = OPcp10_36+ROcp10_66*qdd[17]+qd[17]*(OMcp10_16*ROcp10_56-OMcp10_26*ROcp10_46);
    ACcp10_117 = qdd[1]+OMcp10_26*ORcp10_317-OMcp10_36*ORcp10_217+OPcp10_26*RLcp10_317-OPcp10_36*RLcp10_217;
    ACcp10_217 = qdd[2]-OMcp10_16*ORcp10_317+OMcp10_36*ORcp10_117-OPcp10_16*RLcp10_317+OPcp10_36*RLcp10_117;
    ACcp10_317 = qdd[3]+OMcp10_16*ORcp10_217-OMcp10_26*ORcp10_117+OPcp10_16*RLcp10_217-OPcp10_26*RLcp10_117;

// = = Block_1_0_0_11_1_0 = = 
 
// Symbolic Outputs  

    sens->P[1] = POcp10_117;
    sens->P[2] = POcp10_217;
    sens->P[3] = POcp10_317;
    sens->R[1][1] = ROcp10_117;
    sens->R[1][2] = ROcp10_217;
    sens->R[1][3] = ROcp10_317;
    sens->R[2][1] = ROcp10_46;
    sens->R[2][2] = ROcp10_56;
    sens->R[2][3] = ROcp10_66;
    sens->R[3][1] = ROcp10_717;
    sens->R[3][2] = ROcp10_817;
    sens->R[3][3] = ROcp10_917;
    sens->V[1] = VIcp10_117;
    sens->V[2] = VIcp10_217;
    sens->V[3] = VIcp10_317;
    sens->OM[1] = OMcp10_117;
    sens->OM[2] = OMcp10_217;
    sens->OM[3] = OMcp10_317;
    sens->J[1][1] = (1.0);
    sens->J[1][4] = -RLcp10_217;
    sens->J[1][5] = JTcp10_117_5;
    sens->J[1][6] = JTcp10_117_6;
    sens->J[2][2] = (1.0);
    sens->J[2][4] = RLcp10_117;
    sens->J[2][5] = JTcp10_217_5;
    sens->J[2][6] = JTcp10_217_6;
    sens->J[3][3] = (1.0);
    sens->J[3][5] = JTcp10_317_5;
    sens->J[3][6] = JTcp10_317_6;
    sens->J[4][5] = -S4;
    sens->J[4][6] = ROcp10_15;
    sens->J[4][17] = ROcp10_46;
    sens->J[5][5] = C4;
    sens->J[5][6] = ROcp10_25;
    sens->J[5][17] = ROcp10_56;
    sens->J[6][4] = (1.0);
    sens->J[6][6] = -S5;
    sens->J[6][17] = ROcp10_66;
    sens->A[1] = ACcp10_117;
    sens->A[2] = ACcp10_217;
    sens->A[3] = ACcp10_317;
    sens->OMP[1] = OPcp10_117;
    sens->OMP[2] = OPcp10_217;
    sens->OMP[3] = OPcp10_317;
 
// 
break;
case 12:
 


// = = Block_1_0_0_12_0_1 = = 
 
// Sensor Kinematics 


    ROcp11_15 = C4*C5;
    ROcp11_25 = S4*C5;
    ROcp11_75 = C4*S5;
    ROcp11_85 = S4*S5;
    ROcp11_46 = ROcp11_75*S6-S4*C6;
    ROcp11_56 = ROcp11_85*S6+C4*C6;
    ROcp11_66 = C5*S6;
    ROcp11_76 = ROcp11_75*C6+S4*S6;
    ROcp11_86 = ROcp11_85*C6-C4*S6;
    ROcp11_96 = C5*C6;
    OMcp11_15 = -qd[5]*S4;
    OMcp11_25 = qd[5]*C4;
    OMcp11_16 = OMcp11_15+ROcp11_15*qd[6];
    OMcp11_26 = OMcp11_25+ROcp11_25*qd[6];
    OMcp11_36 = qd[4]-qd[6]*S5;
    OPcp11_16 = ROcp11_15*qdd[6]-qdd[5]*S4-qd[4]*qd[5]*C4-qd[6]*(OMcp11_25*S5+ROcp11_25*qd[4]);
    OPcp11_26 = ROcp11_25*qdd[6]+qdd[5]*C4-qd[4]*qd[5]*S4+qd[6]*(OMcp11_15*S5+ROcp11_15*qd[4]);
    OPcp11_36 = qdd[4]-qdd[6]*S5-qd[5]*qd[6]*C5;

// = = Block_1_0_0_12_0_6 = = 
 
// Sensor Kinematics 


    ROcp11_116 = ROcp11_15*C16-ROcp11_76*S16;
    ROcp11_216 = ROcp11_25*C16-ROcp11_86*S16;
    ROcp11_316 = -(ROcp11_96*S16+C16*S5);
    ROcp11_716 = ROcp11_15*S16+ROcp11_76*C16;
    ROcp11_816 = ROcp11_25*S16+ROcp11_86*C16;
    ROcp11_916 = ROcp11_96*C16-S16*S5;
    RLcp11_116 = ROcp11_46*s->dpt[2][3];
    RLcp11_216 = ROcp11_56*s->dpt[2][3];
    RLcp11_316 = ROcp11_66*s->dpt[2][3];
    POcp11_116 = RLcp11_116+q[1];
    POcp11_216 = RLcp11_216+q[2];
    POcp11_316 = RLcp11_316+q[3];
    OMcp11_116 = OMcp11_16+ROcp11_46*qd[16];
    OMcp11_216 = OMcp11_26+ROcp11_56*qd[16];
    OMcp11_316 = OMcp11_36+ROcp11_66*qd[16];
    ORcp11_116 = OMcp11_26*RLcp11_316-OMcp11_36*RLcp11_216;
    ORcp11_216 = -(OMcp11_16*RLcp11_316-OMcp11_36*RLcp11_116);
    ORcp11_316 = OMcp11_16*RLcp11_216-OMcp11_26*RLcp11_116;
    VIcp11_116 = ORcp11_116+qd[1];
    VIcp11_216 = ORcp11_216+qd[2];
    VIcp11_316 = ORcp11_316+qd[3];
    OPcp11_116 = OPcp11_16+ROcp11_46*qdd[16]+qd[16]*(OMcp11_26*ROcp11_66-OMcp11_36*ROcp11_56);
    OPcp11_216 = OPcp11_26+ROcp11_56*qdd[16]-qd[16]*(OMcp11_16*ROcp11_66-OMcp11_36*ROcp11_46);
    OPcp11_316 = OPcp11_36+ROcp11_66*qdd[16]+qd[16]*(OMcp11_16*ROcp11_56-OMcp11_26*ROcp11_46);
    ACcp11_116 = qdd[1]+OMcp11_26*ORcp11_316-OMcp11_36*ORcp11_216+OPcp11_26*RLcp11_316-OPcp11_36*RLcp11_216;
    ACcp11_216 = qdd[2]-OMcp11_16*ORcp11_316+OMcp11_36*ORcp11_116-OPcp11_16*RLcp11_316+OPcp11_36*RLcp11_116;
    ACcp11_316 = qdd[3]+OMcp11_16*ORcp11_216-OMcp11_26*ORcp11_116+OPcp11_16*RLcp11_216-OPcp11_26*RLcp11_116;

// = = Block_1_0_0_12_1_0 = = 
 
// Symbolic Outputs  

    sens->P[1] = POcp11_116;
    sens->P[2] = POcp11_216;
    sens->P[3] = POcp11_316;
    sens->R[1][1] = ROcp11_116;
    sens->R[1][2] = ROcp11_216;
    sens->R[1][3] = ROcp11_316;
    sens->R[2][1] = ROcp11_46;
    sens->R[2][2] = ROcp11_56;
    sens->R[2][3] = ROcp11_66;
    sens->R[3][1] = ROcp11_716;
    sens->R[3][2] = ROcp11_816;
    sens->R[3][3] = ROcp11_916;
    sens->V[1] = VIcp11_116;
    sens->V[2] = VIcp11_216;
    sens->V[3] = VIcp11_316;
    sens->OM[1] = OMcp11_116;
    sens->OM[2] = OMcp11_216;
    sens->OM[3] = OMcp11_316;
    sens->A[1] = ACcp11_116;
    sens->A[2] = ACcp11_216;
    sens->A[3] = ACcp11_316;
    sens->OMP[1] = OPcp11_116;
    sens->OMP[2] = OPcp11_216;
    sens->OMP[3] = OPcp11_316;
 
// 
break;
case 13:
 


// = = Block_1_0_0_13_0_1 = = 
 
// Sensor Kinematics 


    ROcp12_15 = C4*C5;
    ROcp12_25 = S4*C5;
    ROcp12_75 = C4*S5;
    ROcp12_85 = S4*S5;
    ROcp12_46 = ROcp12_75*S6-S4*C6;
    ROcp12_56 = ROcp12_85*S6+C4*C6;
    ROcp12_66 = C5*S6;
    ROcp12_76 = ROcp12_75*C6+S4*S6;
    ROcp12_86 = ROcp12_85*C6-C4*S6;
    ROcp12_96 = C5*C6;
    OMcp12_15 = -qd[5]*S4;
    OMcp12_25 = qd[5]*C4;
    OMcp12_16 = OMcp12_15+ROcp12_15*qd[6];
    OMcp12_26 = OMcp12_25+ROcp12_25*qd[6];
    OMcp12_36 = qd[4]-qd[6]*S5;
    OPcp12_16 = ROcp12_15*qdd[6]-qdd[5]*S4-qd[4]*qd[5]*C4-qd[6]*(OMcp12_25*S5+ROcp12_25*qd[4]);
    OPcp12_26 = ROcp12_25*qdd[6]+qdd[5]*C4-qd[4]*qd[5]*S4+qd[6]*(OMcp12_15*S5+ROcp12_15*qd[4]);
    OPcp12_36 = qdd[4]-qdd[6]*S5-qd[5]*qd[6]*C5;

// = = Block_1_0_0_13_0_7 = = 
 
// Sensor Kinematics 


    ROcp12_117 = ROcp12_15*C17-ROcp12_76*S17;
    ROcp12_217 = ROcp12_25*C17-ROcp12_86*S17;
    ROcp12_317 = -(ROcp12_96*S17+C17*S5);
    ROcp12_717 = ROcp12_15*S17+ROcp12_76*C17;
    ROcp12_817 = ROcp12_25*S17+ROcp12_86*C17;
    ROcp12_917 = ROcp12_96*C17-S17*S5;
    RLcp12_117 = ROcp12_46*s->dpt[2][4];
    RLcp12_217 = ROcp12_56*s->dpt[2][4];
    RLcp12_317 = ROcp12_66*s->dpt[2][4];
    POcp12_117 = RLcp12_117+q[1];
    POcp12_217 = RLcp12_217+q[2];
    POcp12_317 = RLcp12_317+q[3];
    OMcp12_117 = OMcp12_16+ROcp12_46*qd[17];
    OMcp12_217 = OMcp12_26+ROcp12_56*qd[17];
    OMcp12_317 = OMcp12_36+ROcp12_66*qd[17];
    ORcp12_117 = OMcp12_26*RLcp12_317-OMcp12_36*RLcp12_217;
    ORcp12_217 = -(OMcp12_16*RLcp12_317-OMcp12_36*RLcp12_117);
    ORcp12_317 = OMcp12_16*RLcp12_217-OMcp12_26*RLcp12_117;
    VIcp12_117 = ORcp12_117+qd[1];
    VIcp12_217 = ORcp12_217+qd[2];
    VIcp12_317 = ORcp12_317+qd[3];
    OPcp12_117 = OPcp12_16+ROcp12_46*qdd[17]+qd[17]*(OMcp12_26*ROcp12_66-OMcp12_36*ROcp12_56);
    OPcp12_217 = OPcp12_26+ROcp12_56*qdd[17]-qd[17]*(OMcp12_16*ROcp12_66-OMcp12_36*ROcp12_46);
    OPcp12_317 = OPcp12_36+ROcp12_66*qdd[17]+qd[17]*(OMcp12_16*ROcp12_56-OMcp12_26*ROcp12_46);
    ACcp12_117 = qdd[1]+OMcp12_26*ORcp12_317-OMcp12_36*ORcp12_217+OPcp12_26*RLcp12_317-OPcp12_36*RLcp12_217;
    ACcp12_217 = qdd[2]-OMcp12_16*ORcp12_317+OMcp12_36*ORcp12_117-OPcp12_16*RLcp12_317+OPcp12_36*RLcp12_117;
    ACcp12_317 = qdd[3]+OMcp12_16*ORcp12_217-OMcp12_26*ORcp12_117+OPcp12_16*RLcp12_217-OPcp12_26*RLcp12_117;

// = = Block_1_0_0_13_1_0 = = 
 
// Symbolic Outputs  

    sens->P[1] = POcp12_117;
    sens->P[2] = POcp12_217;
    sens->P[3] = POcp12_317;
    sens->R[1][1] = ROcp12_117;
    sens->R[1][2] = ROcp12_217;
    sens->R[1][3] = ROcp12_317;
    sens->R[2][1] = ROcp12_46;
    sens->R[2][2] = ROcp12_56;
    sens->R[2][3] = ROcp12_66;
    sens->R[3][1] = ROcp12_717;
    sens->R[3][2] = ROcp12_817;
    sens->R[3][3] = ROcp12_917;
    sens->V[1] = VIcp12_117;
    sens->V[2] = VIcp12_217;
    sens->V[3] = VIcp12_317;
    sens->OM[1] = OMcp12_117;
    sens->OM[2] = OMcp12_217;
    sens->OM[3] = OMcp12_317;
    sens->A[1] = ACcp12_117;
    sens->A[2] = ACcp12_217;
    sens->A[3] = ACcp12_317;
    sens->OMP[1] = OPcp12_117;
    sens->OMP[2] = OPcp12_217;
    sens->OMP[3] = OPcp12_317;

break;
default:
break;
}


// ====== END Task 1 ====== 


}
 

