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
//	==> Function : F19 : External Forces
//	==> Flops complexity : 362
//
//	==> Generation Time :  0.000 seconds
//	==> Post-Processing :  0.010 seconds
//
//-------------------------------------------------------------
//
 
#include <math.h> 

#include "mbs_data.h"
#include "mbs_project_interface.h"
 
void mbs_extforces(double **frc,double **trq,
MbsData *s, double tsim)

// double frc[3][17];
// double trq[3][17];
{ 
double PxF1[4]; 
double RxF1[4][4]; 
double VxF1[4]; 
double OMxF1[4]; 
double AxF1[4]; 
double OMPxF1[4]; 
double *SWr1; 
double PxF2[4]; 
double RxF2[4][4]; 
double VxF2[4]; 
double OMxF2[4]; 
double AxF2[4]; 
double OMPxF2[4]; 
double *SWr2; 
 
#include "mbs_extforces_Frank_segway.h" 
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

// = = Block_0_0_0_0_0_6 = = 
 
// Trigonometric Variables  

  C16 = cos(q[16]);
  S16 = sin(q[16]);

// = = Block_0_0_0_0_0_7 = = 
 
// Trigonometric Variables  

  C17 = cos(q[17]);
  S17 = sin(q[17]);

// = = Block_0_0_1_1_0_1 = = 
 
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
  OMcp11_16 = OMcp11_15+qd[6]*ROcp11_15;
  OMcp11_26 = OMcp11_25+qd[6]*ROcp11_25;
  OMcp11_36 = qd[4]-qd[6]*S5;
  OPcp11_16 = -(qd[4]*qd[5]*C4+qd[6]*(qd[4]*ROcp11_25+OMcp11_25*S5)+qdd[5]*S4-qdd[6]*ROcp11_15);
  OPcp11_26 = -(qd[4]*qd[5]*S4-qd[6]*(qd[4]*ROcp11_15+OMcp11_15*S5)-qdd[5]*C4-qdd[6]*ROcp11_25);
  OPcp11_36 = qdd[4]-qd[5]*qd[6]*C5-qdd[6]*S5;

// = = Block_0_0_1_1_0_6 = = 
 
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
  ORcp11_116 = OMcp11_26*RLcp11_316-OMcp11_36*RLcp11_216;
  ORcp11_216 = -(OMcp11_16*RLcp11_316-OMcp11_36*RLcp11_116);
  ORcp11_316 = OMcp11_16*RLcp11_216-OMcp11_26*RLcp11_116;
  PxF1[1] = q[1]+RLcp11_116;
  PxF1[2] = q[2]+RLcp11_216;
  PxF1[3] = q[3]+RLcp11_316;
  RxF1[1][1] = ROcp11_116;
  RxF1[1][2] = ROcp11_216;
  RxF1[1][3] = ROcp11_316;
  RxF1[2][1] = ROcp11_46;
  RxF1[2][2] = ROcp11_56;
  RxF1[2][3] = ROcp11_66;
  RxF1[3][1] = ROcp11_716;
  RxF1[3][2] = ROcp11_816;
  RxF1[3][3] = ROcp11_916;
  VxF1[1] = qd[1]+ORcp11_116;
  VxF1[2] = qd[2]+ORcp11_216;
  VxF1[3] = qd[3]+ORcp11_316;
  OMxF1[1] = OMcp11_16+qd[16]*ROcp11_46;
  OMxF1[2] = OMcp11_26+qd[16]*ROcp11_56;
  OMxF1[3] = OMcp11_36+qd[16]*ROcp11_66;
  AxF1[1] = qdd[1]+OMcp11_26*ORcp11_316-OMcp11_36*ORcp11_216+OPcp11_26*RLcp11_316-OPcp11_36*RLcp11_216;
  AxF1[2] = qdd[2]-OMcp11_16*ORcp11_316+OMcp11_36*ORcp11_116-OPcp11_16*RLcp11_316+OPcp11_36*RLcp11_116;
  AxF1[3] = qdd[3]+OMcp11_16*ORcp11_216-OMcp11_26*ORcp11_116+OPcp11_16*RLcp11_216-OPcp11_26*RLcp11_116;
  OMPxF1[1] = OPcp11_16+qd[16]*(OMcp11_26*ROcp11_66-OMcp11_36*ROcp11_56)+qdd[16]*ROcp11_46;
  OMPxF1[2] = OPcp11_26-qd[16]*(OMcp11_16*ROcp11_66-OMcp11_36*ROcp11_46)+qdd[16]*ROcp11_56;
  OMPxF1[3] = OPcp11_36+qd[16]*(OMcp11_16*ROcp11_56-OMcp11_26*ROcp11_46)+qdd[16]*ROcp11_66;
 
// Sensor Forces Computation 

  SWr1 = user_ExtForces(PxF1,RxF1,VxF1,OMxF1,AxF1,OMPxF1,s,tsim,1);
 
// Sensor Dynamics : Forces projection on body-fixed frames 

  xfrc112 = ROcp11_116*SWr1[1]+ROcp11_216*SWr1[2]+ROcp11_316*SWr1[3];
  xfrc212 = ROcp11_46*SWr1[1]+ROcp11_56*SWr1[2]+ROcp11_66*SWr1[3];
  xfrc312 = ROcp11_716*SWr1[1]+ROcp11_816*SWr1[2]+ROcp11_916*SWr1[3];
  frc[1][16] = s->frc[1][16]+xfrc112;
  frc[2][16] = s->frc[2][16]+xfrc212;
  frc[3][16] = s->frc[3][16]+xfrc312;
  xtrq112 = ROcp11_116*SWr1[4]+ROcp11_216*SWr1[5]+ROcp11_316*SWr1[6];
  xtrq212 = ROcp11_46*SWr1[4]+ROcp11_56*SWr1[5]+ROcp11_66*SWr1[6];
  xtrq312 = ROcp11_716*SWr1[4]+ROcp11_816*SWr1[5]+ROcp11_916*SWr1[6];
  trq[1][16] = s->trq[1][16]+xtrq112-xfrc212*SWr1[9]+xfrc312*SWr1[8];
  trq[2][16] = s->trq[2][16]+xtrq212+xfrc112*SWr1[9]-xfrc312*SWr1[7];
  trq[3][16] = s->trq[3][16]+xtrq312-xfrc112*SWr1[8]+xfrc212*SWr1[7];

// = = Block_0_0_1_2_0_1 = = 
 
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
  OMcp12_16 = OMcp12_15+qd[6]*ROcp12_15;
  OMcp12_26 = OMcp12_25+qd[6]*ROcp12_25;
  OMcp12_36 = qd[4]-qd[6]*S5;
  OPcp12_16 = -(qd[4]*qd[5]*C4+qd[6]*(qd[4]*ROcp12_25+OMcp12_25*S5)+qdd[5]*S4-qdd[6]*ROcp12_15);
  OPcp12_26 = -(qd[4]*qd[5]*S4-qd[6]*(qd[4]*ROcp12_15+OMcp12_15*S5)-qdd[5]*C4-qdd[6]*ROcp12_25);
  OPcp12_36 = qdd[4]-qd[5]*qd[6]*C5-qdd[6]*S5;

// = = Block_0_0_1_2_0_7 = = 
 
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
  ORcp12_117 = OMcp12_26*RLcp12_317-OMcp12_36*RLcp12_217;
  ORcp12_217 = -(OMcp12_16*RLcp12_317-OMcp12_36*RLcp12_117);
  ORcp12_317 = OMcp12_16*RLcp12_217-OMcp12_26*RLcp12_117;
  PxF2[1] = q[1]+RLcp12_117;
  PxF2[2] = q[2]+RLcp12_217;
  PxF2[3] = q[3]+RLcp12_317;
  RxF2[1][1] = ROcp12_117;
  RxF2[1][2] = ROcp12_217;
  RxF2[1][3] = ROcp12_317;
  RxF2[2][1] = ROcp12_46;
  RxF2[2][2] = ROcp12_56;
  RxF2[2][3] = ROcp12_66;
  RxF2[3][1] = ROcp12_717;
  RxF2[3][2] = ROcp12_817;
  RxF2[3][3] = ROcp12_917;
  VxF2[1] = qd[1]+ORcp12_117;
  VxF2[2] = qd[2]+ORcp12_217;
  VxF2[3] = qd[3]+ORcp12_317;
  OMxF2[1] = OMcp12_16+qd[17]*ROcp12_46;
  OMxF2[2] = OMcp12_26+qd[17]*ROcp12_56;
  OMxF2[3] = OMcp12_36+qd[17]*ROcp12_66;
  AxF2[1] = qdd[1]+OMcp12_26*ORcp12_317-OMcp12_36*ORcp12_217+OPcp12_26*RLcp12_317-OPcp12_36*RLcp12_217;
  AxF2[2] = qdd[2]-OMcp12_16*ORcp12_317+OMcp12_36*ORcp12_117-OPcp12_16*RLcp12_317+OPcp12_36*RLcp12_117;
  AxF2[3] = qdd[3]+OMcp12_16*ORcp12_217-OMcp12_26*ORcp12_117+OPcp12_16*RLcp12_217-OPcp12_26*RLcp12_117;
  OMPxF2[1] = OPcp12_16+qd[17]*(OMcp12_26*ROcp12_66-OMcp12_36*ROcp12_56)+qdd[17]*ROcp12_46;
  OMPxF2[2] = OPcp12_26-qd[17]*(OMcp12_16*ROcp12_66-OMcp12_36*ROcp12_46)+qdd[17]*ROcp12_56;
  OMPxF2[3] = OPcp12_36+qd[17]*(OMcp12_16*ROcp12_56-OMcp12_26*ROcp12_46)+qdd[17]*ROcp12_66;
 
// Sensor Forces Computation 

  SWr2 = user_ExtForces(PxF2,RxF2,VxF2,OMxF2,AxF2,OMPxF2,s,tsim,2);
 
// Sensor Dynamics : Forces projection on body-fixed frames 

  xfrc113 = ROcp12_117*SWr2[1]+ROcp12_217*SWr2[2]+ROcp12_317*SWr2[3];
  xfrc213 = ROcp12_46*SWr2[1]+ROcp12_56*SWr2[2]+ROcp12_66*SWr2[3];
  xfrc313 = ROcp12_717*SWr2[1]+ROcp12_817*SWr2[2]+ROcp12_917*SWr2[3];
  frc[1][17] = s->frc[1][17]+xfrc113;
  frc[2][17] = s->frc[2][17]+xfrc213;
  frc[3][17] = s->frc[3][17]+xfrc313;
  xtrq113 = ROcp12_117*SWr2[4]+ROcp12_217*SWr2[5]+ROcp12_317*SWr2[6];
  xtrq213 = ROcp12_46*SWr2[4]+ROcp12_56*SWr2[5]+ROcp12_66*SWr2[6];
  xtrq313 = ROcp12_717*SWr2[4]+ROcp12_817*SWr2[5]+ROcp12_917*SWr2[6];
  trq[1][17] = s->trq[1][17]+xtrq113-xfrc213*SWr2[9]+xfrc313*SWr2[8];
  trq[2][17] = s->trq[2][17]+xtrq213+xfrc113*SWr2[9]-xfrc313*SWr2[7];
  trq[3][17] = s->trq[3][17]+xtrq313-xfrc113*SWr2[8]+xfrc213*SWr2[7];

// = = Block_0_0_1_2_1_0 = = 
 
// Symbolic Outputs  

  frc[1][6] = s->frc[1][6];
  frc[2][6] = s->frc[2][6];
  frc[3][6] = s->frc[3][6];
  frc[1][7] = s->frc[1][7];
  frc[2][7] = s->frc[2][7];
  frc[3][7] = s->frc[3][7];
  frc[1][8] = s->frc[1][8];
  frc[2][8] = s->frc[2][8];
  frc[3][8] = s->frc[3][8];
  frc[1][9] = s->frc[1][9];
  frc[2][9] = s->frc[2][9];
  frc[3][9] = s->frc[3][9];
  frc[1][10] = s->frc[1][10];
  frc[2][10] = s->frc[2][10];
  frc[3][10] = s->frc[3][10];
  frc[1][11] = s->frc[1][11];
  frc[2][11] = s->frc[2][11];
  frc[3][11] = s->frc[3][11];
  frc[1][12] = s->frc[1][12];
  frc[2][12] = s->frc[2][12];
  frc[3][12] = s->frc[3][12];
  frc[1][13] = s->frc[1][13];
  frc[2][13] = s->frc[2][13];
  frc[3][13] = s->frc[3][13];
  frc[1][14] = s->frc[1][14];
  frc[2][14] = s->frc[2][14];
  frc[3][14] = s->frc[3][14];
  frc[1][15] = s->frc[1][15];
  frc[2][15] = s->frc[2][15];
  frc[3][15] = s->frc[3][15];
  trq[1][6] = s->trq[1][6];
  trq[2][6] = s->trq[2][6];
  trq[3][6] = s->trq[3][6];
  trq[1][7] = s->trq[1][7];
  trq[2][7] = s->trq[2][7];
  trq[3][7] = s->trq[3][7];
  trq[1][8] = s->trq[1][8];
  trq[2][8] = s->trq[2][8];
  trq[3][8] = s->trq[3][8];
  trq[1][9] = s->trq[1][9];
  trq[2][9] = s->trq[2][9];
  trq[3][9] = s->trq[3][9];
  trq[1][10] = s->trq[1][10];
  trq[2][10] = s->trq[2][10];
  trq[3][10] = s->trq[3][10];
  trq[1][11] = s->trq[1][11];
  trq[2][11] = s->trq[2][11];
  trq[3][11] = s->trq[3][11];
  trq[1][12] = s->trq[1][12];
  trq[2][12] = s->trq[2][12];
  trq[3][12] = s->trq[3][12];
  trq[1][13] = s->trq[1][13];
  trq[2][13] = s->trq[2][13];
  trq[3][13] = s->trq[3][13];
  trq[1][14] = s->trq[1][14];
  trq[2][14] = s->trq[2][14];
  trq[3][14] = s->trq[3][14];
  trq[1][15] = s->trq[1][15];
  trq[2][15] = s->trq[2][15];
  trq[3][15] = s->trq[3][15];

// ====== END Task 0 ====== 


}
 

