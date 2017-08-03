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
//	==> Function : F 2 : Inverse Dynamics : RNEA
//	==> Flops complexity : 1348
//
//	==> Generation Time :  0.010 seconds
//	==> Post-Processing :  0.020 seconds
//
//-------------------------------------------------------------
//
 
#include <math.h> 

#include "mbs_data.h"
#include "mbs_project_interface.h"
 
void mbs_invdyna(double *Qq,
MbsData *s, double tsim)

// double Qq[17];
{ 
 
#include "mbs_invdyna_Frank_segway.h" 
#define q s->q 
#define qd s->qd 
#define qdd s->qdd 
 
 

// === begin imp_aux === 

// === end imp_aux === 

// ===== BEGIN task 0 ===== 

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
  C15 = cos(q[15]);
  S15 = sin(q[15]);

// = = Block_0_0_0_0_0_6 = = 
 
// Trigonometric Variables  

  C16 = cos(q[16]);
  S16 = sin(q[16]);

// = = Block_0_0_0_0_0_7 = = 
 
// Trigonometric Variables  

  C17 = cos(q[17]);
  S17 = sin(q[17]);

// = = Block_0_1_0_0_0_0 = = 
 
// Forward Kinematics 

  ALPHA33 = qdd[3]-s->g[3];
  ALPHA14 = qdd[1]*C4+qdd[2]*S4;
  ALPHA24 = -(qdd[1]*S4-qdd[2]*C4);
  OM35 = qd[4]*C5;
  OMp35 = -(qd[4]*qd[5]*S5-qdd[4]*C5);
  ALPHA15 = ALPHA14*C5-ALPHA33*S5;
  ALPHA35 = ALPHA14*S5+ALPHA33*C5;
  OM16 = qd[6]-qd[4]*S5;
  OM26 = qd[5]*C6+OM35*S6;
  OM36 = -(qd[5]*S6-OM35*C6);
  OMp16 = qdd[6]-qd[4]*qd[5]*C5-qdd[4]*S5;
  OMp26 = C6*(qdd[5]+qd[6]*OM35)+S6*(OMp35-qd[5]*qd[6]);
  OMp36 = C6*(OMp35-qd[5]*qd[6])-S6*(qdd[5]+qd[6]*OM35);
  BS26 = OM16*OM26;
  BS36 = OM16*OM36;
  BS56 = -(OM16*OM16+OM36*OM36);
  BS66 = OM26*OM36;
  BETA26 = BS26-OMp36;
  BETA86 = BS66+OMp16;
  ALPHA26 = ALPHA24*C6+ALPHA35*S6;
  ALPHA36 = -(ALPHA24*S6-ALPHA35*C6);
  BS57 = -(OM16*OM16+OM36*OM36);
  BS67 = OM26*OM36;
  BETA27 = OM16*OM26-OMp36;
  BETA37 = OMp26+OM16*OM36;
  BETA67 = BS67-OMp16;
  BETA87 = BS67+OMp16;
  ALPHA17 = ALPHA15+(2.0)*qd[7]*OM26+Dz73*(BS36+OMp26);
  ALPHA27 = ALPHA26-(2.0)*qd[7]*OM16+Dz73*(BS66-OMp16);
  ALPHA37 = qdd[7]+ALPHA36-Dz73*(OM16*OM16+OM26*OM26);
  OM18 = OM16*C8-OM36*S8;
  OM28 = qd[8]+OM26;
  OM38 = OM16*S8+OM36*C8;
  OMp18 = C8*(OMp16-qd[8]*OM36)-S8*(OMp36+qd[8]*OM16);
  OMp28 = qdd[8]+OMp26;
  OMp38 = C8*(OMp36+qd[8]*OM16)+S8*(OMp16-qd[8]*OM36);
  BS28 = OM18*OM28;
  ALPHA18 = C8*(ALPHA17+BETA27*s->dpt[2][6])-S8*(ALPHA37+BETA87*s->dpt[2][6]);
  ALPHA28 = ALPHA27+BS57*s->dpt[2][6];
  ALPHA38 = C8*(ALPHA37+BETA87*s->dpt[2][6])+S8*(ALPHA17+BETA27*s->dpt[2][6]);
  OM19 = qd[9]+OM18;
  OM29 = OM28*C9+OM38*S9;
  OM39 = -(OM28*S9-OM38*C9);
  OMp19 = qdd[9]+OMp18;
  OMp29 = C9*(OMp28+qd[9]*OM38)+S9*(OMp38-qd[9]*OM28);
  OMp39 = C9*(OMp38-qd[9]*OM28)-S9*(OMp28+qd[9]*OM38);
  BS39 = OM19*OM39;
  BS99 = -(OM19*OM19+OM29*OM29);
  BETA39 = BS39+OMp29;
  BETA69 = OM29*OM39-OMp19;
  ALPHA29 = ALPHA28*C9+ALPHA38*S9;
  ALPHA39 = -(ALPHA28*S9-ALPHA38*C9);
  OM110 = OM19*C10+OM29*S10;
  OM210 = -(OM19*S10-OM29*C10);
  OM310 = qd[10]+OM39;
  OMp110 = C10*(OMp19+qd[10]*OM29)+S10*(OMp29-qd[10]*OM19);
  OMp210 = C10*(OMp29-qd[10]*OM19)-S10*(OMp19+qd[10]*OM29);
  OMp310 = qdd[10]+OMp39;
  BS610 = OM210*OM310;
  OM111 = OM16*C11-OM36*S11;
  OM211 = qd[11]+OM26;
  OM311 = OM16*S11+OM36*C11;
  OMp111 = C11*(OMp16-qd[11]*OM36)-S11*(OMp36+qd[11]*OM16);
  OMp211 = qdd[11]+OMp26;
  OMp311 = C11*(OMp36+qd[11]*OM16)+S11*(OMp16-qd[11]*OM36);
  BS211 = OM111*OM211;
  ALPHA111 = C11*(ALPHA17+BETA27*s->dpt[2][7])-S11*(ALPHA37+BETA87*s->dpt[2][7]);
  ALPHA211 = ALPHA27+BS57*s->dpt[2][7];
  ALPHA311 = C11*(ALPHA37+BETA87*s->dpt[2][7])+S11*(ALPHA17+BETA27*s->dpt[2][7]);
  OM112 = qd[12]+OM111;
  OM212 = OM211*C12+OM311*S12;
  OM312 = -(OM211*S12-OM311*C12);
  OMp112 = qdd[12]+OMp111;
  OMp212 = C12*(OMp211+qd[12]*OM311)+S12*(OMp311-qd[12]*OM211);
  OMp312 = C12*(OMp311-qd[12]*OM211)-S12*(OMp211+qd[12]*OM311);
  BS312 = OM112*OM312;
  BS912 = -(OM112*OM112+OM212*OM212);
  BETA312 = BS312+OMp212;
  BETA612 = OM212*OM312-OMp112;
  ALPHA212 = ALPHA211*C12+ALPHA311*S12;
  ALPHA312 = -(ALPHA211*S12-ALPHA311*C12);
  OM113 = OM112*C13+OM212*S13;
  OM213 = -(OM112*S13-OM212*C13);
  OM313 = qd[13]+OM312;
  OMp113 = C13*(OMp112+qd[13]*OM212)+S13*(OMp212-qd[13]*OM112);
  OMp213 = C13*(OMp212-qd[13]*OM112)-S13*(OMp112+qd[13]*OM212);
  OMp313 = qdd[13]+OMp312;
  BS613 = OM213*OM313;
  OM114 = OM16*C14+OM26*S14;
  OM214 = -(OM16*S14-OM26*C14);
  OM314 = qd[14]+OM36;
  OMp114 = C14*(OMp16+qd[14]*OM26)+S14*(OMp26-qd[14]*OM16);
  OMp214 = C14*(OMp26-qd[14]*OM16)-S14*(OMp16+qd[14]*OM26);
  OMp314 = qdd[14]+OMp36;
  ALPHA114 = C14*(ALPHA17+BETA37*s->dpt[3][8])+S14*(ALPHA27+BETA67*s->dpt[3][8]);
  ALPHA214 = C14*(ALPHA27+BETA67*s->dpt[3][8])-S14*(ALPHA17+BETA37*s->dpt[3][8]);
  ALPHA314 = ALPHA37-s->dpt[3][8]*(OM16*OM16+OM26*OM26);
  OM115 = OM114*C15-OM314*S15;
  OM215 = qd[15]+OM214;
  OM315 = OM114*S15+OM314*C15;
  OM116 = OM16*C16-OM36*S16;
  OM216 = qd[16]+OM26;
  OM316 = OM16*S16+OM36*C16;
  OM117 = OM16*C17-OM36*S17;
  OM217 = qd[17]+OM26;
  OM317 = OM16*S17+OM36*C17;
 
// Backward Dynamics 

  Fs117 = -(s->frc[1][17]-s->m[17]*(C17*(ALPHA15+BETA26*s->dpt[2][4])-S17*(ALPHA36+BETA86*s->dpt[2][4])));
  Fs317 = -(s->frc[3][17]-s->m[17]*(C17*(ALPHA36+BETA86*s->dpt[2][4])+S17*(ALPHA15+BETA26*s->dpt[2][4])));
  Cq117 = -(s->trq[1][17]-s->In[1][17]*(C17*(OMp16-qd[17]*OM36)-S17*(OMp36+qd[17]*OM16))+OM217*OM317*(s->In[5][17]-
 s->In[9][17]));
  Cq217 = -(s->trq[2][17]-s->In[5][17]*(qdd[17]+OMp26)-OM117*OM317*(s->In[1][17]-s->In[9][17]));
  Cq317 = -(s->trq[3][17]-s->In[9][17]*(C17*(OMp36+qd[17]*OM16)+S17*(OMp16-qd[17]*OM36))+OM117*OM217*(s->In[1][17]-
 s->In[5][17]));
  Fs116 = -(s->frc[1][16]-s->m[16]*(C16*(ALPHA15+BETA26*s->dpt[2][3])-S16*(ALPHA36+BETA86*s->dpt[2][3])));
  Fs316 = -(s->frc[3][16]-s->m[16]*(C16*(ALPHA36+BETA86*s->dpt[2][3])+S16*(ALPHA15+BETA26*s->dpt[2][3])));
  Cq116 = -(s->trq[1][16]-s->In[1][16]*(C16*(OMp16-qd[16]*OM36)-S16*(OMp36+qd[16]*OM16))+OM216*OM316*(s->In[5][16]-
 s->In[9][16]));
  Cq216 = -(s->trq[2][16]-s->In[5][16]*(qdd[16]+OMp26)-OM116*OM316*(s->In[1][16]-s->In[9][16]));
  Cq316 = -(s->trq[3][16]-s->In[9][16]*(C16*(OMp36+qd[16]*OM16)+S16*(OMp16-qd[16]*OM36))+OM116*OM216*(s->In[1][16]-
 s->In[5][16]));
  Fs115 = -(s->frc[1][15]-s->m[15]*(ALPHA114*C15-ALPHA314*S15));
  Fs315 = -(s->frc[3][15]-s->m[15]*(ALPHA114*S15+ALPHA314*C15));
  Cq115 = -(s->trq[1][15]-s->In[1][15]*(C15*(OMp114-qd[15]*OM314)-S15*(OMp314+qd[15]*OM114))+OM215*OM315*(s->In[5][15]-
 s->In[9][15]));
  Cq215 = -(s->trq[2][15]-s->In[5][15]*(qdd[15]+OMp214)-OM115*OM315*(s->In[1][15]-s->In[9][15]));
  Cq315 = -(s->trq[3][15]-s->In[9][15]*(C15*(OMp314+qd[15]*OM114)+S15*(OMp114-qd[15]*OM314))+OM115*OM215*(s->In[1][15]-
 s->In[5][15]));
  Fq114 = -(s->frc[1][14]-s->m[14]*ALPHA114-Fs115*C15-Fs315*S15);
  Fq214 = -(s->frc[2][14]+s->frc[2][15]-s->m[14]*ALPHA214-s->m[15]*ALPHA214);
  Cq114 = -(s->trq[1][14]-s->In[1][14]*OMp114-Cq115*C15-Cq315*S15+OM214*OM314*(s->In[5][14]-s->In[9][14]));
  Cq214 = -(s->trq[2][14]-Cq215-s->In[5][14]*OMp214-OM114*OM314*(s->In[1][14]-s->In[9][14]));
  Cq314 = -(s->trq[3][14]-s->In[9][14]*OMp314+Cq115*S15-Cq315*C15+OM114*OM214*(s->In[1][14]-s->In[5][14]));
  Fs113 = -(s->frc[1][13]+s->m[13]*(s->l[2][13]*(OMp313-OM113*OM213)-s->l[3][13]*(OMp213+OM113*OM313)-C13*(ALPHA111+
 BETA312*s->dpt[3][15])-S13*(ALPHA212+BETA612*s->dpt[3][15])));
  Fs213 = -(s->frc[2][13]+s->m[13]*(s->l[2][13]*(OM113*OM113+OM313*OM313)-s->l[3][13]*(BS613-OMp113)-C13*(ALPHA212+
 BETA612*s->dpt[3][15])+S13*(ALPHA111+BETA312*s->dpt[3][15])));
  Fs313 = -(s->frc[3][13]-s->m[13]*(ALPHA312+BS912*s->dpt[3][15]+s->l[2][13]*(BS613+OMp113)-s->l[3][13]*(OM113*OM113+
 OM213*OM213)));
  Cq113 = -(s->trq[1][13]-s->In[1][13]*OMp113+Fs213*s->l[3][13]-Fs313*s->l[2][13]+OM213*OM313*(s->In[5][13]-s->In[9][13]
 ));
  Cq213 = -(s->trq[2][13]-s->In[5][13]*OMp213-Fs113*s->l[3][13]-OM113*OM313*(s->In[1][13]-s->In[9][13]));
  Cq313 = -(s->trq[3][13]-s->In[9][13]*OMp313+Fs113*s->l[2][13]+OM113*OM213*(s->In[1][13]-s->In[5][13]));
  Fs112 = -(s->frc[1][12]-s->m[12]*(ALPHA111+BETA312*s->l[3][12]-s->l[1][12]*(OM212*OM212+OM312*OM312)));
  Fs212 = -(s->frc[2][12]-s->m[12]*(ALPHA212+BETA612*s->l[3][12]+s->l[1][12]*(OMp312+OM112*OM212)));
  Fs312 = -(s->frc[3][12]-s->m[12]*(ALPHA312+BS912*s->l[3][12]+s->l[1][12]*(BS312-OMp212)));
  Fq212 = Fs212+Fs113*S13+Fs213*C13;
  Fq312 = Fs312+Fs313;
  Cq112 = -(s->trq[1][12]-s->In[1][12]*OMp112-Cq113*C13+Cq213*S13+Fs212*s->l[3][12]+OM212*OM312*(s->In[5][12]-
 s->In[9][12])+s->dpt[3][15]*(Fs113*S13+Fs213*C13));
  Cq212 = -(s->trq[2][12]-s->In[5][12]*OMp212-Cq113*S13-Cq213*C13-Fs112*s->l[3][12]+Fs312*s->l[1][12]-OM112*OM312*(
 s->In[1][12]-s->In[9][12])-s->dpt[3][15]*(Fs113*C13-Fs213*S13));
  Cq312 = -(s->trq[3][12]-Cq313-s->In[9][12]*OMp312-Fs212*s->l[1][12]+OM112*OM212*(s->In[1][12]-s->In[5][12]));
  Fs111 = -(s->frc[1][11]-s->m[11]*(ALPHA111-s->l[1][11]*(OM211*OM211+OM311*OM311)+s->l[2][11]*(BS211-OMp311)));
  Fs211 = -(s->frc[2][11]-s->m[11]*(ALPHA211+s->l[1][11]*(BS211+OMp311)-s->l[2][11]*(OM111*OM111+OM311*OM311)));
  Fs311 = -(s->frc[3][11]-s->m[11]*(ALPHA311-s->l[1][11]*(OMp211-OM111*OM311)+s->l[2][11]*(OMp111+OM211*OM311)));
  Fq111 = Fs111+Fs112+Fs113*C13-Fs213*S13;
  Fq311 = Fs311+Fq212*S12+Fq312*C12;
  Cq111 = -(s->trq[1][11]-Cq112-s->In[1][11]*OMp111-Fs311*s->l[2][11]+OM211*OM311*(s->In[5][11]-s->In[9][11]));
  Cq211 = -(s->trq[2][11]-s->In[5][11]*OMp211-Cq212*C12+Cq312*S12+Fs311*s->l[1][11]-OM111*OM311*(s->In[1][11]-
 s->In[9][11]));
  Cq311 = -(s->trq[3][11]-s->In[9][11]*OMp311-Cq212*S12-Cq312*C12+Fs111*s->l[2][11]-Fs211*s->l[1][11]+OM111*OM211*(
 s->In[1][11]-s->In[5][11]));
  Fs110 = -(s->frc[1][10]+s->m[10]*(s->l[2][10]*(OMp310-OM110*OM210)-s->l[3][10]*(OMp210+OM110*OM310)-C10*(ALPHA18+
 BETA39*s->dpt[3][11])-S10*(ALPHA29+BETA69*s->dpt[3][11])));
  Fs210 = -(s->frc[2][10]+s->m[10]*(s->l[2][10]*(OM110*OM110+OM310*OM310)-s->l[3][10]*(BS610-OMp110)-C10*(ALPHA29+BETA69
 *s->dpt[3][11])+S10*(ALPHA18+BETA39*s->dpt[3][11])));
  Fs310 = -(s->frc[3][10]-s->m[10]*(ALPHA39+BS99*s->dpt[3][11]+s->l[2][10]*(BS610+OMp110)-s->l[3][10]*(OM110*OM110+OM210
 *OM210)));
  Cq110 = -(s->trq[1][10]-s->In[1][10]*OMp110+Fs210*s->l[3][10]-Fs310*s->l[2][10]+OM210*OM310*(s->In[5][10]-s->In[9][10]
 ));
  Cq210 = -(s->trq[2][10]-s->In[5][10]*OMp210-Fs110*s->l[3][10]-OM110*OM310*(s->In[1][10]-s->In[9][10]));
  Cq310 = -(s->trq[3][10]-s->In[9][10]*OMp310+Fs110*s->l[2][10]+OM110*OM210*(s->In[1][10]-s->In[5][10]));
  Fs19 = -(s->frc[1][9]-s->m[9]*(ALPHA18+BETA39*s->l[3][9]-s->l[1][9]*(OM29*OM29+OM39*OM39)));
  Fs29 = -(s->frc[2][9]-s->m[9]*(ALPHA29+BETA69*s->l[3][9]+s->l[1][9]*(OMp39+OM19*OM29)));
  Fs39 = -(s->frc[3][9]-s->m[9]*(ALPHA39+BS99*s->l[3][9]+s->l[1][9]*(BS39-OMp29)));
  Fq29 = Fs29+Fs110*S10+Fs210*C10;
  Fq39 = Fs310+Fs39;
  Cq19 = -(s->trq[1][9]-s->In[1][9]*OMp19-Cq110*C10+Cq210*S10+Fs29*s->l[3][9]+OM29*OM39*(s->In[5][9]-s->In[9][9])+
 s->dpt[3][11]*(Fs110*S10+Fs210*C10));
  Cq29 = -(s->trq[2][9]-s->In[5][9]*OMp29-Cq110*S10-Cq210*C10-Fs19*s->l[3][9]+Fs39*s->l[1][9]-OM19*OM39*(s->In[1][9]-
 s->In[9][9])-s->dpt[3][11]*(Fs110*C10-Fs210*S10));
  Cq39 = -(s->trq[3][9]-Cq310-s->In[9][9]*OMp39-Fs29*s->l[1][9]+OM19*OM29*(s->In[1][9]-s->In[5][9]));
  Fs18 = -(s->frc[1][8]-s->m[8]*(ALPHA18-s->l[1][8]*(OM28*OM28+OM38*OM38)+s->l[2][8]*(BS28-OMp38)));
  Fs28 = -(s->frc[2][8]-s->m[8]*(ALPHA28+s->l[1][8]*(BS28+OMp38)-s->l[2][8]*(OM18*OM18+OM38*OM38)));
  Fs38 = -(s->frc[3][8]-s->m[8]*(ALPHA38-s->l[1][8]*(OMp28-OM18*OM38)+s->l[2][8]*(OMp18+OM28*OM38)));
  Fq18 = Fs18+Fs19+Fs110*C10-Fs210*S10;
  Fq38 = Fs38+Fq29*S9+Fq39*C9;
  Cq18 = -(s->trq[1][8]-Cq19-s->In[1][8]*OMp18-Fs38*s->l[2][8]+OM28*OM38*(s->In[5][8]-s->In[9][8]));
  Cq28 = -(s->trq[2][8]-s->In[5][8]*OMp28-Cq29*C9+Cq39*S9+Fs38*s->l[1][8]-OM18*OM38*(s->In[1][8]-s->In[9][8]));
  Cq38 = -(s->trq[3][8]-s->In[9][8]*OMp38-Cq29*S9-Cq39*C9+Fs18*s->l[2][8]-Fs28*s->l[1][8]+OM18*OM28*(s->In[1][8]-
 s->In[5][8]));
  Fq17 = -(s->frc[1][7]-s->m[7]*ALPHA17-Fq111*C11-Fq114*C14-Fq18*C8+Fq214*S14-Fq311*S11-Fq38*S8);
  Fq27 = -(s->frc[2][7]-Fs211-Fs28-s->m[7]*ALPHA27-Fq114*S14-Fq212*C12-Fq214*C14-Fq29*C9+Fq312*S12+Fq39*S9);
  Fq37 = -(s->frc[3][14]+s->frc[3][7]-s->m[14]*ALPHA314-s->m[7]*ALPHA37+Fq111*S11+Fq18*S8-Fq311*C11-Fq38*C8+Fs115*S15-
 Fs315*C15);
  Fs16 = -(s->frc[1][6]-s->m[6]*(ALPHA15+BETA26*s->l[2][6]-s->l[1][6]*(OM26*OM26+OM36*OM36)));
  Fs26 = -(s->frc[2][6]-s->m[6]*(ALPHA26+BS56*s->l[2][6]+s->l[1][6]*(BS26+OMp36)));
  Fs36 = -(s->frc[3][6]-s->m[6]*(ALPHA36+BETA86*s->l[2][6]+s->l[1][6]*(BS36-OMp26)));
  Fq16 = Fq17+Fs16+Fs116*C16+Fs117*C17+Fs316*S16+Fs317*S17;
  Fq26 = -(s->frc[2][16]+s->frc[2][17]-Fq27-Fs26-s->m[16]*(ALPHA26+BS56*s->dpt[2][3])-s->m[17]*(ALPHA26+BS56*
 s->dpt[2][4]));
  Fq36 = Fq37+Fs36-Fs116*S16-Fs117*S17+Fs316*C16+Fs317*C17;
  Cq16 = -(s->trq[1][6]+s->trq[1][7]-s->In[1][6]*OMp16-s->In[1][7]*OMp16-Cq111*C11-Cq114*C14-Cq116*C16-Cq117*C17-Cq18*C8
 +Cq214*S14-Cq311*S11-Cq316*S16-Cq317*S17-Cq38*S8+Dz73*Fq27-Fs36*s->l[2][6]+OM26*OM36*(s->In[5][6]-s->In[9][6])+OM26*OM36*(
 s->In[5][7]-s->In[9][7])+s->dpt[2][3]*(Fs116*S16-Fs316*C16)+s->dpt[2][4]*(Fs117*S17-Fs317*C17)+s->dpt[2][6]*(Fq18*S8-Fq38*C8
 )+s->dpt[2][7]*(Fq111*S11-Fq311*C11)+s->dpt[3][8]*(Fq114*S14+Fq214*C14));
  Cq26 = -(s->trq[2][6]+s->trq[2][7]-Cq211-Cq216-Cq217-Cq28-s->In[5][6]*OMp26-s->In[5][7]*OMp26-Cq114*S14-Cq214*C14-Dz73
 *Fq17+Fs36*s->l[1][6]-OM16*OM36*(s->In[1][6]-s->In[9][6])-OM16*OM36*(s->In[1][7]-s->In[9][7])-s->dpt[3][8]*(Fq114*C14-Fq214*
 S14));
  Cq36 = -(s->trq[3][6]+s->trq[3][7]-Cq314-s->In[9][6]*OMp36-s->In[9][7]*OMp36+Cq111*S11+Cq116*S16+Cq117*S17+Cq18*S8-
 Cq311*C11-Cq316*C16-Cq317*C17-Cq38*C8+Fs16*s->l[2][6]-Fs26*s->l[1][6]+OM16*OM26*(s->In[1][6]-s->In[5][6])+OM16*OM26*(
 s->In[1][7]-s->In[5][7])+s->dpt[2][3]*(Fs116*C16+Fs316*S16)+s->dpt[2][4]*(Fs117*C17+Fs317*S17)+s->dpt[2][6]*(Fq18*C8+Fq38*S8
 )+s->dpt[2][7]*(Fq111*C11+Fq311*S11));
  Fq25 = Fq26*C6-Fq36*S6;
  Fq35 = Fq26*S6+Fq36*C6;
  Cq25 = Cq26*C6-Cq36*S6;
  Fq14 = Fq16*C5+Fq35*S5;
  Fq34 = -(Fq16*S5-Fq35*C5);
  Cq34 = -(Cq16*S5-C5*(Cq26*S6+Cq36*C6));
  Fq13 = Fq14*C4-Fq25*S4;
  Fq23 = Fq14*S4+Fq25*C4;

// = = Block_0_2_0_0_0_0 = = 
 
// Symbolic Outputs  

  Qq[1] = Fq13;
  Qq[2] = Fq23;
  Qq[3] = Fq34;
  Qq[4] = Cq34;
  Qq[5] = Cq25;
  Qq[6] = Cq16;
  Qq[7] = Fq37;
  Qq[8] = Cq28;
  Qq[9] = Cq19;
  Qq[10] = Cq310;
  Qq[11] = Cq211;
  Qq[12] = Cq112;
  Qq[13] = Cq313;
  Qq[14] = Cq314;
  Qq[15] = Cq215;
  Qq[16] = Cq216;
  Qq[17] = Cq217;

// ====== END Task 0 ====== 


}
 

