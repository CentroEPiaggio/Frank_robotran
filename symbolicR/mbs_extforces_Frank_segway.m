%
%-------------------------------------------------------------
%
%	ROBOTRAN - Version 6.6 (build : february 22, 2008)
%
%	Copyright 
%	Universite catholique de Louvain 
%	Departement de Mecanique 
%	Unite de Production Mecanique et Machines 
%	2, Place du Levant 
%	1348 Louvain-la-Neuve 
%	http://www.robotran.be// 
%
%	==> Generation Date : Sat Jul 15 16:15:34 2017
%
%	==> Project name : Frank_segway
%	==> using XML input file 
%
%	==> Number of joints : 17
%
%	==> Function : F19 : External Forces
%	==> Flops complexity : 800
%
%	==> Generation Time :  0.010 seconds
%	==> Post-Processing :  0.020 seconds
%
%-------------------------------------------------------------
%
function [frc,trq] = extforces(s,tsim,usrfun)

 frc = zeros(3,17);
 trq = zeros(3,17);

q = s.q; 
qd = s.qd; 
qdd = s.qdd; 
frc = s.frc; 
trq = s.trq; 

% === begin imp_aux === 

% === end imp_aux === 

% ===== BEGIN task 0 ===== 
 
% Sensor Kinematics 



% = = Block_0_0_0_0_0_1 = = 
 
% Trigonometric Variables  

  C4 = cos(q(4));
  S4 = sin(q(4));
  C5 = cos(q(5));
  S5 = sin(q(5));
  C6 = cos(q(6));
  S6 = sin(q(6));

% = = Block_0_0_0_0_0_6 = = 
 
% Trigonometric Variables  

  C16 = cos(q(16));
  S16 = sin(q(16));

% = = Block_0_0_0_0_0_7 = = 
 
% Trigonometric Variables  

  C17 = cos(q(17));
  S17 = sin(q(17));

% = = Block_0_0_1_1_0_1 = = 
 
% Sensor Kinematics 


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
  OMcp11_15 = -qd(5)*S4;
  OMcp11_25 = qd(5)*C4;
  OMcp11_16 = OMcp11_15+qd(6)*ROcp11_15;
  OMcp11_26 = OMcp11_25+qd(6)*ROcp11_25;
  OMcp11_36 = qd(4)-qd(6)*S5;
  OPcp11_16 = -(qd(4)*qd(5)*C4+qd(6)*(qd(4)*ROcp11_25+OMcp11_25*S5)+qdd(5)*S4-qdd(6)*ROcp11_15);
  OPcp11_26 = -(qd(4)*qd(5)*S4-qd(6)*(qd(4)*ROcp11_15+OMcp11_15*S5)-qdd(5)*C4-qdd(6)*ROcp11_25);
  OPcp11_36 = qdd(4)-qd(5)*qd(6)*C5-qdd(6)*S5;

% = = Block_0_0_1_1_0_6 = = 
 
% Sensor Kinematics 


  ROcp11_116 = ROcp11_15*C16-ROcp11_76*S16;
  ROcp11_216 = ROcp11_25*C16-ROcp11_86*S16;
  ROcp11_316 = -(ROcp11_96*S16+C16*S5);
  ROcp11_716 = ROcp11_15*S16+ROcp11_76*C16;
  ROcp11_816 = ROcp11_25*S16+ROcp11_86*C16;
  ROcp11_916 = ROcp11_96*C16-S16*S5;
  RLcp11_116 = ROcp11_46*s.dpt(2,3);
  RLcp11_216 = ROcp11_56*s.dpt(2,3);
  RLcp11_316 = ROcp11_66*s.dpt(2,3);
  ORcp11_116 = OMcp11_26*RLcp11_316-OMcp11_36*RLcp11_216;
  ORcp11_216 = -(OMcp11_16*RLcp11_316-OMcp11_36*RLcp11_116);
  ORcp11_316 = OMcp11_16*RLcp11_216-OMcp11_26*RLcp11_116;
  PxF1(1) = q(1)+RLcp11_116;
  PxF1(2) = q(2)+RLcp11_216;
  PxF1(3) = q(3)+RLcp11_316;
  RxF1(1,1) = ROcp11_116;
  RxF1(1,2) = ROcp11_216;
  RxF1(1,3) = ROcp11_316;
  RxF1(2,1) = ROcp11_46;
  RxF1(2,2) = ROcp11_56;
  RxF1(2,3) = ROcp11_66;
  RxF1(3,1) = ROcp11_716;
  RxF1(3,2) = ROcp11_816;
  RxF1(3,3) = ROcp11_916;
  VxF1(1) = qd(1)+ORcp11_116;
  VxF1(2) = qd(2)+ORcp11_216;
  VxF1(3) = qd(3)+ORcp11_316;
  OMxF1(1) = OMcp11_16+qd(16)*ROcp11_46;
  OMxF1(2) = OMcp11_26+qd(16)*ROcp11_56;
  OMxF1(3) = OMcp11_36+qd(16)*ROcp11_66;
  AxF1(1) = qdd(1)+OMcp11_26*ORcp11_316-OMcp11_36*ORcp11_216+OPcp11_26*RLcp11_316-OPcp11_36*RLcp11_216;
  AxF1(2) = qdd(2)-OMcp11_16*ORcp11_316+OMcp11_36*ORcp11_116-OPcp11_16*RLcp11_316+OPcp11_36*RLcp11_116;
  AxF1(3) = qdd(3)+OMcp11_16*ORcp11_216-OMcp11_26*ORcp11_116+OPcp11_16*RLcp11_216-OPcp11_26*RLcp11_116;
  OMPxF1(1) = OPcp11_16+qd(16)*(OMcp11_26*ROcp11_66-OMcp11_36*ROcp11_56)+qdd(16)*ROcp11_46;
  OMPxF1(2) = OPcp11_26-qd(16)*(OMcp11_16*ROcp11_66-OMcp11_36*ROcp11_46)+qdd(16)*ROcp11_56;
  OMPxF1(3) = OPcp11_36+qd(16)*(OMcp11_16*ROcp11_56-OMcp11_26*ROcp11_46)+qdd(16)*ROcp11_66;
 
% Sensor Forces Computation 

  SWr1 = usrfun.fext(PxF1,RxF1,VxF1,OMxF1,AxF1,OMPxF1,s,tsim,1);
 
% Sensor Dynamics : Forces projection on body-fixed frames 

  xfrc112 = ROcp11_116*SWr1(1)+ROcp11_216*SWr1(2)+ROcp11_316*SWr1(3);
  xfrc212 = ROcp11_46*SWr1(1)+ROcp11_56*SWr1(2)+ROcp11_66*SWr1(3);
  xfrc312 = ROcp11_716*SWr1(1)+ROcp11_816*SWr1(2)+ROcp11_916*SWr1(3);
  s.frc(1,16) = s.frc(1,16)+xfrc112;
  s.frc(2,16) = s.frc(2,16)+xfrc212;
  s.frc(3,16) = s.frc(3,16)+xfrc312;
  xtrq112 = ROcp11_116*SWr1(4)+ROcp11_216*SWr1(5)+ROcp11_316*SWr1(6);
  xtrq212 = ROcp11_46*SWr1(4)+ROcp11_56*SWr1(5)+ROcp11_66*SWr1(6);
  xtrq312 = ROcp11_716*SWr1(4)+ROcp11_816*SWr1(5)+ROcp11_916*SWr1(6);
  s.trq(1,16) = s.trq(1,16)+xtrq112-xfrc212*SWr1(9)+xfrc312*SWr1(8);
  s.trq(2,16) = s.trq(2,16)+xtrq212+xfrc112*SWr1(9)-xfrc312*SWr1(7);
  s.trq(3,16) = s.trq(3,16)+xtrq312-xfrc112*SWr1(8)+xfrc212*SWr1(7);

% = = Block_0_0_1_2_0_1 = = 
 
% Sensor Kinematics 


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
  OMcp12_15 = -qd(5)*S4;
  OMcp12_25 = qd(5)*C4;
  OMcp12_16 = OMcp12_15+qd(6)*ROcp12_15;
  OMcp12_26 = OMcp12_25+qd(6)*ROcp12_25;
  OMcp12_36 = qd(4)-qd(6)*S5;
  OPcp12_16 = -(qd(4)*qd(5)*C4+qd(6)*(qd(4)*ROcp12_25+OMcp12_25*S5)+qdd(5)*S4-qdd(6)*ROcp12_15);
  OPcp12_26 = -(qd(4)*qd(5)*S4-qd(6)*(qd(4)*ROcp12_15+OMcp12_15*S5)-qdd(5)*C4-qdd(6)*ROcp12_25);
  OPcp12_36 = qdd(4)-qd(5)*qd(6)*C5-qdd(6)*S5;

% = = Block_0_0_1_2_0_6 = = 
 
% Sensor Kinematics 


  ROcp12_116 = ROcp12_15*C16-ROcp12_76*S16;
  ROcp12_216 = ROcp12_25*C16-ROcp12_86*S16;
  ROcp12_316 = -(ROcp12_96*S16+C16*S5);
  ROcp12_716 = ROcp12_15*S16+ROcp12_76*C16;
  ROcp12_816 = ROcp12_25*S16+ROcp12_86*C16;
  ROcp12_916 = ROcp12_96*C16-S16*S5;
  RLcp12_116 = ROcp12_46*s.dpt(2,3);
  RLcp12_216 = ROcp12_56*s.dpt(2,3);
  RLcp12_316 = ROcp12_66*s.dpt(2,3);
  OMcp12_116 = OMcp12_16+qd(16)*ROcp12_46;
  OMcp12_216 = OMcp12_26+qd(16)*ROcp12_56;
  OMcp12_316 = OMcp12_36+qd(16)*ROcp12_66;
  ORcp12_116 = OMcp12_26*RLcp12_316-OMcp12_36*RLcp12_216;
  ORcp12_216 = -(OMcp12_16*RLcp12_316-OMcp12_36*RLcp12_116);
  ORcp12_316 = OMcp12_16*RLcp12_216-OMcp12_26*RLcp12_116;
  OPcp12_116 = OPcp12_16+qd(16)*(OMcp12_26*ROcp12_66-OMcp12_36*ROcp12_56)+qdd(16)*ROcp12_46;
  OPcp12_216 = OPcp12_26-qd(16)*(OMcp12_16*ROcp12_66-OMcp12_36*ROcp12_46)+qdd(16)*ROcp12_56;
  OPcp12_316 = OPcp12_36+qd(16)*(OMcp12_16*ROcp12_56-OMcp12_26*ROcp12_46)+qdd(16)*ROcp12_66;
  RLcp12_130 = ROcp12_716*s.dpt(3,19);
  RLcp12_230 = ROcp12_816*s.dpt(3,19);
  RLcp12_330 = ROcp12_916*s.dpt(3,19);
  ORcp12_130 = OMcp12_216*RLcp12_330-OMcp12_316*RLcp12_230;
  ORcp12_230 = -(OMcp12_116*RLcp12_330-OMcp12_316*RLcp12_130);
  ORcp12_330 = OMcp12_116*RLcp12_230-OMcp12_216*RLcp12_130;
  PxF2(1) = q(1)+RLcp12_116+RLcp12_130;
  PxF2(2) = q(2)+RLcp12_216+RLcp12_230;
  PxF2(3) = q(3)+RLcp12_316+RLcp12_330;
  RxF2(1,1) = ROcp12_116;
  RxF2(1,2) = ROcp12_216;
  RxF2(1,3) = ROcp12_316;
  RxF2(2,1) = ROcp12_46;
  RxF2(2,2) = ROcp12_56;
  RxF2(2,3) = ROcp12_66;
  RxF2(3,1) = ROcp12_716;
  RxF2(3,2) = ROcp12_816;
  RxF2(3,3) = ROcp12_916;
  VxF2(1) = qd(1)+ORcp12_116+ORcp12_130;
  VxF2(2) = qd(2)+ORcp12_216+ORcp12_230;
  VxF2(3) = qd(3)+ORcp12_316+ORcp12_330;
  OMxF2(1) = OMcp12_116;
  OMxF2(2) = OMcp12_216;
  OMxF2(3) = OMcp12_316;
  AxF2(1) = qdd(1)+OMcp12_216*ORcp12_330+OMcp12_26*ORcp12_316-OMcp12_316*ORcp12_230-OMcp12_36*ORcp12_216+OPcp12_216*RLcp12_330+OPcp12_26*...
 RLcp12_316-OPcp12_316*RLcp12_230-OPcp12_36*RLcp12_216;
  AxF2(2) = qdd(2)-OMcp12_116*ORcp12_330-OMcp12_16*ORcp12_316+OMcp12_316*ORcp12_130+OMcp12_36*ORcp12_116-OPcp12_116*RLcp12_330-OPcp12_16*...
 RLcp12_316+OPcp12_316*RLcp12_130+OPcp12_36*RLcp12_116;
  AxF2(3) = qdd(3)+OMcp12_116*ORcp12_230+OMcp12_16*ORcp12_216-OMcp12_216*ORcp12_130-OMcp12_26*ORcp12_116+OPcp12_116*RLcp12_230+OPcp12_16*...
 RLcp12_216-OPcp12_216*RLcp12_130-OPcp12_26*RLcp12_116;
  OMPxF2(1) = OPcp12_116;
  OMPxF2(2) = OPcp12_216;
  OMPxF2(3) = OPcp12_316;
 
% Sensor Forces Computation 

  SWr2 = usrfun.fext(PxF2,RxF2,VxF2,OMxF2,AxF2,OMPxF2,s,tsim,2);
 
% Sensor Dynamics : Forces projection on body-fixed frames 

  xfrc113 = ROcp12_116*SWr2(1)+ROcp12_216*SWr2(2)+ROcp12_316*SWr2(3);
  xfrc213 = ROcp12_46*SWr2(1)+ROcp12_56*SWr2(2)+ROcp12_66*SWr2(3);
  xfrc313 = ROcp12_716*SWr2(1)+ROcp12_816*SWr2(2)+ROcp12_916*SWr2(3);
  frc(1,16) = s.frc(1,16)+xfrc113;
  frc(2,16) = s.frc(2,16)+xfrc213;
  frc(3,16) = s.frc(3,16)+xfrc313;
  xtrq113 = ROcp12_116*SWr2(4)+ROcp12_216*SWr2(5)+ROcp12_316*SWr2(6);
  xtrq213 = ROcp12_46*SWr2(4)+ROcp12_56*SWr2(5)+ROcp12_66*SWr2(6);
  xtrq313 = ROcp12_716*SWr2(4)+ROcp12_816*SWr2(5)+ROcp12_916*SWr2(6);
  trq(1,16) = s.trq(1,16)+xtrq113-xfrc213*SWr2(9)+xfrc313*SWr2(8);
  trq(2,16) = s.trq(2,16)+xtrq213+xfrc113*SWr2(9)-xfrc313*SWr2(7);
  trq(3,16) = s.trq(3,16)+xtrq313-xfrc113*SWr2(8)+xfrc213*SWr2(7);

% = = Block_0_0_1_3_0_1 = = 
 
% Sensor Kinematics 


  ROcp13_15 = C4*C5;
  ROcp13_25 = S4*C5;
  ROcp13_75 = C4*S5;
  ROcp13_85 = S4*S5;
  ROcp13_46 = ROcp13_75*S6-S4*C6;
  ROcp13_56 = ROcp13_85*S6+C4*C6;
  ROcp13_66 = C5*S6;
  ROcp13_76 = ROcp13_75*C6+S4*S6;
  ROcp13_86 = ROcp13_85*C6-C4*S6;
  ROcp13_96 = C5*C6;
  OMcp13_15 = -qd(5)*S4;
  OMcp13_25 = qd(5)*C4;
  OMcp13_16 = OMcp13_15+qd(6)*ROcp13_15;
  OMcp13_26 = OMcp13_25+qd(6)*ROcp13_25;
  OMcp13_36 = qd(4)-qd(6)*S5;
  OPcp13_16 = -(qd(4)*qd(5)*C4+qd(6)*(qd(4)*ROcp13_25+OMcp13_25*S5)+qdd(5)*S4-qdd(6)*ROcp13_15);
  OPcp13_26 = -(qd(4)*qd(5)*S4-qd(6)*(qd(4)*ROcp13_15+OMcp13_15*S5)-qdd(5)*C4-qdd(6)*ROcp13_25);
  OPcp13_36 = qdd(4)-qd(5)*qd(6)*C5-qdd(6)*S5;

% = = Block_0_0_1_3_0_7 = = 
 
% Sensor Kinematics 


  ROcp13_117 = ROcp13_15*C17-ROcp13_76*S17;
  ROcp13_217 = ROcp13_25*C17-ROcp13_86*S17;
  ROcp13_317 = -(ROcp13_96*S17+C17*S5);
  ROcp13_717 = ROcp13_15*S17+ROcp13_76*C17;
  ROcp13_817 = ROcp13_25*S17+ROcp13_86*C17;
  ROcp13_917 = ROcp13_96*C17-S17*S5;
  RLcp13_117 = ROcp13_46*s.dpt(2,4);
  RLcp13_217 = ROcp13_56*s.dpt(2,4);
  RLcp13_317 = ROcp13_66*s.dpt(2,4);
  ORcp13_117 = OMcp13_26*RLcp13_317-OMcp13_36*RLcp13_217;
  ORcp13_217 = -(OMcp13_16*RLcp13_317-OMcp13_36*RLcp13_117);
  ORcp13_317 = OMcp13_16*RLcp13_217-OMcp13_26*RLcp13_117;
  PxF3(1) = q(1)+RLcp13_117;
  PxF3(2) = q(2)+RLcp13_217;
  PxF3(3) = q(3)+RLcp13_317;
  RxF3(1,1) = ROcp13_117;
  RxF3(1,2) = ROcp13_217;
  RxF3(1,3) = ROcp13_317;
  RxF3(2,1) = ROcp13_46;
  RxF3(2,2) = ROcp13_56;
  RxF3(2,3) = ROcp13_66;
  RxF3(3,1) = ROcp13_717;
  RxF3(3,2) = ROcp13_817;
  RxF3(3,3) = ROcp13_917;
  VxF3(1) = qd(1)+ORcp13_117;
  VxF3(2) = qd(2)+ORcp13_217;
  VxF3(3) = qd(3)+ORcp13_317;
  OMxF3(1) = OMcp13_16+qd(17)*ROcp13_46;
  OMxF3(2) = OMcp13_26+qd(17)*ROcp13_56;
  OMxF3(3) = OMcp13_36+qd(17)*ROcp13_66;
  AxF3(1) = qdd(1)+OMcp13_26*ORcp13_317-OMcp13_36*ORcp13_217+OPcp13_26*RLcp13_317-OPcp13_36*RLcp13_217;
  AxF3(2) = qdd(2)-OMcp13_16*ORcp13_317+OMcp13_36*ORcp13_117-OPcp13_16*RLcp13_317+OPcp13_36*RLcp13_117;
  AxF3(3) = qdd(3)+OMcp13_16*ORcp13_217-OMcp13_26*ORcp13_117+OPcp13_16*RLcp13_217-OPcp13_26*RLcp13_117;
  OMPxF3(1) = OPcp13_16+qd(17)*(OMcp13_26*ROcp13_66-OMcp13_36*ROcp13_56)+qdd(17)*ROcp13_46;
  OMPxF3(2) = OPcp13_26-qd(17)*(OMcp13_16*ROcp13_66-OMcp13_36*ROcp13_46)+qdd(17)*ROcp13_56;
  OMPxF3(3) = OPcp13_36+qd(17)*(OMcp13_16*ROcp13_56-OMcp13_26*ROcp13_46)+qdd(17)*ROcp13_66;
 
% Sensor Forces Computation 

  SWr3 = usrfun.fext(PxF3,RxF3,VxF3,OMxF3,AxF3,OMPxF3,s,tsim,3);
 
% Sensor Dynamics : Forces projection on body-fixed frames 

  xfrc114 = ROcp13_117*SWr3(1)+ROcp13_217*SWr3(2)+ROcp13_317*SWr3(3);
  xfrc214 = ROcp13_46*SWr3(1)+ROcp13_56*SWr3(2)+ROcp13_66*SWr3(3);
  xfrc314 = ROcp13_717*SWr3(1)+ROcp13_817*SWr3(2)+ROcp13_917*SWr3(3);
  s.frc(1,17) = s.frc(1,17)+xfrc114;
  s.frc(2,17) = s.frc(2,17)+xfrc214;
  s.frc(3,17) = s.frc(3,17)+xfrc314;
  xtrq114 = ROcp13_117*SWr3(4)+ROcp13_217*SWr3(5)+ROcp13_317*SWr3(6);
  xtrq214 = ROcp13_46*SWr3(4)+ROcp13_56*SWr3(5)+ROcp13_66*SWr3(6);
  xtrq314 = ROcp13_717*SWr3(4)+ROcp13_817*SWr3(5)+ROcp13_917*SWr3(6);
  s.trq(1,17) = s.trq(1,17)+xtrq114-xfrc214*SWr3(9)+xfrc314*SWr3(8);
  s.trq(2,17) = s.trq(2,17)+xtrq214+xfrc114*SWr3(9)-xfrc314*SWr3(7);
  s.trq(3,17) = s.trq(3,17)+xtrq314-xfrc114*SWr3(8)+xfrc214*SWr3(7);

% = = Block_0_0_1_4_0_1 = = 
 
% Sensor Kinematics 


  ROcp14_15 = C4*C5;
  ROcp14_25 = S4*C5;
  ROcp14_75 = C4*S5;
  ROcp14_85 = S4*S5;
  ROcp14_46 = ROcp14_75*S6-S4*C6;
  ROcp14_56 = ROcp14_85*S6+C4*C6;
  ROcp14_66 = C5*S6;
  ROcp14_76 = ROcp14_75*C6+S4*S6;
  ROcp14_86 = ROcp14_85*C6-C4*S6;
  ROcp14_96 = C5*C6;
  OMcp14_15 = -qd(5)*S4;
  OMcp14_25 = qd(5)*C4;
  OMcp14_16 = OMcp14_15+qd(6)*ROcp14_15;
  OMcp14_26 = OMcp14_25+qd(6)*ROcp14_25;
  OMcp14_36 = qd(4)-qd(6)*S5;
  OPcp14_16 = -(qd(4)*qd(5)*C4+qd(6)*(qd(4)*ROcp14_25+OMcp14_25*S5)+qdd(5)*S4-qdd(6)*ROcp14_15);
  OPcp14_26 = -(qd(4)*qd(5)*S4-qd(6)*(qd(4)*ROcp14_15+OMcp14_15*S5)-qdd(5)*C4-qdd(6)*ROcp14_25);
  OPcp14_36 = qdd(4)-qd(5)*qd(6)*C5-qdd(6)*S5;

% = = Block_0_0_1_4_0_7 = = 
 
% Sensor Kinematics 


  ROcp14_117 = ROcp14_15*C17-ROcp14_76*S17;
  ROcp14_217 = ROcp14_25*C17-ROcp14_86*S17;
  ROcp14_317 = -(ROcp14_96*S17+C17*S5);
  ROcp14_717 = ROcp14_15*S17+ROcp14_76*C17;
  ROcp14_817 = ROcp14_25*S17+ROcp14_86*C17;
  ROcp14_917 = ROcp14_96*C17-S17*S5;
  RLcp14_117 = ROcp14_46*s.dpt(2,4);
  RLcp14_217 = ROcp14_56*s.dpt(2,4);
  RLcp14_317 = ROcp14_66*s.dpt(2,4);
  OMcp14_117 = OMcp14_16+qd(17)*ROcp14_46;
  OMcp14_217 = OMcp14_26+qd(17)*ROcp14_56;
  OMcp14_317 = OMcp14_36+qd(17)*ROcp14_66;
  ORcp14_117 = OMcp14_26*RLcp14_317-OMcp14_36*RLcp14_217;
  ORcp14_217 = -(OMcp14_16*RLcp14_317-OMcp14_36*RLcp14_117);
  ORcp14_317 = OMcp14_16*RLcp14_217-OMcp14_26*RLcp14_117;
  OPcp14_117 = OPcp14_16+qd(17)*(OMcp14_26*ROcp14_66-OMcp14_36*ROcp14_56)+qdd(17)*ROcp14_46;
  OPcp14_217 = OPcp14_26-qd(17)*(OMcp14_16*ROcp14_66-OMcp14_36*ROcp14_46)+qdd(17)*ROcp14_56;
  OPcp14_317 = OPcp14_36+qd(17)*(OMcp14_16*ROcp14_56-OMcp14_26*ROcp14_46)+qdd(17)*ROcp14_66;
  RLcp14_132 = ROcp14_717*s.dpt(3,21);
  RLcp14_232 = ROcp14_817*s.dpt(3,21);
  RLcp14_332 = ROcp14_917*s.dpt(3,21);
  ORcp14_132 = OMcp14_217*RLcp14_332-OMcp14_317*RLcp14_232;
  ORcp14_232 = -(OMcp14_117*RLcp14_332-OMcp14_317*RLcp14_132);
  ORcp14_332 = OMcp14_117*RLcp14_232-OMcp14_217*RLcp14_132;
  PxF4(1) = q(1)+RLcp14_117+RLcp14_132;
  PxF4(2) = q(2)+RLcp14_217+RLcp14_232;
  PxF4(3) = q(3)+RLcp14_317+RLcp14_332;
  RxF4(1,1) = ROcp14_117;
  RxF4(1,2) = ROcp14_217;
  RxF4(1,3) = ROcp14_317;
  RxF4(2,1) = ROcp14_46;
  RxF4(2,2) = ROcp14_56;
  RxF4(2,3) = ROcp14_66;
  RxF4(3,1) = ROcp14_717;
  RxF4(3,2) = ROcp14_817;
  RxF4(3,3) = ROcp14_917;
  VxF4(1) = qd(1)+ORcp14_117+ORcp14_132;
  VxF4(2) = qd(2)+ORcp14_217+ORcp14_232;
  VxF4(3) = qd(3)+ORcp14_317+ORcp14_332;
  OMxF4(1) = OMcp14_117;
  OMxF4(2) = OMcp14_217;
  OMxF4(3) = OMcp14_317;
  AxF4(1) = qdd(1)+OMcp14_217*ORcp14_332+OMcp14_26*ORcp14_317-OMcp14_317*ORcp14_232-OMcp14_36*ORcp14_217+OPcp14_217*RLcp14_332+OPcp14_26*...
 RLcp14_317-OPcp14_317*RLcp14_232-OPcp14_36*RLcp14_217;
  AxF4(2) = qdd(2)-OMcp14_117*ORcp14_332-OMcp14_16*ORcp14_317+OMcp14_317*ORcp14_132+OMcp14_36*ORcp14_117-OPcp14_117*RLcp14_332-OPcp14_16*...
 RLcp14_317+OPcp14_317*RLcp14_132+OPcp14_36*RLcp14_117;
  AxF4(3) = qdd(3)+OMcp14_117*ORcp14_232+OMcp14_16*ORcp14_217-OMcp14_217*ORcp14_132-OMcp14_26*ORcp14_117+OPcp14_117*RLcp14_232+OPcp14_16*...
 RLcp14_217-OPcp14_217*RLcp14_132-OPcp14_26*RLcp14_117;
  OMPxF4(1) = OPcp14_117;
  OMPxF4(2) = OPcp14_217;
  OMPxF4(3) = OPcp14_317;
 
% Sensor Forces Computation 

  SWr4 = usrfun.fext(PxF4,RxF4,VxF4,OMxF4,AxF4,OMPxF4,s,tsim,4);
 
% Sensor Dynamics : Forces projection on body-fixed frames 

  xfrc115 = ROcp14_117*SWr4(1)+ROcp14_217*SWr4(2)+ROcp14_317*SWr4(3);
  xfrc215 = ROcp14_46*SWr4(1)+ROcp14_56*SWr4(2)+ROcp14_66*SWr4(3);
  xfrc315 = ROcp14_717*SWr4(1)+ROcp14_817*SWr4(2)+ROcp14_917*SWr4(3);
  frc(1,17) = s.frc(1,17)+xfrc115;
  frc(2,17) = s.frc(2,17)+xfrc215;
  frc(3,17) = s.frc(3,17)+xfrc315;
  xtrq115 = ROcp14_117*SWr4(4)+ROcp14_217*SWr4(5)+ROcp14_317*SWr4(6);
  xtrq215 = ROcp14_46*SWr4(4)+ROcp14_56*SWr4(5)+ROcp14_66*SWr4(6);
  xtrq315 = ROcp14_717*SWr4(4)+ROcp14_817*SWr4(5)+ROcp14_917*SWr4(6);
  trq(1,17) = s.trq(1,17)+xtrq115-xfrc215*SWr4(9)+xfrc315*SWr4(8);
  trq(2,17) = s.trq(2,17)+xtrq215+xfrc115*SWr4(9)-xfrc315*SWr4(7);
  trq(3,17) = s.trq(3,17)+xtrq315-xfrc115*SWr4(8)+xfrc215*SWr4(7);

% = = Block_0_0_1_4_1_0 = = 
 
% Symbolic Outputs  

  frc(1,6) = s.frc(1,6);
  frc(2,6) = s.frc(2,6);
  frc(3,6) = s.frc(3,6);
  frc(1,7) = s.frc(1,7);
  frc(2,7) = s.frc(2,7);
  frc(3,7) = s.frc(3,7);
  frc(1,8) = s.frc(1,8);
  frc(2,8) = s.frc(2,8);
  frc(3,8) = s.frc(3,8);
  frc(1,9) = s.frc(1,9);
  frc(2,9) = s.frc(2,9);
  frc(3,9) = s.frc(3,9);
  frc(1,10) = s.frc(1,10);
  frc(2,10) = s.frc(2,10);
  frc(3,10) = s.frc(3,10);
  frc(1,11) = s.frc(1,11);
  frc(2,11) = s.frc(2,11);
  frc(3,11) = s.frc(3,11);
  frc(1,12) = s.frc(1,12);
  frc(2,12) = s.frc(2,12);
  frc(3,12) = s.frc(3,12);
  frc(1,13) = s.frc(1,13);
  frc(2,13) = s.frc(2,13);
  frc(3,13) = s.frc(3,13);
  frc(1,14) = s.frc(1,14);
  frc(2,14) = s.frc(2,14);
  frc(3,14) = s.frc(3,14);
  frc(1,15) = s.frc(1,15);
  frc(2,15) = s.frc(2,15);
  frc(3,15) = s.frc(3,15);
  trq(1,6) = s.trq(1,6);
  trq(2,6) = s.trq(2,6);
  trq(3,6) = s.trq(3,6);
  trq(1,7) = s.trq(1,7);
  trq(2,7) = s.trq(2,7);
  trq(3,7) = s.trq(3,7);
  trq(1,8) = s.trq(1,8);
  trq(2,8) = s.trq(2,8);
  trq(3,8) = s.trq(3,8);
  trq(1,9) = s.trq(1,9);
  trq(2,9) = s.trq(2,9);
  trq(3,9) = s.trq(3,9);
  trq(1,10) = s.trq(1,10);
  trq(2,10) = s.trq(2,10);
  trq(3,10) = s.trq(3,10);
  trq(1,11) = s.trq(1,11);
  trq(2,11) = s.trq(2,11);
  trq(3,11) = s.trq(3,11);
  trq(1,12) = s.trq(1,12);
  trq(2,12) = s.trq(2,12);
  trq(3,12) = s.trq(3,12);
  trq(1,13) = s.trq(1,13);
  trq(2,13) = s.trq(2,13);
  trq(3,13) = s.trq(3,13);
  trq(1,14) = s.trq(1,14);
  trq(2,14) = s.trq(2,14);
  trq(3,14) = s.trq(3,14);
  trq(1,15) = s.trq(1,15);
  trq(2,15) = s.trq(2,15);
  trq(3,15) = s.trq(3,15);

% ====== END Task 0 ====== 

  

