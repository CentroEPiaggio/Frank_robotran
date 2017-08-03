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
%	==> Function : F 6 : Sensors Kinematical Informations (sens) 
%	==> Flops complexity : 2659
%
%	==> Generation Time :  0.040 seconds
%	==> Post-Processing :  0.040 seconds
%
%-------------------------------------------------------------
%
function [sens] = gensensor(s,tsim,usrfun,isens)

 sens.P = zeros(3,1);
 sens.R = zeros(3,3);
 sens.V = zeros(3,1);
 sens.OM = zeros(3,1);
 sens.A = zeros(3,1);
 sens.OMP = zeros(3,1);
 sens.J = zeros(6,17);

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

% = = Block_0_0_0_0_0_2 = = 
 
% Augmented Joint Position Vectors   

  Dz73 = q(7)+s.dpt(3,2);

% = = Block_0_0_0_0_0_3 = = 
 
% Trigonometric Variables  

  C8 = cos(q(8));
  S8 = sin(q(8));
  C9 = cos(q(9));
  S9 = sin(q(9));
  C10 = cos(q(10));
  S10 = sin(q(10));

% = = Block_0_0_0_0_0_4 = = 
 
% Trigonometric Variables  

  C11 = cos(q(11));
  S11 = sin(q(11));
  C12 = cos(q(12));
  S12 = sin(q(12));
  C13 = cos(q(13));
  S13 = sin(q(13));

% = = Block_0_0_0_0_0_5 = = 
 
% Trigonometric Variables  

  C14 = cos(q(14));
  S14 = sin(q(14));
  C15 = cos(q(15));
  S15 = sin(q(15));

% = = Block_0_0_0_0_0_6 = = 
 
% Trigonometric Variables  

  C16 = cos(q(16));
  S16 = sin(q(16));

% = = Block_0_0_0_0_0_7 = = 
 
% Trigonometric Variables  

  C17 = cos(q(17));
  S17 = sin(q(17));

% ====== END Task 0 ====== 

% ===== BEGIN task 1 ===== 
 
switch isens

 
% 
case 1, 


% = = Block_1_0_0_1_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = q(1);
    sens.R(1,1) = (1.0);
    sens.R(2,2) = (1.0);
    sens.R(3,3) = (1.0);
    sens.V(1) = qd(1);
    sens.A(1) = qdd(1);
 
% 
case 2, 


% = = Block_1_0_0_2_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = q(1);
    sens.P(2) = q(2);
    sens.R(1,1) = (1.0);
    sens.R(2,2) = (1.0);
    sens.R(3,3) = (1.0);
    sens.V(1) = qd(1);
    sens.V(2) = qd(2);
    sens.A(1) = qdd(1);
    sens.A(2) = qdd(2);
 
% 
case 3, 


% = = Block_1_0_0_3_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = q(1);
    sens.P(2) = q(2);
    sens.P(3) = q(3);
    sens.R(1,1) = (1.0);
    sens.R(2,2) = (1.0);
    sens.R(3,3) = (1.0);
    sens.V(1) = qd(1);
    sens.V(2) = qd(2);
    sens.V(3) = qd(3);
    sens.A(1) = qdd(1);
    sens.A(2) = qdd(2);
    sens.A(3) = qdd(3);
 
% 
case 4, 


% = = Block_1_0_0_4_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = q(1);
    sens.P(2) = q(2);
    sens.P(3) = q(3);
    sens.R(1,1) = C4;
    sens.R(1,2) = S4;
    sens.R(2,1) = -S4;
    sens.R(2,2) = C4;
    sens.R(3,3) = (1.0);
    sens.V(1) = qd(1);
    sens.V(2) = qd(2);
    sens.V(3) = qd(3);
    sens.OM(3) = qd(4);
    sens.A(1) = qdd(1);
    sens.A(2) = qdd(2);
    sens.A(3) = qdd(3);
    sens.OMP(3) = qdd(4);
 
% 
case 5, 


% = = Block_1_0_0_5_0_1 = = 
 
% Sensor Kinematics 


    ROcp4_15 = C4*C5;
    ROcp4_25 = S4*C5;
    ROcp4_75 = C4*S5;
    ROcp4_85 = S4*S5;
    OMcp4_15 = -qd(5)*S4;
    OMcp4_25 = qd(5)*C4;
    OPcp4_15 = -(qdd(5)*S4+qd(4)*qd(5)*C4);
    OPcp4_25 = qdd(5)*C4-qd(4)*qd(5)*S4;

% = = Block_1_0_0_5_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = q(1);
    sens.P(2) = q(2);
    sens.P(3) = q(3);
    sens.R(1,1) = ROcp4_15;
    sens.R(1,2) = ROcp4_25;
    sens.R(1,3) = -S5;
    sens.R(2,1) = -S4;
    sens.R(2,2) = C4;
    sens.R(3,1) = ROcp4_75;
    sens.R(3,2) = ROcp4_85;
    sens.R(3,3) = C5;
    sens.V(1) = qd(1);
    sens.V(2) = qd(2);
    sens.V(3) = qd(3);
    sens.OM(1) = OMcp4_15;
    sens.OM(2) = OMcp4_25;
    sens.OM(3) = qd(4);
    sens.A(1) = qdd(1);
    sens.A(2) = qdd(2);
    sens.A(3) = qdd(3);
    sens.OMP(1) = OPcp4_15;
    sens.OMP(2) = OPcp4_25;
    sens.OMP(3) = qdd(4);
 
% 
case 6, 


% = = Block_1_0_0_6_0_1 = = 
 
% Sensor Kinematics 


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
    OMcp5_15 = -qd(5)*S4;
    OMcp5_25 = qd(5)*C4;
    OMcp5_16 = OMcp5_15+ROcp5_15*qd(6);
    OMcp5_26 = OMcp5_25+ROcp5_25*qd(6);
    OMcp5_36 = qd(4)-qd(6)*S5;
    OPcp5_16 = ROcp5_15*qdd(6)-qdd(5)*S4-qd(4)*qd(5)*C4-qd(6)*(OMcp5_25*S5+ROcp5_25*qd(4));
    OPcp5_26 = ROcp5_25*qdd(6)+qdd(5)*C4-qd(4)*qd(5)*S4+qd(6)*(OMcp5_15*S5+ROcp5_15*qd(4));
    OPcp5_36 = qdd(4)-qdd(6)*S5-qd(5)*qd(6)*C5;

% = = Block_1_0_0_6_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = q(1);
    sens.P(2) = q(2);
    sens.P(3) = q(3);
    sens.R(1,1) = ROcp5_15;
    sens.R(1,2) = ROcp5_25;
    sens.R(1,3) = -S5;
    sens.R(2,1) = ROcp5_46;
    sens.R(2,2) = ROcp5_56;
    sens.R(2,3) = ROcp5_66;
    sens.R(3,1) = ROcp5_76;
    sens.R(3,2) = ROcp5_86;
    sens.R(3,3) = ROcp5_96;
    sens.V(1) = qd(1);
    sens.V(2) = qd(2);
    sens.V(3) = qd(3);
    sens.OM(1) = OMcp5_16;
    sens.OM(2) = OMcp5_26;
    sens.OM(3) = OMcp5_36;
    sens.A(1) = qdd(1);
    sens.A(2) = qdd(2);
    sens.A(3) = qdd(3);
    sens.OMP(1) = OPcp5_16;
    sens.OMP(2) = OPcp5_26;
    sens.OMP(3) = OPcp5_36;
 
% 
case 7, 


% = = Block_1_0_0_7_0_1 = = 
 
% Sensor Kinematics 


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
    OMcp6_15 = -qd(5)*S4;
    OMcp6_25 = qd(5)*C4;
    OMcp6_16 = OMcp6_15+ROcp6_15*qd(6);
    OMcp6_26 = OMcp6_25+ROcp6_25*qd(6);
    OMcp6_36 = qd(4)-qd(6)*S5;
    OPcp6_16 = ROcp6_15*qdd(6)-qdd(5)*S4-qd(4)*qd(5)*C4-qd(6)*(OMcp6_25*S5+ROcp6_25*qd(4));
    OPcp6_26 = ROcp6_25*qdd(6)+qdd(5)*C4-qd(4)*qd(5)*S4+qd(6)*(OMcp6_15*S5+ROcp6_15*qd(4));
    OPcp6_36 = qdd(4)-qdd(6)*S5-qd(5)*qd(6)*C5;

% = = Block_1_0_0_7_0_2 = = 
 
% Sensor Kinematics 


    RLcp6_17 = Dz73*ROcp6_76;
    RLcp6_27 = Dz73*ROcp6_86;
    RLcp6_37 = Dz73*ROcp6_96;
    POcp6_17 = RLcp6_17+q(1);
    POcp6_27 = RLcp6_27+q(2);
    POcp6_37 = RLcp6_37+q(3);
    ORcp6_17 = OMcp6_26*RLcp6_37-OMcp6_36*RLcp6_27;
    ORcp6_27 = -(OMcp6_16*RLcp6_37-OMcp6_36*RLcp6_17);
    ORcp6_37 = OMcp6_16*RLcp6_27-OMcp6_26*RLcp6_17;
    VIcp6_17 = ORcp6_17+qd(1)+ROcp6_76*qd(7);
    VIcp6_27 = ORcp6_27+qd(2)+ROcp6_86*qd(7);
    VIcp6_37 = ORcp6_37+qd(3)+ROcp6_96*qd(7);
    ACcp6_17 = qdd(1)+OMcp6_26*ORcp6_37-OMcp6_36*ORcp6_27+OPcp6_26*RLcp6_37-OPcp6_36*RLcp6_27+ROcp6_76*qdd(7)+(2.0)*qd(7)*(OMcp6_26*ROcp6_96-OMcp6_36*...
 ROcp6_86);
    ACcp6_27 = qdd(2)-OMcp6_16*ORcp6_37+OMcp6_36*ORcp6_17-OPcp6_16*RLcp6_37+OPcp6_36*RLcp6_17+ROcp6_86*qdd(7)-(2.0)*qd(7)*(OMcp6_16*ROcp6_96-OMcp6_36*...
 ROcp6_76);
    ACcp6_37 = qdd(3)+OMcp6_16*ORcp6_27-OMcp6_26*ORcp6_17+OPcp6_16*RLcp6_27-OPcp6_26*RLcp6_17+ROcp6_96*qdd(7)+(2.0)*qd(7)*(OMcp6_16*ROcp6_86-OMcp6_26*...
 ROcp6_76);

% = = Block_1_0_0_7_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = POcp6_17;
    sens.P(2) = POcp6_27;
    sens.P(3) = POcp6_37;
    sens.R(1,1) = ROcp6_15;
    sens.R(1,2) = ROcp6_25;
    sens.R(1,3) = -S5;
    sens.R(2,1) = ROcp6_46;
    sens.R(2,2) = ROcp6_56;
    sens.R(2,3) = ROcp6_66;
    sens.R(3,1) = ROcp6_76;
    sens.R(3,2) = ROcp6_86;
    sens.R(3,3) = ROcp6_96;
    sens.V(1) = VIcp6_17;
    sens.V(2) = VIcp6_27;
    sens.V(3) = VIcp6_37;
    sens.OM(1) = OMcp6_16;
    sens.OM(2) = OMcp6_26;
    sens.OM(3) = OMcp6_36;
    sens.A(1) = ACcp6_17;
    sens.A(2) = ACcp6_27;
    sens.A(3) = ACcp6_37;
    sens.OMP(1) = OPcp6_16;
    sens.OMP(2) = OPcp6_26;
    sens.OMP(3) = OPcp6_36;
 
% 
case 8, 


% = = Block_1_0_0_8_0_1 = = 
 
% Sensor Kinematics 


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
    OMcp7_15 = -qd(5)*S4;
    OMcp7_25 = qd(5)*C4;
    OMcp7_16 = OMcp7_15+ROcp7_15*qd(6);
    OMcp7_26 = OMcp7_25+ROcp7_25*qd(6);
    OMcp7_36 = qd(4)-qd(6)*S5;
    OPcp7_16 = ROcp7_15*qdd(6)-qdd(5)*S4-qd(4)*qd(5)*C4-qd(6)*(OMcp7_25*S5+ROcp7_25*qd(4));
    OPcp7_26 = ROcp7_25*qdd(6)+qdd(5)*C4-qd(4)*qd(5)*S4+qd(6)*(OMcp7_15*S5+ROcp7_15*qd(4));
    OPcp7_36 = qdd(4)-qdd(6)*S5-qd(5)*qd(6)*C5;

% = = Block_1_0_0_8_0_2 = = 
 
% Sensor Kinematics 


    RLcp7_17 = Dz73*ROcp7_76;
    RLcp7_27 = Dz73*ROcp7_86;
    RLcp7_37 = Dz73*ROcp7_96;
    ORcp7_17 = OMcp7_26*RLcp7_37-OMcp7_36*RLcp7_27;
    ORcp7_27 = -(OMcp7_16*RLcp7_37-OMcp7_36*RLcp7_17);
    ORcp7_37 = OMcp7_16*RLcp7_27-OMcp7_26*RLcp7_17;

% = = Block_1_0_0_8_0_3 = = 
 
% Sensor Kinematics 


    ROcp7_18 = ROcp7_15*C8-ROcp7_76*S8;
    ROcp7_28 = ROcp7_25*C8-ROcp7_86*S8;
    ROcp7_38 = -(ROcp7_96*S8+S5*C8);
    ROcp7_78 = ROcp7_15*S8+ROcp7_76*C8;
    ROcp7_88 = ROcp7_25*S8+ROcp7_86*C8;
    ROcp7_98 = ROcp7_96*C8-S5*S8;
    RLcp7_18 = ROcp7_46*s.dpt(2,6);
    RLcp7_28 = ROcp7_56*s.dpt(2,6);
    RLcp7_38 = ROcp7_66*s.dpt(2,6);
    POcp7_18 = RLcp7_17+RLcp7_18+q(1);
    POcp7_28 = RLcp7_27+RLcp7_28+q(2);
    POcp7_38 = RLcp7_37+RLcp7_38+q(3);
    OMcp7_18 = OMcp7_16+ROcp7_46*qd(8);
    OMcp7_28 = OMcp7_26+ROcp7_56*qd(8);
    OMcp7_38 = OMcp7_36+ROcp7_66*qd(8);
    ORcp7_18 = OMcp7_26*RLcp7_38-OMcp7_36*RLcp7_28;
    ORcp7_28 = -(OMcp7_16*RLcp7_38-OMcp7_36*RLcp7_18);
    ORcp7_38 = OMcp7_16*RLcp7_28-OMcp7_26*RLcp7_18;
    VIcp7_18 = ORcp7_17+ORcp7_18+qd(1)+ROcp7_76*qd(7);
    VIcp7_28 = ORcp7_27+ORcp7_28+qd(2)+ROcp7_86*qd(7);
    VIcp7_38 = ORcp7_37+ORcp7_38+qd(3)+ROcp7_96*qd(7);
    OPcp7_18 = OPcp7_16+ROcp7_46*qdd(8)+qd(8)*(OMcp7_26*ROcp7_66-OMcp7_36*ROcp7_56);
    OPcp7_28 = OPcp7_26+ROcp7_56*qdd(8)-qd(8)*(OMcp7_16*ROcp7_66-OMcp7_36*ROcp7_46);
    OPcp7_38 = OPcp7_36+ROcp7_66*qdd(8)+qd(8)*(OMcp7_16*ROcp7_56-OMcp7_26*ROcp7_46);
    ACcp7_18 = qdd(1)+OMcp7_26*ORcp7_37+OMcp7_26*ORcp7_38-OMcp7_36*ORcp7_27-OMcp7_36*ORcp7_28+OPcp7_26*RLcp7_37+OPcp7_26*RLcp7_38-OPcp7_36*RLcp7_27...
 -OPcp7_36*RLcp7_28+ROcp7_76*qdd(7)+(2.0)*qd(7)*(OMcp7_26*ROcp7_96-OMcp7_36*ROcp7_86);
    ACcp7_28 = qdd(2)-OMcp7_16*ORcp7_37-OMcp7_16*ORcp7_38+OMcp7_36*ORcp7_17+OMcp7_36*ORcp7_18-OPcp7_16*RLcp7_37-OPcp7_16*RLcp7_38+OPcp7_36*RLcp7_17...
 +OPcp7_36*RLcp7_18+ROcp7_86*qdd(7)-(2.0)*qd(7)*(OMcp7_16*ROcp7_96-OMcp7_36*ROcp7_76);
    ACcp7_38 = qdd(3)+OMcp7_16*ORcp7_27+OMcp7_16*ORcp7_28-OMcp7_26*ORcp7_17-OMcp7_26*ORcp7_18+OPcp7_16*RLcp7_27+OPcp7_16*RLcp7_28-OPcp7_26*RLcp7_17...
 -OPcp7_26*RLcp7_18+ROcp7_96*qdd(7)+(2.0)*qd(7)*(OMcp7_16*ROcp7_86-OMcp7_26*ROcp7_76);

% = = Block_1_0_0_8_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = POcp7_18;
    sens.P(2) = POcp7_28;
    sens.P(3) = POcp7_38;
    sens.R(1,1) = ROcp7_18;
    sens.R(1,2) = ROcp7_28;
    sens.R(1,3) = ROcp7_38;
    sens.R(2,1) = ROcp7_46;
    sens.R(2,2) = ROcp7_56;
    sens.R(2,3) = ROcp7_66;
    sens.R(3,1) = ROcp7_78;
    sens.R(3,2) = ROcp7_88;
    sens.R(3,3) = ROcp7_98;
    sens.V(1) = VIcp7_18;
    sens.V(2) = VIcp7_28;
    sens.V(3) = VIcp7_38;
    sens.OM(1) = OMcp7_18;
    sens.OM(2) = OMcp7_28;
    sens.OM(3) = OMcp7_38;
    sens.A(1) = ACcp7_18;
    sens.A(2) = ACcp7_28;
    sens.A(3) = ACcp7_38;
    sens.OMP(1) = OPcp7_18;
    sens.OMP(2) = OPcp7_28;
    sens.OMP(3) = OPcp7_38;
 
% 
case 9, 


% = = Block_1_0_0_9_0_1 = = 
 
% Sensor Kinematics 


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
    OMcp8_15 = -qd(5)*S4;
    OMcp8_25 = qd(5)*C4;
    OMcp8_16 = OMcp8_15+ROcp8_15*qd(6);
    OMcp8_26 = OMcp8_25+ROcp8_25*qd(6);
    OMcp8_36 = qd(4)-qd(6)*S5;
    OPcp8_16 = ROcp8_15*qdd(6)-qdd(5)*S4-qd(4)*qd(5)*C4-qd(6)*(OMcp8_25*S5+ROcp8_25*qd(4));
    OPcp8_26 = ROcp8_25*qdd(6)+qdd(5)*C4-qd(4)*qd(5)*S4+qd(6)*(OMcp8_15*S5+ROcp8_15*qd(4));
    OPcp8_36 = qdd(4)-qdd(6)*S5-qd(5)*qd(6)*C5;

% = = Block_1_0_0_9_0_2 = = 
 
% Sensor Kinematics 


    RLcp8_17 = Dz73*ROcp8_76;
    RLcp8_27 = Dz73*ROcp8_86;
    RLcp8_37 = Dz73*ROcp8_96;
    ORcp8_17 = OMcp8_26*RLcp8_37-OMcp8_36*RLcp8_27;
    ORcp8_27 = -(OMcp8_16*RLcp8_37-OMcp8_36*RLcp8_17);
    ORcp8_37 = OMcp8_16*RLcp8_27-OMcp8_26*RLcp8_17;

% = = Block_1_0_0_9_0_3 = = 
 
% Sensor Kinematics 


    ROcp8_18 = ROcp8_15*C8-ROcp8_76*S8;
    ROcp8_28 = ROcp8_25*C8-ROcp8_86*S8;
    ROcp8_38 = -(ROcp8_96*S8+S5*C8);
    ROcp8_78 = ROcp8_15*S8+ROcp8_76*C8;
    ROcp8_88 = ROcp8_25*S8+ROcp8_86*C8;
    ROcp8_98 = ROcp8_96*C8-S5*S8;
    ROcp8_49 = ROcp8_46*C9+ROcp8_78*S9;
    ROcp8_59 = ROcp8_56*C9+ROcp8_88*S9;
    ROcp8_69 = ROcp8_66*C9+ROcp8_98*S9;
    ROcp8_79 = -(ROcp8_46*S9-ROcp8_78*C9);
    ROcp8_89 = -(ROcp8_56*S9-ROcp8_88*C9);
    ROcp8_99 = -(ROcp8_66*S9-ROcp8_98*C9);
    RLcp8_18 = ROcp8_46*s.dpt(2,6);
    RLcp8_28 = ROcp8_56*s.dpt(2,6);
    RLcp8_38 = ROcp8_66*s.dpt(2,6);
    POcp8_18 = RLcp8_17+RLcp8_18+q(1);
    POcp8_28 = RLcp8_27+RLcp8_28+q(2);
    POcp8_38 = RLcp8_37+RLcp8_38+q(3);
    OMcp8_18 = OMcp8_16+ROcp8_46*qd(8);
    OMcp8_28 = OMcp8_26+ROcp8_56*qd(8);
    OMcp8_38 = OMcp8_36+ROcp8_66*qd(8);
    ORcp8_18 = OMcp8_26*RLcp8_38-OMcp8_36*RLcp8_28;
    ORcp8_28 = -(OMcp8_16*RLcp8_38-OMcp8_36*RLcp8_18);
    ORcp8_38 = OMcp8_16*RLcp8_28-OMcp8_26*RLcp8_18;
    VIcp8_18 = ORcp8_17+ORcp8_18+qd(1)+ROcp8_76*qd(7);
    VIcp8_28 = ORcp8_27+ORcp8_28+qd(2)+ROcp8_86*qd(7);
    VIcp8_38 = ORcp8_37+ORcp8_38+qd(3)+ROcp8_96*qd(7);
    ACcp8_18 = qdd(1)+OMcp8_26*ORcp8_37+OMcp8_26*ORcp8_38-OMcp8_36*ORcp8_27-OMcp8_36*ORcp8_28+OPcp8_26*RLcp8_37+OPcp8_26*RLcp8_38-OPcp8_36*RLcp8_27...
 -OPcp8_36*RLcp8_28+ROcp8_76*qdd(7)+(2.0)*qd(7)*(OMcp8_26*ROcp8_96-OMcp8_36*ROcp8_86);
    ACcp8_28 = qdd(2)-OMcp8_16*ORcp8_37-OMcp8_16*ORcp8_38+OMcp8_36*ORcp8_17+OMcp8_36*ORcp8_18-OPcp8_16*RLcp8_37-OPcp8_16*RLcp8_38+OPcp8_36*RLcp8_17...
 +OPcp8_36*RLcp8_18+ROcp8_86*qdd(7)-(2.0)*qd(7)*(OMcp8_16*ROcp8_96-OMcp8_36*ROcp8_76);
    ACcp8_38 = qdd(3)+OMcp8_16*ORcp8_27+OMcp8_16*ORcp8_28-OMcp8_26*ORcp8_17-OMcp8_26*ORcp8_18+OPcp8_16*RLcp8_27+OPcp8_16*RLcp8_28-OPcp8_26*RLcp8_17...
 -OPcp8_26*RLcp8_18+ROcp8_96*qdd(7)+(2.0)*qd(7)*(OMcp8_16*ROcp8_86-OMcp8_26*ROcp8_76);
    OMcp8_19 = OMcp8_18+ROcp8_18*qd(9);
    OMcp8_29 = OMcp8_28+ROcp8_28*qd(9);
    OMcp8_39 = OMcp8_38+ROcp8_38*qd(9);
    OPcp8_19 = OPcp8_16+ROcp8_18*qdd(9)+ROcp8_46*qdd(8)+qd(8)*(OMcp8_26*ROcp8_66-OMcp8_36*ROcp8_56)+qd(9)*(OMcp8_28*ROcp8_38-OMcp8_38*ROcp8_28);
    OPcp8_29 = OPcp8_26+ROcp8_28*qdd(9)+ROcp8_56*qdd(8)-qd(8)*(OMcp8_16*ROcp8_66-OMcp8_36*ROcp8_46)-qd(9)*(OMcp8_18*ROcp8_38-OMcp8_38*ROcp8_18);
    OPcp8_39 = OPcp8_36+ROcp8_38*qdd(9)+ROcp8_66*qdd(8)+qd(8)*(OMcp8_16*ROcp8_56-OMcp8_26*ROcp8_46)+qd(9)*(OMcp8_18*ROcp8_28-OMcp8_28*ROcp8_18);

% = = Block_1_0_0_9_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = POcp8_18;
    sens.P(2) = POcp8_28;
    sens.P(3) = POcp8_38;
    sens.R(1,1) = ROcp8_18;
    sens.R(1,2) = ROcp8_28;
    sens.R(1,3) = ROcp8_38;
    sens.R(2,1) = ROcp8_49;
    sens.R(2,2) = ROcp8_59;
    sens.R(2,3) = ROcp8_69;
    sens.R(3,1) = ROcp8_79;
    sens.R(3,2) = ROcp8_89;
    sens.R(3,3) = ROcp8_99;
    sens.V(1) = VIcp8_18;
    sens.V(2) = VIcp8_28;
    sens.V(3) = VIcp8_38;
    sens.OM(1) = OMcp8_19;
    sens.OM(2) = OMcp8_29;
    sens.OM(3) = OMcp8_39;
    sens.A(1) = ACcp8_18;
    sens.A(2) = ACcp8_28;
    sens.A(3) = ACcp8_38;
    sens.OMP(1) = OPcp8_19;
    sens.OMP(2) = OPcp8_29;
    sens.OMP(3) = OPcp8_39;
 
% 
case 10, 


% = = Block_1_0_0_10_0_1 = = 
 
% Sensor Kinematics 


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
    OMcp9_15 = -qd(5)*S4;
    OMcp9_25 = qd(5)*C4;
    OMcp9_16 = OMcp9_15+ROcp9_15*qd(6);
    OMcp9_26 = OMcp9_25+ROcp9_25*qd(6);
    OMcp9_36 = qd(4)-qd(6)*S5;
    OPcp9_16 = ROcp9_15*qdd(6)-qdd(5)*S4-qd(4)*qd(5)*C4-qd(6)*(OMcp9_25*S5+ROcp9_25*qd(4));
    OPcp9_26 = ROcp9_25*qdd(6)+qdd(5)*C4-qd(4)*qd(5)*S4+qd(6)*(OMcp9_15*S5+ROcp9_15*qd(4));
    OPcp9_36 = qdd(4)-qdd(6)*S5-qd(5)*qd(6)*C5;

% = = Block_1_0_0_10_0_2 = = 
 
% Sensor Kinematics 


    RLcp9_17 = Dz73*ROcp9_76;
    RLcp9_27 = Dz73*ROcp9_86;
    RLcp9_37 = Dz73*ROcp9_96;
    ORcp9_17 = OMcp9_26*RLcp9_37-OMcp9_36*RLcp9_27;
    ORcp9_27 = -(OMcp9_16*RLcp9_37-OMcp9_36*RLcp9_17);
    ORcp9_37 = OMcp9_16*RLcp9_27-OMcp9_26*RLcp9_17;

% = = Block_1_0_0_10_0_3 = = 
 
% Sensor Kinematics 


    ROcp9_18 = ROcp9_15*C8-ROcp9_76*S8;
    ROcp9_28 = ROcp9_25*C8-ROcp9_86*S8;
    ROcp9_38 = -(ROcp9_96*S8+S5*C8);
    ROcp9_78 = ROcp9_15*S8+ROcp9_76*C8;
    ROcp9_88 = ROcp9_25*S8+ROcp9_86*C8;
    ROcp9_98 = ROcp9_96*C8-S5*S8;
    ROcp9_49 = ROcp9_46*C9+ROcp9_78*S9;
    ROcp9_59 = ROcp9_56*C9+ROcp9_88*S9;
    ROcp9_69 = ROcp9_66*C9+ROcp9_98*S9;
    ROcp9_79 = -(ROcp9_46*S9-ROcp9_78*C9);
    ROcp9_89 = -(ROcp9_56*S9-ROcp9_88*C9);
    ROcp9_99 = -(ROcp9_66*S9-ROcp9_98*C9);
    ROcp9_110 = ROcp9_18*C10+ROcp9_49*S10;
    ROcp9_210 = ROcp9_28*C10+ROcp9_59*S10;
    ROcp9_310 = ROcp9_38*C10+ROcp9_69*S10;
    ROcp9_410 = -(ROcp9_18*S10-ROcp9_49*C10);
    ROcp9_510 = -(ROcp9_28*S10-ROcp9_59*C10);
    ROcp9_610 = -(ROcp9_38*S10-ROcp9_69*C10);
    RLcp9_18 = ROcp9_46*s.dpt(2,6);
    RLcp9_28 = ROcp9_56*s.dpt(2,6);
    RLcp9_38 = ROcp9_66*s.dpt(2,6);
    OMcp9_18 = OMcp9_16+ROcp9_46*qd(8);
    OMcp9_28 = OMcp9_26+ROcp9_56*qd(8);
    OMcp9_38 = OMcp9_36+ROcp9_66*qd(8);
    ORcp9_18 = OMcp9_26*RLcp9_38-OMcp9_36*RLcp9_28;
    ORcp9_28 = -(OMcp9_16*RLcp9_38-OMcp9_36*RLcp9_18);
    ORcp9_38 = OMcp9_16*RLcp9_28-OMcp9_26*RLcp9_18;
    OMcp9_19 = OMcp9_18+ROcp9_18*qd(9);
    OMcp9_29 = OMcp9_28+ROcp9_28*qd(9);
    OMcp9_39 = OMcp9_38+ROcp9_38*qd(9);
    OPcp9_19 = OPcp9_16+ROcp9_18*qdd(9)+ROcp9_46*qdd(8)+qd(8)*(OMcp9_26*ROcp9_66-OMcp9_36*ROcp9_56)+qd(9)*(OMcp9_28*ROcp9_38-OMcp9_38*ROcp9_28);
    OPcp9_29 = OPcp9_26+ROcp9_28*qdd(9)+ROcp9_56*qdd(8)-qd(8)*(OMcp9_16*ROcp9_66-OMcp9_36*ROcp9_46)-qd(9)*(OMcp9_18*ROcp9_38-OMcp9_38*ROcp9_18);
    OPcp9_39 = OPcp9_36+ROcp9_38*qdd(9)+ROcp9_66*qdd(8)+qd(8)*(OMcp9_16*ROcp9_56-OMcp9_26*ROcp9_46)+qd(9)*(OMcp9_18*ROcp9_28-OMcp9_28*ROcp9_18);
    RLcp9_110 = ROcp9_79*s.dpt(3,11);
    RLcp9_210 = ROcp9_89*s.dpt(3,11);
    RLcp9_310 = ROcp9_99*s.dpt(3,11);
    POcp9_110 = RLcp9_110+RLcp9_17+RLcp9_18+q(1);
    POcp9_210 = RLcp9_210+RLcp9_27+RLcp9_28+q(2);
    POcp9_310 = RLcp9_310+RLcp9_37+RLcp9_38+q(3);
    OMcp9_110 = OMcp9_19+ROcp9_79*qd(10);
    OMcp9_210 = OMcp9_29+ROcp9_89*qd(10);
    OMcp9_310 = OMcp9_39+ROcp9_99*qd(10);
    ORcp9_110 = OMcp9_29*RLcp9_310-OMcp9_39*RLcp9_210;
    ORcp9_210 = -(OMcp9_19*RLcp9_310-OMcp9_39*RLcp9_110);
    ORcp9_310 = OMcp9_19*RLcp9_210-OMcp9_29*RLcp9_110;
    VIcp9_110 = ORcp9_110+ORcp9_17+ORcp9_18+qd(1)+ROcp9_76*qd(7);
    VIcp9_210 = ORcp9_210+ORcp9_27+ORcp9_28+qd(2)+ROcp9_86*qd(7);
    VIcp9_310 = ORcp9_310+ORcp9_37+ORcp9_38+qd(3)+ROcp9_96*qd(7);
    OPcp9_110 = OPcp9_19+ROcp9_79*qdd(10)+qd(10)*(OMcp9_29*ROcp9_99-OMcp9_39*ROcp9_89);
    OPcp9_210 = OPcp9_29+ROcp9_89*qdd(10)-qd(10)*(OMcp9_19*ROcp9_99-OMcp9_39*ROcp9_79);
    OPcp9_310 = OPcp9_39+ROcp9_99*qdd(10)+qd(10)*(OMcp9_19*ROcp9_89-OMcp9_29*ROcp9_79);
    ACcp9_110 = qdd(1)+OMcp9_26*ORcp9_37+OMcp9_26*ORcp9_38+OMcp9_29*ORcp9_310-OMcp9_36*ORcp9_27-OMcp9_36*ORcp9_28-OMcp9_39*ORcp9_210+OPcp9_26*...
 RLcp9_37+OPcp9_26*RLcp9_38+OPcp9_29*RLcp9_310-OPcp9_36*RLcp9_27-OPcp9_36*RLcp9_28-OPcp9_39*RLcp9_210+ROcp9_76*qdd(7)+(2.0)*qd(7)*(OMcp9_26*ROcp9_96-...
 OMcp9_36*ROcp9_86);
    ACcp9_210 = qdd(2)-OMcp9_16*ORcp9_37-OMcp9_16*ORcp9_38-OMcp9_19*ORcp9_310+OMcp9_36*ORcp9_17+OMcp9_36*ORcp9_18+OMcp9_39*ORcp9_110-OPcp9_16*...
 RLcp9_37-OPcp9_16*RLcp9_38-OPcp9_19*RLcp9_310+OPcp9_36*RLcp9_17+OPcp9_36*RLcp9_18+OPcp9_39*RLcp9_110+ROcp9_86*qdd(7)-(2.0)*qd(7)*(OMcp9_16*ROcp9_96-...
 OMcp9_36*ROcp9_76);
    ACcp9_310 = qdd(3)+OMcp9_16*ORcp9_27+OMcp9_16*ORcp9_28+OMcp9_19*ORcp9_210-OMcp9_26*ORcp9_17-OMcp9_26*ORcp9_18-OMcp9_29*ORcp9_110+OPcp9_16*...
 RLcp9_27+OPcp9_16*RLcp9_28+OPcp9_19*RLcp9_210-OPcp9_26*RLcp9_17-OPcp9_26*RLcp9_18-OPcp9_29*RLcp9_110+ROcp9_96*qdd(7)+(2.0)*qd(7)*(OMcp9_16*ROcp9_86-...
 OMcp9_26*ROcp9_76);

% = = Block_1_0_0_10_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = POcp9_110;
    sens.P(2) = POcp9_210;
    sens.P(3) = POcp9_310;
    sens.R(1,1) = ROcp9_110;
    sens.R(1,2) = ROcp9_210;
    sens.R(1,3) = ROcp9_310;
    sens.R(2,1) = ROcp9_410;
    sens.R(2,2) = ROcp9_510;
    sens.R(2,3) = ROcp9_610;
    sens.R(3,1) = ROcp9_79;
    sens.R(3,2) = ROcp9_89;
    sens.R(3,3) = ROcp9_99;
    sens.V(1) = VIcp9_110;
    sens.V(2) = VIcp9_210;
    sens.V(3) = VIcp9_310;
    sens.OM(1) = OMcp9_110;
    sens.OM(2) = OMcp9_210;
    sens.OM(3) = OMcp9_310;
    sens.A(1) = ACcp9_110;
    sens.A(2) = ACcp9_210;
    sens.A(3) = ACcp9_310;
    sens.OMP(1) = OPcp9_110;
    sens.OMP(2) = OPcp9_210;
    sens.OMP(3) = OPcp9_310;
 
% 
case 11, 


% = = Block_1_0_0_11_0_1 = = 
 
% Sensor Kinematics 


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
    OMcp10_15 = -qd(5)*S4;
    OMcp10_25 = qd(5)*C4;
    OMcp10_16 = OMcp10_15+ROcp10_15*qd(6);
    OMcp10_26 = OMcp10_25+ROcp10_25*qd(6);
    OMcp10_36 = qd(4)-qd(6)*S5;
    OPcp10_16 = ROcp10_15*qdd(6)-qdd(5)*S4-qd(4)*qd(5)*C4-qd(6)*(OMcp10_25*S5+ROcp10_25*qd(4));
    OPcp10_26 = ROcp10_25*qdd(6)+qdd(5)*C4-qd(4)*qd(5)*S4+qd(6)*(OMcp10_15*S5+ROcp10_15*qd(4));
    OPcp10_36 = qdd(4)-qdd(6)*S5-qd(5)*qd(6)*C5;

% = = Block_1_0_0_11_0_2 = = 
 
% Sensor Kinematics 


    RLcp10_17 = Dz73*ROcp10_76;
    RLcp10_27 = Dz73*ROcp10_86;
    RLcp10_37 = Dz73*ROcp10_96;
    ORcp10_17 = OMcp10_26*RLcp10_37-OMcp10_36*RLcp10_27;
    ORcp10_27 = -(OMcp10_16*RLcp10_37-OMcp10_36*RLcp10_17);
    ORcp10_37 = OMcp10_16*RLcp10_27-OMcp10_26*RLcp10_17;

% = = Block_1_0_0_11_0_4 = = 
 
% Sensor Kinematics 


    ROcp10_111 = ROcp10_15*C11-ROcp10_76*S11;
    ROcp10_211 = ROcp10_25*C11-ROcp10_86*S11;
    ROcp10_311 = -(ROcp10_96*S11+C11*S5);
    ROcp10_711 = ROcp10_15*S11+ROcp10_76*C11;
    ROcp10_811 = ROcp10_25*S11+ROcp10_86*C11;
    ROcp10_911 = ROcp10_96*C11-S11*S5;
    RLcp10_111 = ROcp10_46*s.dpt(2,7);
    RLcp10_211 = ROcp10_56*s.dpt(2,7);
    RLcp10_311 = ROcp10_66*s.dpt(2,7);
    POcp10_111 = RLcp10_111+RLcp10_17+q(1);
    POcp10_211 = RLcp10_211+RLcp10_27+q(2);
    POcp10_311 = RLcp10_311+RLcp10_37+q(3);
    OMcp10_111 = OMcp10_16+ROcp10_46*qd(11);
    OMcp10_211 = OMcp10_26+ROcp10_56*qd(11);
    OMcp10_311 = OMcp10_36+ROcp10_66*qd(11);
    ORcp10_111 = OMcp10_26*RLcp10_311-OMcp10_36*RLcp10_211;
    ORcp10_211 = -(OMcp10_16*RLcp10_311-OMcp10_36*RLcp10_111);
    ORcp10_311 = OMcp10_16*RLcp10_211-OMcp10_26*RLcp10_111;
    VIcp10_111 = ORcp10_111+ORcp10_17+qd(1)+ROcp10_76*qd(7);
    VIcp10_211 = ORcp10_211+ORcp10_27+qd(2)+ROcp10_86*qd(7);
    VIcp10_311 = ORcp10_311+ORcp10_37+qd(3)+ROcp10_96*qd(7);
    OPcp10_111 = OPcp10_16+ROcp10_46*qdd(11)+qd(11)*(OMcp10_26*ROcp10_66-OMcp10_36*ROcp10_56);
    OPcp10_211 = OPcp10_26+ROcp10_56*qdd(11)-qd(11)*(OMcp10_16*ROcp10_66-OMcp10_36*ROcp10_46);
    OPcp10_311 = OPcp10_36+ROcp10_66*qdd(11)+qd(11)*(OMcp10_16*ROcp10_56-OMcp10_26*ROcp10_46);
    ACcp10_111 = qdd(1)+OMcp10_26*ORcp10_311+OMcp10_26*ORcp10_37-OMcp10_36*ORcp10_211-OMcp10_36*ORcp10_27+OPcp10_26*RLcp10_311+OPcp10_26*RLcp10_37-...
 OPcp10_36*RLcp10_211-OPcp10_36*RLcp10_27+ROcp10_76*qdd(7)+(2.0)*qd(7)*(OMcp10_26*ROcp10_96-OMcp10_36*ROcp10_86);
    ACcp10_211 = qdd(2)-OMcp10_16*ORcp10_311-OMcp10_16*ORcp10_37+OMcp10_36*ORcp10_111+OMcp10_36*ORcp10_17-OPcp10_16*RLcp10_311-OPcp10_16*RLcp10_37+...
 OPcp10_36*RLcp10_111+OPcp10_36*RLcp10_17+ROcp10_86*qdd(7)-(2.0)*qd(7)*(OMcp10_16*ROcp10_96-OMcp10_36*ROcp10_76);
    ACcp10_311 = qdd(3)+OMcp10_16*ORcp10_211+OMcp10_16*ORcp10_27-OMcp10_26*ORcp10_111-OMcp10_26*ORcp10_17+OPcp10_16*RLcp10_211+OPcp10_16*RLcp10_27-...
 OPcp10_26*RLcp10_111-OPcp10_26*RLcp10_17+ROcp10_96*qdd(7)+(2.0)*qd(7)*(OMcp10_16*ROcp10_86-OMcp10_26*ROcp10_76);

% = = Block_1_0_0_11_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = POcp10_111;
    sens.P(2) = POcp10_211;
    sens.P(3) = POcp10_311;
    sens.R(1,1) = ROcp10_111;
    sens.R(1,2) = ROcp10_211;
    sens.R(1,3) = ROcp10_311;
    sens.R(2,1) = ROcp10_46;
    sens.R(2,2) = ROcp10_56;
    sens.R(2,3) = ROcp10_66;
    sens.R(3,1) = ROcp10_711;
    sens.R(3,2) = ROcp10_811;
    sens.R(3,3) = ROcp10_911;
    sens.V(1) = VIcp10_111;
    sens.V(2) = VIcp10_211;
    sens.V(3) = VIcp10_311;
    sens.OM(1) = OMcp10_111;
    sens.OM(2) = OMcp10_211;
    sens.OM(3) = OMcp10_311;
    sens.A(1) = ACcp10_111;
    sens.A(2) = ACcp10_211;
    sens.A(3) = ACcp10_311;
    sens.OMP(1) = OPcp10_111;
    sens.OMP(2) = OPcp10_211;
    sens.OMP(3) = OPcp10_311;
 
% 
case 12, 


% = = Block_1_0_0_12_0_1 = = 
 
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
    OMcp11_16 = OMcp11_15+ROcp11_15*qd(6);
    OMcp11_26 = OMcp11_25+ROcp11_25*qd(6);
    OMcp11_36 = qd(4)-qd(6)*S5;
    OPcp11_16 = ROcp11_15*qdd(6)-qdd(5)*S4-qd(4)*qd(5)*C4-qd(6)*(OMcp11_25*S5+ROcp11_25*qd(4));
    OPcp11_26 = ROcp11_25*qdd(6)+qdd(5)*C4-qd(4)*qd(5)*S4+qd(6)*(OMcp11_15*S5+ROcp11_15*qd(4));
    OPcp11_36 = qdd(4)-qdd(6)*S5-qd(5)*qd(6)*C5;

% = = Block_1_0_0_12_0_2 = = 
 
% Sensor Kinematics 


    RLcp11_17 = Dz73*ROcp11_76;
    RLcp11_27 = Dz73*ROcp11_86;
    RLcp11_37 = Dz73*ROcp11_96;
    ORcp11_17 = OMcp11_26*RLcp11_37-OMcp11_36*RLcp11_27;
    ORcp11_27 = -(OMcp11_16*RLcp11_37-OMcp11_36*RLcp11_17);
    ORcp11_37 = OMcp11_16*RLcp11_27-OMcp11_26*RLcp11_17;

% = = Block_1_0_0_12_0_4 = = 
 
% Sensor Kinematics 


    ROcp11_111 = ROcp11_15*C11-ROcp11_76*S11;
    ROcp11_211 = ROcp11_25*C11-ROcp11_86*S11;
    ROcp11_311 = -(ROcp11_96*S11+C11*S5);
    ROcp11_711 = ROcp11_15*S11+ROcp11_76*C11;
    ROcp11_811 = ROcp11_25*S11+ROcp11_86*C11;
    ROcp11_911 = ROcp11_96*C11-S11*S5;
    ROcp11_412 = ROcp11_46*C12+ROcp11_711*S12;
    ROcp11_512 = ROcp11_56*C12+ROcp11_811*S12;
    ROcp11_612 = ROcp11_66*C12+ROcp11_911*S12;
    ROcp11_712 = -(ROcp11_46*S12-ROcp11_711*C12);
    ROcp11_812 = -(ROcp11_56*S12-ROcp11_811*C12);
    ROcp11_912 = -(ROcp11_66*S12-ROcp11_911*C12);
    RLcp11_111 = ROcp11_46*s.dpt(2,7);
    RLcp11_211 = ROcp11_56*s.dpt(2,7);
    RLcp11_311 = ROcp11_66*s.dpt(2,7);
    POcp11_111 = RLcp11_111+RLcp11_17+q(1);
    POcp11_211 = RLcp11_211+RLcp11_27+q(2);
    POcp11_311 = RLcp11_311+RLcp11_37+q(3);
    OMcp11_111 = OMcp11_16+ROcp11_46*qd(11);
    OMcp11_211 = OMcp11_26+ROcp11_56*qd(11);
    OMcp11_311 = OMcp11_36+ROcp11_66*qd(11);
    ORcp11_111 = OMcp11_26*RLcp11_311-OMcp11_36*RLcp11_211;
    ORcp11_211 = -(OMcp11_16*RLcp11_311-OMcp11_36*RLcp11_111);
    ORcp11_311 = OMcp11_16*RLcp11_211-OMcp11_26*RLcp11_111;
    VIcp11_111 = ORcp11_111+ORcp11_17+qd(1)+ROcp11_76*qd(7);
    VIcp11_211 = ORcp11_211+ORcp11_27+qd(2)+ROcp11_86*qd(7);
    VIcp11_311 = ORcp11_311+ORcp11_37+qd(3)+ROcp11_96*qd(7);
    ACcp11_111 = qdd(1)+OMcp11_26*ORcp11_311+OMcp11_26*ORcp11_37-OMcp11_36*ORcp11_211-OMcp11_36*ORcp11_27+OPcp11_26*RLcp11_311+OPcp11_26*RLcp11_37-...
 OPcp11_36*RLcp11_211-OPcp11_36*RLcp11_27+ROcp11_76*qdd(7)+(2.0)*qd(7)*(OMcp11_26*ROcp11_96-OMcp11_36*ROcp11_86);
    ACcp11_211 = qdd(2)-OMcp11_16*ORcp11_311-OMcp11_16*ORcp11_37+OMcp11_36*ORcp11_111+OMcp11_36*ORcp11_17-OPcp11_16*RLcp11_311-OPcp11_16*RLcp11_37+...
 OPcp11_36*RLcp11_111+OPcp11_36*RLcp11_17+ROcp11_86*qdd(7)-(2.0)*qd(7)*(OMcp11_16*ROcp11_96-OMcp11_36*ROcp11_76);
    ACcp11_311 = qdd(3)+OMcp11_16*ORcp11_211+OMcp11_16*ORcp11_27-OMcp11_26*ORcp11_111-OMcp11_26*ORcp11_17+OPcp11_16*RLcp11_211+OPcp11_16*RLcp11_27-...
 OPcp11_26*RLcp11_111-OPcp11_26*RLcp11_17+ROcp11_96*qdd(7)+(2.0)*qd(7)*(OMcp11_16*ROcp11_86-OMcp11_26*ROcp11_76);
    OMcp11_112 = OMcp11_111+ROcp11_111*qd(12);
    OMcp11_212 = OMcp11_211+ROcp11_211*qd(12);
    OMcp11_312 = OMcp11_311+ROcp11_311*qd(12);
    OPcp11_112 = OPcp11_16+ROcp11_111*qdd(12)+ROcp11_46*qdd(11)+qd(11)*(OMcp11_26*ROcp11_66-OMcp11_36*ROcp11_56)+qd(12)*(OMcp11_211*ROcp11_311-...
 OMcp11_311*ROcp11_211);
    OPcp11_212 = OPcp11_26+ROcp11_211*qdd(12)+ROcp11_56*qdd(11)-qd(11)*(OMcp11_16*ROcp11_66-OMcp11_36*ROcp11_46)-qd(12)*(OMcp11_111*ROcp11_311-...
 OMcp11_311*ROcp11_111);
    OPcp11_312 = OPcp11_36+ROcp11_311*qdd(12)+ROcp11_66*qdd(11)+qd(11)*(OMcp11_16*ROcp11_56-OMcp11_26*ROcp11_46)+qd(12)*(OMcp11_111*ROcp11_211-...
 OMcp11_211*ROcp11_111);

% = = Block_1_0_0_12_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = POcp11_111;
    sens.P(2) = POcp11_211;
    sens.P(3) = POcp11_311;
    sens.R(1,1) = ROcp11_111;
    sens.R(1,2) = ROcp11_211;
    sens.R(1,3) = ROcp11_311;
    sens.R(2,1) = ROcp11_412;
    sens.R(2,2) = ROcp11_512;
    sens.R(2,3) = ROcp11_612;
    sens.R(3,1) = ROcp11_712;
    sens.R(3,2) = ROcp11_812;
    sens.R(3,3) = ROcp11_912;
    sens.V(1) = VIcp11_111;
    sens.V(2) = VIcp11_211;
    sens.V(3) = VIcp11_311;
    sens.OM(1) = OMcp11_112;
    sens.OM(2) = OMcp11_212;
    sens.OM(3) = OMcp11_312;
    sens.A(1) = ACcp11_111;
    sens.A(2) = ACcp11_211;
    sens.A(3) = ACcp11_311;
    sens.OMP(1) = OPcp11_112;
    sens.OMP(2) = OPcp11_212;
    sens.OMP(3) = OPcp11_312;
 
% 
case 13, 


% = = Block_1_0_0_13_0_1 = = 
 
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
    OMcp12_16 = OMcp12_15+ROcp12_15*qd(6);
    OMcp12_26 = OMcp12_25+ROcp12_25*qd(6);
    OMcp12_36 = qd(4)-qd(6)*S5;
    OPcp12_16 = ROcp12_15*qdd(6)-qdd(5)*S4-qd(4)*qd(5)*C4-qd(6)*(OMcp12_25*S5+ROcp12_25*qd(4));
    OPcp12_26 = ROcp12_25*qdd(6)+qdd(5)*C4-qd(4)*qd(5)*S4+qd(6)*(OMcp12_15*S5+ROcp12_15*qd(4));
    OPcp12_36 = qdd(4)-qdd(6)*S5-qd(5)*qd(6)*C5;

% = = Block_1_0_0_13_0_2 = = 
 
% Sensor Kinematics 


    RLcp12_17 = Dz73*ROcp12_76;
    RLcp12_27 = Dz73*ROcp12_86;
    RLcp12_37 = Dz73*ROcp12_96;
    ORcp12_17 = OMcp12_26*RLcp12_37-OMcp12_36*RLcp12_27;
    ORcp12_27 = -(OMcp12_16*RLcp12_37-OMcp12_36*RLcp12_17);
    ORcp12_37 = OMcp12_16*RLcp12_27-OMcp12_26*RLcp12_17;

% = = Block_1_0_0_13_0_4 = = 
 
% Sensor Kinematics 


    ROcp12_111 = ROcp12_15*C11-ROcp12_76*S11;
    ROcp12_211 = ROcp12_25*C11-ROcp12_86*S11;
    ROcp12_311 = -(ROcp12_96*S11+C11*S5);
    ROcp12_711 = ROcp12_15*S11+ROcp12_76*C11;
    ROcp12_811 = ROcp12_25*S11+ROcp12_86*C11;
    ROcp12_911 = ROcp12_96*C11-S11*S5;
    ROcp12_412 = ROcp12_46*C12+ROcp12_711*S12;
    ROcp12_512 = ROcp12_56*C12+ROcp12_811*S12;
    ROcp12_612 = ROcp12_66*C12+ROcp12_911*S12;
    ROcp12_712 = -(ROcp12_46*S12-ROcp12_711*C12);
    ROcp12_812 = -(ROcp12_56*S12-ROcp12_811*C12);
    ROcp12_912 = -(ROcp12_66*S12-ROcp12_911*C12);
    ROcp12_113 = ROcp12_111*C13+ROcp12_412*S13;
    ROcp12_213 = ROcp12_211*C13+ROcp12_512*S13;
    ROcp12_313 = ROcp12_311*C13+ROcp12_612*S13;
    ROcp12_413 = -(ROcp12_111*S13-ROcp12_412*C13);
    ROcp12_513 = -(ROcp12_211*S13-ROcp12_512*C13);
    ROcp12_613 = -(ROcp12_311*S13-ROcp12_612*C13);
    RLcp12_111 = ROcp12_46*s.dpt(2,7);
    RLcp12_211 = ROcp12_56*s.dpt(2,7);
    RLcp12_311 = ROcp12_66*s.dpt(2,7);
    OMcp12_111 = OMcp12_16+ROcp12_46*qd(11);
    OMcp12_211 = OMcp12_26+ROcp12_56*qd(11);
    OMcp12_311 = OMcp12_36+ROcp12_66*qd(11);
    ORcp12_111 = OMcp12_26*RLcp12_311-OMcp12_36*RLcp12_211;
    ORcp12_211 = -(OMcp12_16*RLcp12_311-OMcp12_36*RLcp12_111);
    ORcp12_311 = OMcp12_16*RLcp12_211-OMcp12_26*RLcp12_111;
    OMcp12_112 = OMcp12_111+ROcp12_111*qd(12);
    OMcp12_212 = OMcp12_211+ROcp12_211*qd(12);
    OMcp12_312 = OMcp12_311+ROcp12_311*qd(12);
    OPcp12_112 = OPcp12_16+ROcp12_111*qdd(12)+ROcp12_46*qdd(11)+qd(11)*(OMcp12_26*ROcp12_66-OMcp12_36*ROcp12_56)+qd(12)*(OMcp12_211*ROcp12_311-...
 OMcp12_311*ROcp12_211);
    OPcp12_212 = OPcp12_26+ROcp12_211*qdd(12)+ROcp12_56*qdd(11)-qd(11)*(OMcp12_16*ROcp12_66-OMcp12_36*ROcp12_46)-qd(12)*(OMcp12_111*ROcp12_311-...
 OMcp12_311*ROcp12_111);
    OPcp12_312 = OPcp12_36+ROcp12_311*qdd(12)+ROcp12_66*qdd(11)+qd(11)*(OMcp12_16*ROcp12_56-OMcp12_26*ROcp12_46)+qd(12)*(OMcp12_111*ROcp12_211-...
 OMcp12_211*ROcp12_111);
    RLcp12_113 = ROcp12_712*s.dpt(3,15);
    RLcp12_213 = ROcp12_812*s.dpt(3,15);
    RLcp12_313 = ROcp12_912*s.dpt(3,15);
    POcp12_113 = RLcp12_111+RLcp12_113+RLcp12_17+q(1);
    POcp12_213 = RLcp12_211+RLcp12_213+RLcp12_27+q(2);
    POcp12_313 = RLcp12_311+RLcp12_313+RLcp12_37+q(3);
    OMcp12_113 = OMcp12_112+ROcp12_712*qd(13);
    OMcp12_213 = OMcp12_212+ROcp12_812*qd(13);
    OMcp12_313 = OMcp12_312+ROcp12_912*qd(13);
    ORcp12_113 = OMcp12_212*RLcp12_313-OMcp12_312*RLcp12_213;
    ORcp12_213 = -(OMcp12_112*RLcp12_313-OMcp12_312*RLcp12_113);
    ORcp12_313 = OMcp12_112*RLcp12_213-OMcp12_212*RLcp12_113;
    VIcp12_113 = ORcp12_111+ORcp12_113+ORcp12_17+qd(1)+ROcp12_76*qd(7);
    VIcp12_213 = ORcp12_211+ORcp12_213+ORcp12_27+qd(2)+ROcp12_86*qd(7);
    VIcp12_313 = ORcp12_311+ORcp12_313+ORcp12_37+qd(3)+ROcp12_96*qd(7);
    OPcp12_113 = OPcp12_112+ROcp12_712*qdd(13)+qd(13)*(OMcp12_212*ROcp12_912-OMcp12_312*ROcp12_812);
    OPcp12_213 = OPcp12_212+ROcp12_812*qdd(13)-qd(13)*(OMcp12_112*ROcp12_912-OMcp12_312*ROcp12_712);
    OPcp12_313 = OPcp12_312+ROcp12_912*qdd(13)+qd(13)*(OMcp12_112*ROcp12_812-OMcp12_212*ROcp12_712);
    ACcp12_113 = qdd(1)+OMcp12_212*ORcp12_313+OMcp12_26*ORcp12_311+OMcp12_26*ORcp12_37-OMcp12_312*ORcp12_213-OMcp12_36*ORcp12_211-OMcp12_36*...
 ORcp12_27+OPcp12_212*RLcp12_313+OPcp12_26*RLcp12_311+OPcp12_26*RLcp12_37-OPcp12_312*RLcp12_213-OPcp12_36*RLcp12_211-OPcp12_36*RLcp12_27+ROcp12_76*...
 qdd(7)+(2.0)*qd(7)*(OMcp12_26*ROcp12_96-OMcp12_36*ROcp12_86);
    ACcp12_213 = qdd(2)-OMcp12_112*ORcp12_313-OMcp12_16*ORcp12_311-OMcp12_16*ORcp12_37+OMcp12_312*ORcp12_113+OMcp12_36*ORcp12_111+OMcp12_36*...
 ORcp12_17-OPcp12_112*RLcp12_313-OPcp12_16*RLcp12_311-OPcp12_16*RLcp12_37+OPcp12_312*RLcp12_113+OPcp12_36*RLcp12_111+OPcp12_36*RLcp12_17+ROcp12_86*...
 qdd(7)-(2.0)*qd(7)*(OMcp12_16*ROcp12_96-OMcp12_36*ROcp12_76);
    ACcp12_313 = qdd(3)+OMcp12_112*ORcp12_213+OMcp12_16*ORcp12_211+OMcp12_16*ORcp12_27-OMcp12_212*ORcp12_113-OMcp12_26*ORcp12_111-OMcp12_26*...
 ORcp12_17+OPcp12_112*RLcp12_213+OPcp12_16*RLcp12_211+OPcp12_16*RLcp12_27-OPcp12_212*RLcp12_113-OPcp12_26*RLcp12_111-OPcp12_26*RLcp12_17+ROcp12_96*...
 qdd(7)+(2.0)*qd(7)*(OMcp12_16*ROcp12_86-OMcp12_26*ROcp12_76);

% = = Block_1_0_0_13_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = POcp12_113;
    sens.P(2) = POcp12_213;
    sens.P(3) = POcp12_313;
    sens.R(1,1) = ROcp12_113;
    sens.R(1,2) = ROcp12_213;
    sens.R(1,3) = ROcp12_313;
    sens.R(2,1) = ROcp12_413;
    sens.R(2,2) = ROcp12_513;
    sens.R(2,3) = ROcp12_613;
    sens.R(3,1) = ROcp12_712;
    sens.R(3,2) = ROcp12_812;
    sens.R(3,3) = ROcp12_912;
    sens.V(1) = VIcp12_113;
    sens.V(2) = VIcp12_213;
    sens.V(3) = VIcp12_313;
    sens.OM(1) = OMcp12_113;
    sens.OM(2) = OMcp12_213;
    sens.OM(3) = OMcp12_313;
    sens.A(1) = ACcp12_113;
    sens.A(2) = ACcp12_213;
    sens.A(3) = ACcp12_313;
    sens.OMP(1) = OPcp12_113;
    sens.OMP(2) = OPcp12_213;
    sens.OMP(3) = OPcp12_313;
 
% 
case 14, 


% = = Block_1_0_0_14_0_1 = = 
 
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
    OMcp13_16 = OMcp13_15+ROcp13_15*qd(6);
    OMcp13_26 = OMcp13_25+ROcp13_25*qd(6);
    OMcp13_36 = qd(4)-qd(6)*S5;
    OPcp13_16 = ROcp13_15*qdd(6)-qdd(5)*S4-qd(4)*qd(5)*C4-qd(6)*(OMcp13_25*S5+ROcp13_25*qd(4));
    OPcp13_26 = ROcp13_25*qdd(6)+qdd(5)*C4-qd(4)*qd(5)*S4+qd(6)*(OMcp13_15*S5+ROcp13_15*qd(4));
    OPcp13_36 = qdd(4)-qdd(6)*S5-qd(5)*qd(6)*C5;

% = = Block_1_0_0_14_0_2 = = 
 
% Sensor Kinematics 


    RLcp13_17 = Dz73*ROcp13_76;
    RLcp13_27 = Dz73*ROcp13_86;
    RLcp13_37 = Dz73*ROcp13_96;
    ORcp13_17 = OMcp13_26*RLcp13_37-OMcp13_36*RLcp13_27;
    ORcp13_27 = -(OMcp13_16*RLcp13_37-OMcp13_36*RLcp13_17);
    ORcp13_37 = OMcp13_16*RLcp13_27-OMcp13_26*RLcp13_17;

% = = Block_1_0_0_14_0_5 = = 
 
% Sensor Kinematics 


    ROcp13_114 = ROcp13_15*C14+ROcp13_46*S14;
    ROcp13_214 = ROcp13_25*C14+ROcp13_56*S14;
    ROcp13_314 = ROcp13_66*S14-C14*S5;
    ROcp13_414 = -(ROcp13_15*S14-ROcp13_46*C14);
    ROcp13_514 = -(ROcp13_25*S14-ROcp13_56*C14);
    ROcp13_614 = ROcp13_66*C14+S14*S5;
    RLcp13_114 = ROcp13_76*s.dpt(3,8);
    RLcp13_214 = ROcp13_86*s.dpt(3,8);
    RLcp13_314 = ROcp13_96*s.dpt(3,8);
    POcp13_114 = RLcp13_114+RLcp13_17+q(1);
    POcp13_214 = RLcp13_214+RLcp13_27+q(2);
    POcp13_314 = RLcp13_314+RLcp13_37+q(3);
    OMcp13_114 = OMcp13_16+ROcp13_76*qd(14);
    OMcp13_214 = OMcp13_26+ROcp13_86*qd(14);
    OMcp13_314 = OMcp13_36+ROcp13_96*qd(14);
    ORcp13_114 = OMcp13_26*RLcp13_314-OMcp13_36*RLcp13_214;
    ORcp13_214 = -(OMcp13_16*RLcp13_314-OMcp13_36*RLcp13_114);
    ORcp13_314 = OMcp13_16*RLcp13_214-OMcp13_26*RLcp13_114;
    VIcp13_114 = ORcp13_114+ORcp13_17+qd(1)+ROcp13_76*qd(7);
    VIcp13_214 = ORcp13_214+ORcp13_27+qd(2)+ROcp13_86*qd(7);
    VIcp13_314 = ORcp13_314+ORcp13_37+qd(3)+ROcp13_96*qd(7);
    OPcp13_114 = OPcp13_16+ROcp13_76*qdd(14)+qd(14)*(OMcp13_26*ROcp13_96-OMcp13_36*ROcp13_86);
    OPcp13_214 = OPcp13_26+ROcp13_86*qdd(14)-qd(14)*(OMcp13_16*ROcp13_96-OMcp13_36*ROcp13_76);
    OPcp13_314 = OPcp13_36+ROcp13_96*qdd(14)+qd(14)*(OMcp13_16*ROcp13_86-OMcp13_26*ROcp13_76);
    ACcp13_114 = qdd(1)+OMcp13_26*ORcp13_314+OMcp13_26*ORcp13_37-OMcp13_36*ORcp13_214-OMcp13_36*ORcp13_27+OPcp13_26*RLcp13_314+OPcp13_26*RLcp13_37-...
 OPcp13_36*RLcp13_214-OPcp13_36*RLcp13_27+ROcp13_76*qdd(7)+(2.0)*qd(7)*(OMcp13_26*ROcp13_96-OMcp13_36*ROcp13_86);
    ACcp13_214 = qdd(2)-OMcp13_16*ORcp13_314-OMcp13_16*ORcp13_37+OMcp13_36*ORcp13_114+OMcp13_36*ORcp13_17-OPcp13_16*RLcp13_314-OPcp13_16*RLcp13_37+...
 OPcp13_36*RLcp13_114+OPcp13_36*RLcp13_17+ROcp13_86*qdd(7)-(2.0)*qd(7)*(OMcp13_16*ROcp13_96-OMcp13_36*ROcp13_76);
    ACcp13_314 = qdd(3)+OMcp13_16*ORcp13_214+OMcp13_16*ORcp13_27-OMcp13_26*ORcp13_114-OMcp13_26*ORcp13_17+OPcp13_16*RLcp13_214+OPcp13_16*RLcp13_27-...
 OPcp13_26*RLcp13_114-OPcp13_26*RLcp13_17+ROcp13_96*qdd(7)+(2.0)*qd(7)*(OMcp13_16*ROcp13_86-OMcp13_26*ROcp13_76);

% = = Block_1_0_0_14_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = POcp13_114;
    sens.P(2) = POcp13_214;
    sens.P(3) = POcp13_314;
    sens.R(1,1) = ROcp13_114;
    sens.R(1,2) = ROcp13_214;
    sens.R(1,3) = ROcp13_314;
    sens.R(2,1) = ROcp13_414;
    sens.R(2,2) = ROcp13_514;
    sens.R(2,3) = ROcp13_614;
    sens.R(3,1) = ROcp13_76;
    sens.R(3,2) = ROcp13_86;
    sens.R(3,3) = ROcp13_96;
    sens.V(1) = VIcp13_114;
    sens.V(2) = VIcp13_214;
    sens.V(3) = VIcp13_314;
    sens.OM(1) = OMcp13_114;
    sens.OM(2) = OMcp13_214;
    sens.OM(3) = OMcp13_314;
    sens.A(1) = ACcp13_114;
    sens.A(2) = ACcp13_214;
    sens.A(3) = ACcp13_314;
    sens.OMP(1) = OPcp13_114;
    sens.OMP(2) = OPcp13_214;
    sens.OMP(3) = OPcp13_314;
 
% 
case 15, 


% = = Block_1_0_0_15_0_1 = = 
 
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
    OMcp14_16 = OMcp14_15+ROcp14_15*qd(6);
    OMcp14_26 = OMcp14_25+ROcp14_25*qd(6);
    OMcp14_36 = qd(4)-qd(6)*S5;
    OPcp14_16 = ROcp14_15*qdd(6)-qdd(5)*S4-qd(4)*qd(5)*C4-qd(6)*(OMcp14_25*S5+ROcp14_25*qd(4));
    OPcp14_26 = ROcp14_25*qdd(6)+qdd(5)*C4-qd(4)*qd(5)*S4+qd(6)*(OMcp14_15*S5+ROcp14_15*qd(4));
    OPcp14_36 = qdd(4)-qdd(6)*S5-qd(5)*qd(6)*C5;

% = = Block_1_0_0_15_0_2 = = 
 
% Sensor Kinematics 


    RLcp14_17 = Dz73*ROcp14_76;
    RLcp14_27 = Dz73*ROcp14_86;
    RLcp14_37 = Dz73*ROcp14_96;
    ORcp14_17 = OMcp14_26*RLcp14_37-OMcp14_36*RLcp14_27;
    ORcp14_27 = -(OMcp14_16*RLcp14_37-OMcp14_36*RLcp14_17);
    ORcp14_37 = OMcp14_16*RLcp14_27-OMcp14_26*RLcp14_17;

% = = Block_1_0_0_15_0_5 = = 
 
% Sensor Kinematics 


    ROcp14_114 = ROcp14_15*C14+ROcp14_46*S14;
    ROcp14_214 = ROcp14_25*C14+ROcp14_56*S14;
    ROcp14_314 = ROcp14_66*S14-C14*S5;
    ROcp14_414 = -(ROcp14_15*S14-ROcp14_46*C14);
    ROcp14_514 = -(ROcp14_25*S14-ROcp14_56*C14);
    ROcp14_614 = ROcp14_66*C14+S14*S5;
    ROcp14_115 = ROcp14_114*C15-ROcp14_76*S15;
    ROcp14_215 = ROcp14_214*C15-ROcp14_86*S15;
    ROcp14_315 = ROcp14_314*C15-ROcp14_96*S15;
    ROcp14_715 = ROcp14_114*S15+ROcp14_76*C15;
    ROcp14_815 = ROcp14_214*S15+ROcp14_86*C15;
    ROcp14_915 = ROcp14_314*S15+ROcp14_96*C15;
    RLcp14_114 = ROcp14_76*s.dpt(3,8);
    RLcp14_214 = ROcp14_86*s.dpt(3,8);
    RLcp14_314 = ROcp14_96*s.dpt(3,8);
    POcp14_114 = RLcp14_114+RLcp14_17+q(1);
    POcp14_214 = RLcp14_214+RLcp14_27+q(2);
    POcp14_314 = RLcp14_314+RLcp14_37+q(3);
    OMcp14_114 = OMcp14_16+ROcp14_76*qd(14);
    OMcp14_214 = OMcp14_26+ROcp14_86*qd(14);
    OMcp14_314 = OMcp14_36+ROcp14_96*qd(14);
    ORcp14_114 = OMcp14_26*RLcp14_314-OMcp14_36*RLcp14_214;
    ORcp14_214 = -(OMcp14_16*RLcp14_314-OMcp14_36*RLcp14_114);
    ORcp14_314 = OMcp14_16*RLcp14_214-OMcp14_26*RLcp14_114;
    VIcp14_114 = ORcp14_114+ORcp14_17+qd(1)+ROcp14_76*qd(7);
    VIcp14_214 = ORcp14_214+ORcp14_27+qd(2)+ROcp14_86*qd(7);
    VIcp14_314 = ORcp14_314+ORcp14_37+qd(3)+ROcp14_96*qd(7);
    ACcp14_114 = qdd(1)+OMcp14_26*ORcp14_314+OMcp14_26*ORcp14_37-OMcp14_36*ORcp14_214-OMcp14_36*ORcp14_27+OPcp14_26*RLcp14_314+OPcp14_26*RLcp14_37-...
 OPcp14_36*RLcp14_214-OPcp14_36*RLcp14_27+ROcp14_76*qdd(7)+(2.0)*qd(7)*(OMcp14_26*ROcp14_96-OMcp14_36*ROcp14_86);
    ACcp14_214 = qdd(2)-OMcp14_16*ORcp14_314-OMcp14_16*ORcp14_37+OMcp14_36*ORcp14_114+OMcp14_36*ORcp14_17-OPcp14_16*RLcp14_314-OPcp14_16*RLcp14_37+...
 OPcp14_36*RLcp14_114+OPcp14_36*RLcp14_17+ROcp14_86*qdd(7)-(2.0)*qd(7)*(OMcp14_16*ROcp14_96-OMcp14_36*ROcp14_76);
    ACcp14_314 = qdd(3)+OMcp14_16*ORcp14_214+OMcp14_16*ORcp14_27-OMcp14_26*ORcp14_114-OMcp14_26*ORcp14_17+OPcp14_16*RLcp14_214+OPcp14_16*RLcp14_27-...
 OPcp14_26*RLcp14_114-OPcp14_26*RLcp14_17+ROcp14_96*qdd(7)+(2.0)*qd(7)*(OMcp14_16*ROcp14_86-OMcp14_26*ROcp14_76);
    OMcp14_115 = OMcp14_114+ROcp14_414*qd(15);
    OMcp14_215 = OMcp14_214+ROcp14_514*qd(15);
    OMcp14_315 = OMcp14_314+ROcp14_614*qd(15);
    OPcp14_115 = OPcp14_16+ROcp14_414*qdd(15)+ROcp14_76*qdd(14)+qd(14)*(OMcp14_26*ROcp14_96-OMcp14_36*ROcp14_86)+qd(15)*(OMcp14_214*ROcp14_614-...
 OMcp14_314*ROcp14_514);
    OPcp14_215 = OPcp14_26+ROcp14_514*qdd(15)+ROcp14_86*qdd(14)-qd(14)*(OMcp14_16*ROcp14_96-OMcp14_36*ROcp14_76)-qd(15)*(OMcp14_114*ROcp14_614-...
 OMcp14_314*ROcp14_414);
    OPcp14_315 = OPcp14_36+ROcp14_614*qdd(15)+ROcp14_96*qdd(14)+qd(14)*(OMcp14_16*ROcp14_86-OMcp14_26*ROcp14_76)+qd(15)*(OMcp14_114*ROcp14_514-...
 OMcp14_214*ROcp14_414);

% = = Block_1_0_0_15_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = POcp14_114;
    sens.P(2) = POcp14_214;
    sens.P(3) = POcp14_314;
    sens.R(1,1) = ROcp14_115;
    sens.R(1,2) = ROcp14_215;
    sens.R(1,3) = ROcp14_315;
    sens.R(2,1) = ROcp14_414;
    sens.R(2,2) = ROcp14_514;
    sens.R(2,3) = ROcp14_614;
    sens.R(3,1) = ROcp14_715;
    sens.R(3,2) = ROcp14_815;
    sens.R(3,3) = ROcp14_915;
    sens.V(1) = VIcp14_114;
    sens.V(2) = VIcp14_214;
    sens.V(3) = VIcp14_314;
    sens.OM(1) = OMcp14_115;
    sens.OM(2) = OMcp14_215;
    sens.OM(3) = OMcp14_315;
    sens.A(1) = ACcp14_114;
    sens.A(2) = ACcp14_214;
    sens.A(3) = ACcp14_314;
    sens.OMP(1) = OPcp14_115;
    sens.OMP(2) = OPcp14_215;
    sens.OMP(3) = OPcp14_315;
 
% 
case 16, 


% = = Block_1_0_0_16_0_1 = = 
 
% Sensor Kinematics 


    ROcp15_15 = C4*C5;
    ROcp15_25 = S4*C5;
    ROcp15_75 = C4*S5;
    ROcp15_85 = S4*S5;
    ROcp15_46 = ROcp15_75*S6-S4*C6;
    ROcp15_56 = ROcp15_85*S6+C4*C6;
    ROcp15_66 = C5*S6;
    ROcp15_76 = ROcp15_75*C6+S4*S6;
    ROcp15_86 = ROcp15_85*C6-C4*S6;
    ROcp15_96 = C5*C6;
    OMcp15_15 = -qd(5)*S4;
    OMcp15_25 = qd(5)*C4;
    OMcp15_16 = OMcp15_15+ROcp15_15*qd(6);
    OMcp15_26 = OMcp15_25+ROcp15_25*qd(6);
    OMcp15_36 = qd(4)-qd(6)*S5;
    OPcp15_16 = ROcp15_15*qdd(6)-qdd(5)*S4-qd(4)*qd(5)*C4-qd(6)*(OMcp15_25*S5+ROcp15_25*qd(4));
    OPcp15_26 = ROcp15_25*qdd(6)+qdd(5)*C4-qd(4)*qd(5)*S4+qd(6)*(OMcp15_15*S5+ROcp15_15*qd(4));
    OPcp15_36 = qdd(4)-qdd(6)*S5-qd(5)*qd(6)*C5;

% = = Block_1_0_0_16_0_6 = = 
 
% Sensor Kinematics 


    ROcp15_116 = ROcp15_15*C16-ROcp15_76*S16;
    ROcp15_216 = ROcp15_25*C16-ROcp15_86*S16;
    ROcp15_316 = -(ROcp15_96*S16+C16*S5);
    ROcp15_716 = ROcp15_15*S16+ROcp15_76*C16;
    ROcp15_816 = ROcp15_25*S16+ROcp15_86*C16;
    ROcp15_916 = ROcp15_96*C16-S16*S5;
    RLcp15_116 = ROcp15_46*s.dpt(2,3);
    RLcp15_216 = ROcp15_56*s.dpt(2,3);
    RLcp15_316 = ROcp15_66*s.dpt(2,3);
    POcp15_116 = RLcp15_116+q(1);
    POcp15_216 = RLcp15_216+q(2);
    POcp15_316 = RLcp15_316+q(3);
    OMcp15_116 = OMcp15_16+ROcp15_46*qd(16);
    OMcp15_216 = OMcp15_26+ROcp15_56*qd(16);
    OMcp15_316 = OMcp15_36+ROcp15_66*qd(16);
    ORcp15_116 = OMcp15_26*RLcp15_316-OMcp15_36*RLcp15_216;
    ORcp15_216 = -(OMcp15_16*RLcp15_316-OMcp15_36*RLcp15_116);
    ORcp15_316 = OMcp15_16*RLcp15_216-OMcp15_26*RLcp15_116;
    VIcp15_116 = ORcp15_116+qd(1);
    VIcp15_216 = ORcp15_216+qd(2);
    VIcp15_316 = ORcp15_316+qd(3);
    OPcp15_116 = OPcp15_16+ROcp15_46*qdd(16)+qd(16)*(OMcp15_26*ROcp15_66-OMcp15_36*ROcp15_56);
    OPcp15_216 = OPcp15_26+ROcp15_56*qdd(16)-qd(16)*(OMcp15_16*ROcp15_66-OMcp15_36*ROcp15_46);
    OPcp15_316 = OPcp15_36+ROcp15_66*qdd(16)+qd(16)*(OMcp15_16*ROcp15_56-OMcp15_26*ROcp15_46);
    ACcp15_116 = qdd(1)+OMcp15_26*ORcp15_316-OMcp15_36*ORcp15_216+OPcp15_26*RLcp15_316-OPcp15_36*RLcp15_216;
    ACcp15_216 = qdd(2)-OMcp15_16*ORcp15_316+OMcp15_36*ORcp15_116-OPcp15_16*RLcp15_316+OPcp15_36*RLcp15_116;
    ACcp15_316 = qdd(3)+OMcp15_16*ORcp15_216-OMcp15_26*ORcp15_116+OPcp15_16*RLcp15_216-OPcp15_26*RLcp15_116;

% = = Block_1_0_0_16_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = POcp15_116;
    sens.P(2) = POcp15_216;
    sens.P(3) = POcp15_316;
    sens.R(1,1) = ROcp15_116;
    sens.R(1,2) = ROcp15_216;
    sens.R(1,3) = ROcp15_316;
    sens.R(2,1) = ROcp15_46;
    sens.R(2,2) = ROcp15_56;
    sens.R(2,3) = ROcp15_66;
    sens.R(3,1) = ROcp15_716;
    sens.R(3,2) = ROcp15_816;
    sens.R(3,3) = ROcp15_916;
    sens.V(1) = VIcp15_116;
    sens.V(2) = VIcp15_216;
    sens.V(3) = VIcp15_316;
    sens.OM(1) = OMcp15_116;
    sens.OM(2) = OMcp15_216;
    sens.OM(3) = OMcp15_316;
    sens.A(1) = ACcp15_116;
    sens.A(2) = ACcp15_216;
    sens.A(3) = ACcp15_316;
    sens.OMP(1) = OPcp15_116;
    sens.OMP(2) = OPcp15_216;
    sens.OMP(3) = OPcp15_316;
 
% 
case 17, 


% = = Block_1_0_0_17_0_1 = = 
 
% Sensor Kinematics 


    ROcp16_15 = C4*C5;
    ROcp16_25 = S4*C5;
    ROcp16_75 = C4*S5;
    ROcp16_85 = S4*S5;
    ROcp16_46 = ROcp16_75*S6-S4*C6;
    ROcp16_56 = ROcp16_85*S6+C4*C6;
    ROcp16_66 = C5*S6;
    ROcp16_76 = ROcp16_75*C6+S4*S6;
    ROcp16_86 = ROcp16_85*C6-C4*S6;
    ROcp16_96 = C5*C6;
    OMcp16_15 = -qd(5)*S4;
    OMcp16_25 = qd(5)*C4;
    OMcp16_16 = OMcp16_15+ROcp16_15*qd(6);
    OMcp16_26 = OMcp16_25+ROcp16_25*qd(6);
    OMcp16_36 = qd(4)-qd(6)*S5;
    OPcp16_16 = ROcp16_15*qdd(6)-qdd(5)*S4-qd(4)*qd(5)*C4-qd(6)*(OMcp16_25*S5+ROcp16_25*qd(4));
    OPcp16_26 = ROcp16_25*qdd(6)+qdd(5)*C4-qd(4)*qd(5)*S4+qd(6)*(OMcp16_15*S5+ROcp16_15*qd(4));
    OPcp16_36 = qdd(4)-qdd(6)*S5-qd(5)*qd(6)*C5;

% = = Block_1_0_0_17_0_7 = = 
 
% Sensor Kinematics 


    ROcp16_117 = ROcp16_15*C17-ROcp16_76*S17;
    ROcp16_217 = ROcp16_25*C17-ROcp16_86*S17;
    ROcp16_317 = -(ROcp16_96*S17+C17*S5);
    ROcp16_717 = ROcp16_15*S17+ROcp16_76*C17;
    ROcp16_817 = ROcp16_25*S17+ROcp16_86*C17;
    ROcp16_917 = ROcp16_96*C17-S17*S5;
    RLcp16_117 = ROcp16_46*s.dpt(2,4);
    RLcp16_217 = ROcp16_56*s.dpt(2,4);
    RLcp16_317 = ROcp16_66*s.dpt(2,4);
    POcp16_117 = RLcp16_117+q(1);
    POcp16_217 = RLcp16_217+q(2);
    POcp16_317 = RLcp16_317+q(3);
    OMcp16_117 = OMcp16_16+ROcp16_46*qd(17);
    OMcp16_217 = OMcp16_26+ROcp16_56*qd(17);
    OMcp16_317 = OMcp16_36+ROcp16_66*qd(17);
    ORcp16_117 = OMcp16_26*RLcp16_317-OMcp16_36*RLcp16_217;
    ORcp16_217 = -(OMcp16_16*RLcp16_317-OMcp16_36*RLcp16_117);
    ORcp16_317 = OMcp16_16*RLcp16_217-OMcp16_26*RLcp16_117;
    VIcp16_117 = ORcp16_117+qd(1);
    VIcp16_217 = ORcp16_217+qd(2);
    VIcp16_317 = ORcp16_317+qd(3);
    OPcp16_117 = OPcp16_16+ROcp16_46*qdd(17)+qd(17)*(OMcp16_26*ROcp16_66-OMcp16_36*ROcp16_56);
    OPcp16_217 = OPcp16_26+ROcp16_56*qdd(17)-qd(17)*(OMcp16_16*ROcp16_66-OMcp16_36*ROcp16_46);
    OPcp16_317 = OPcp16_36+ROcp16_66*qdd(17)+qd(17)*(OMcp16_16*ROcp16_56-OMcp16_26*ROcp16_46);
    ACcp16_117 = qdd(1)+OMcp16_26*ORcp16_317-OMcp16_36*ORcp16_217+OPcp16_26*RLcp16_317-OPcp16_36*RLcp16_217;
    ACcp16_217 = qdd(2)-OMcp16_16*ORcp16_317+OMcp16_36*ORcp16_117-OPcp16_16*RLcp16_317+OPcp16_36*RLcp16_117;
    ACcp16_317 = qdd(3)+OMcp16_16*ORcp16_217-OMcp16_26*ORcp16_117+OPcp16_16*RLcp16_217-OPcp16_26*RLcp16_117;

% = = Block_1_0_0_17_1_0 = = 
 
% Symbolic Outputs  

    sens.P(1) = POcp16_117;
    sens.P(2) = POcp16_217;
    sens.P(3) = POcp16_317;
    sens.R(1,1) = ROcp16_117;
    sens.R(1,2) = ROcp16_217;
    sens.R(1,3) = ROcp16_317;
    sens.R(2,1) = ROcp16_46;
    sens.R(2,2) = ROcp16_56;
    sens.R(2,3) = ROcp16_66;
    sens.R(3,1) = ROcp16_717;
    sens.R(3,2) = ROcp16_817;
    sens.R(3,3) = ROcp16_917;
    sens.V(1) = VIcp16_117;
    sens.V(2) = VIcp16_217;
    sens.V(3) = VIcp16_317;
    sens.OM(1) = OMcp16_117;
    sens.OM(2) = OMcp16_217;
    sens.OM(3) = OMcp16_317;
    sens.A(1) = ACcp16_117;
    sens.A(2) = ACcp16_217;
    sens.A(3) = ACcp16_317;
    sens.OMP(1) = OPcp16_117;
    sens.OMP(2) = OPcp16_217;
    sens.OMP(3) = OPcp16_317;

end


% ====== END Task 1 ====== 

  

