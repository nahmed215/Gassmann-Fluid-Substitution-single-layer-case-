close all, clear all, 
clc,
%% INPUT Parameters for FSM, these are derived from well logs (Well_16_1_13_Composite)

% load FSM_data.mat D d T P Vcl Vp1 Vs1 Rho1 Vp2 Vs2 Rho2 PhiT PhiE Sw; 

D = 1917; d = 50; Vp2 = 2.923694; Vs2 = 1.613505; Rho2 = 2.169049;
T = 78.186244; P = 26.0; PhiT = 0.264803; PhiE = 0.253118;
Sw = 0.333908; Vcl = 0.125913; Kcl = 21; Kqrt = 37; rh_cl = 2600; 
rh_qrt = 2650; S = 0.003; API = 35; SG = 0.6;


%% Function OUTPUTS Ksat, VP, VS, rh_eff after Gassmann FSM.
[Ksat, K_Gass,K_GH,VP_GH,VP,VS,rh_eff,K_wet,K_hyd,Kd,Ks] = ahmedFSM(T, P, S, PhiE, Vcl, Kcl, Kqrt, API,...
    SG, Sw, rh_cl, rh_qrt, Vs2, Rho2);

Vp1 = 3.835670; Vs1 = 1.957859; Rho1 = 2.418910; % Upper medium properties
fprintf('VP1 = %f\n', Vp1); fprintf('VS1 = %f\n', Vs1); fprintf('Rho1 = %f\n', Rho1);
fprintf('VP2 = %f\n', Vp2); fprintf('VS2 = %f\n', Vs2); fprintf('Rho2 = %f\n', Rho2);
fprintf('Gass VP = %f\n', VP); fprintf('Gass VS = %f\n', VS); fprintf('Gass Rho = %f\n', rh_eff);

fprintf('Patchy VP = %f\n', VP_GH); 

