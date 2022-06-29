close all, clear all, 
clc,
%% INPUT Parameters for FSM, these are derived from well logs (Well_16_1_13_Composite)

% load FSM_data.mat D d T P Vcl Vp1 Vs1 Rho1 Vp2 Vs2 Rho2 PhiT PhiE Sw; 

D = 1917; d = 50; Vp2 = 2.923694; Vs2 = 1.613505; Rho2 = 2.169049;
T = 78.186244; P = 26.0; PhiT = 0.264803; PhiE = 0.253118;
Vcl = 0.125913; Kcl = 21; Kqrt = 37; rh_cl = 2600; 
rh_qrt = 2650; S = 0.003; API = 35; SG = 0.6; Mu_cl = 7; Mu_qrt = 44;


% %% Function OUTPUTS Ksat, VP, VS, rh_eff after Gassmann FSM.
% [Ksat, K_Gass,K_GH,VP_GH,VP,VS,rh_eff,K_wet,K_hyd,Kd,Ks] = ahmedFSM_sw_iter(T, P, S, PhiE, Vcl, Kcl, Kqrt, API,...
%     SG, rh_cl, rh_qrt, Vs2, Rho2);

Kd = 38.18*(1-3.39*PhiE+1.95*PhiE^2); 
Vqrt = 1-Vcl;
Kv = Kqrt*Vqrt+Kcl*Vcl;
Kr = 1/((Vqrt/Kqrt)+(Vcl/Kcl));
Ks = (Kv+Kr)/2;                      
alf = 1-Kd/Ks;
% Fluid properties for 'brine fluid' starts here
r_w = 1+1*10^-6*(-80*T-3.3*T^2+0.00175*T^3+489*P-2*T*P+0.016*T^2*P-1.3*...  % water density
    10^-5*T^3*P-0.333*P^2-0.002*T*P^2);
r_b = r_w+S*(0.668+0.44*S+1*10^-6*(300*P-2400*P*S+T*(80+3*T-3300*...        % brine density
    S-13*P+47*P*S)));
r_br = r_b*1000;

w(1,1) = 1402.85;                  w(1,3) = 3.4370*10^(-3);   
w(2,1) = 4.871;                    w(2,3) = 1.7390*10^(-4);
w(3,1) = -0.04783;                 w(3,3) = -2.135*10^(-6);
w(4,1) = 1.487*10^(-4);            w(4,3) = -1.455*10^(-8);
w(5,1) = -2.197*10^(-7);           w(5,3) = 5.230*10^(-11);

w(1,2) = 1.524;                    w(1,4) = -1.197*10^(-5);
w(2,2) = -0.0111;                  w(2,4) = -1.628*10^(-6);
w(3,2) = 2.747*10^(-4);            w(3,4) = 1.2370*10^(-8);
w(4,2) = -6.503*10^(-7);           w(4,4) = 1.327*10^(-10);
w(5,2) = 7.987*10^(-10);           w(5,4) = -4.614*10^(-13);
sum=0;
for i = 1:5
    for j = 1:4
        sum = sum+w(i,j)*T^(i-1)*P^(j-1); 
    end
end
vp_w = sum;                                                                 % water velocity

vp_b = vp_w+S*(1170-9.6*T+0.055*T^2-8.5*10^-5*T^3+2.6*P-0.0029*T*...        % Brine velocity
    P-0.476*P^2)+S^1.5*(780-10*P+0.16*P^2)-1820*S^2;
Kw = vp_b^2*r_br/10^9;
% Fluid properties for 'oil fluid' starts here
r_o = 141.5/(API+131.5);                                                    % Batzle & Wang 1992 Eq. 14 
r_p = r_o+(0.00277*P-1.71*10^-7*P^3)*(r_o-1.15)^2+3.49*10^-4*P;             % Batzle & Wang 1992 Eq. 18
rh_o = (r_p/(0.972+3.81*10^-4*(T+17.78)^1.175))*10^3;                       % Batzle & Wang 1992 Eq. 19 (density in Kg/cm3)
Vp_o = 2096*(r_o/(2.6-r_o))^0.5-3.7*T+4.64*P+0.0115*(4.12*...
    (1.08*r_o^-1-1)^0.5-1)*T*P;                                             % Batzle & Wang 1992 Eq. 20a (VP based on r_o and T&P)
Vpo = 15450*(77.1+API)^-0.5-3.7*T+4.64*P+0.0115*(0.36*API^0.5-1)*T*P;       % Batzle & Wang 1992 Eq. 20b (VP based on API and T&P)
RG = 0.02123*SG*(P*exp((4.072/r_o)-0.00377*T))^1.205;                       % Batzle & Wang 1992 Eq. 21a (GOR based on r_o and T&P)
Rg = 2.03*SG*(P*exp(0.02878*API-0.00377*T))^1.205;                          % Batzle & Wang 1992 Eq. 21b (GOR based on API and T&P)
Ko = rh_o*Vpo^2*10^-9;

Sw_i = 0;
inter = 0.1;
Sw_f = 1;
Sw = Sw_i:inter:Sw_f;

for ii = 0:1
    
  So = 1-Sw;
  
Kfl = 1./((Sw./Kw)+(So./Ko));                                                  % Wood's Equation fluid modulus
rh_fl = r_br.*Sw+rh_o.*So;
rh_mat = Vcl*rh_cl+Vqrt*rh_qrt;
rh_eff = rh_mat*(1-PhiE)+PhiE.*rh_fl;
Vs2 =Vs2*1000; 
Rho2 = Rho2*1000;
mu = (Vs2^2*Rho2)/10^9;

M = ((PhiE./Kfl)+((alf-PhiE)/Ks)).^(-1);
Ksat = Kd+(alf^2).*M;                                                % Gassmann Equation for Saturated Rock Bulk Modulus
VP = sqrt((Ksat.*10^9+(4/3).*mu)./rh_eff); 
VS = sqrt(mu./rh_eff);
rh_eff = rh_eff./1000;
  
  
end

figure, plot(Sw,VP./1000,'r','displayname', 'V_P (km/s)','linewidth', 1.5); hold on; 
plot(Sw,VS./1000,'--b','displayname','V_S (km/s)','linewidth', 1.5);
plot(Sw,rh_eff,'-.k','displayname','rho (g/cc)','linewidth', 1.5);
title('Elastic properties vs. saturation');
xlabel('S_w (v/v)'); ylabel('Elastic properties'); 
xlim([0 1]); xticks([0:0.2:1]); yticks([1.4:0.4:3.2]);
legend

%% For Effective stress

Pk = 33.82; % Max pressure in MPA % Depeth = 1917 m = 6289 ft, Fracture gradient = 1.8 gr/cc = 0.78 psi/ft, Fracture Pressure = 0.780*6289 = 4905.42 psi = 33.8216 MPa
Ek = 0.1;
P_mu = 40.82; E_mu =1.8;
Mu_s = Mu_qrt.*Vqrt + Mu_cl.*Vcl;
Cphi = 0.40; 

Sw1 = 0.33;
So1 = 1-Sw1;
Kfl_p = 1/((Sw1/Kw)+(So1/Ko));  

rh_flp = r_br*Sw1+rh_o*So1;
rh_effp = rh_mat*(1-PhiE)+PhiE.*rh_flp;

VES_i = 22;
int = 0.5;
VES_f = 28;

VES = VES_i:int:VES_f;

for jj = 22:28
    
Kdry = (Ks.*(1-PhiE./Cphi))./(1+Ek.*exp(-VES./Pk)); % option 1
M_dry = (Mu_s.*(1-PhiE./Cphi))./(1+E_mu.*exp(-VES./P_mu));   

alf2 = 1-Kdry./Ks;
M = ((PhiE./Kfl_p)+((alf2-PhiE)./Ks)).^(-1);
Ksat_P = Kdry+(alf2.^2).*M;                                                % Gassmann Equation for Saturated Rock Bulk Modulus
VP = sqrt((Ksat_P.*10^9+(4/3).*M_dry)./rh_effp); 
VS = sqrt(M_dry.*10^9./rh_effp);


end


figure, subplot(211); plot(VES,VP./1000,'r','displayname', 'V_P (km/s)','linewidth', 1.5);
title('V_P vs. effective pressure');
xlabel('P_eff (MPa)'); ylabel('P wave velocity (km/s'); 
xlim([22 28]); xticks([22:2.28]);
legend

subplot(212); plot(VES,VS./1000,'--b','displayname','V_S (km/s)','linewidth', 1.5);
title('V_S vs. effective pressure');
xlabel('P_eff (MPa)'); ylabel('S wave velocity (km/s)'); 
xlim([22 28]); xticks([22:2.28]); 
legend


