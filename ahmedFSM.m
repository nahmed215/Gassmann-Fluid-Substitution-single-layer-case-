function [Ksat, K_Gass,K_GH,VP_GH,VP,VS,rh_eff,K_wet,K_hyd,Kd,Ks] = ahmedFSM(T, P, S, PhiE, Vcl, Kcl, Kqrt,...
    API, SG, Sw, rh_cl, rh_qrt, Vs2, Rho2)
% function[Ksat, VP, VS, rh_eff] = ahmedFSM(T, P, S, PhiE, Vcl, Kcl, Kqrt, API,...
%    SG, Sw, rh_cl, rh_qrt, Vs2, Rho2)

%   This function retuns saturated rock bulk modulus (Ksat), P and S waves
%   velocities (VP & VS) and effective density as outputs after performing
%   Gassmann-1951 fluid substitution at desired fluid saturation (oil/brine)
%   levels. The inputs required for this modeling are tempreature (T in C),
%   pressure (P in MPa), salinty (S in ppm), effective porosity (PhiE in
%   frac.), clay volume (Vcl in frac.), bulk modulus of clay (Kcl in GPa),
%   bulk modulus of quartz (Kqrt in GPa), API gravity of oil, specific
%   gravity of gas dissolved in oil (SG), initial water saturation (Sw in
%   frac.), density of clay mineral (rh_cl in kg/m³), density of quartz (
%   rh_qrt in kg/m³), S wave velocity from logs (Vs2 in Km/s) and logs bulk
%   density of rock (Rho2 in g/cm3).
% Note: This code only consider two fluids (oil/brine) and rock matrix is
% comprises on quartz and clay minerals.

% The other mathematical notations are:
%  Kd = dry rock modulus or Bulk modulus of rock skeleton (Murphy's Eq)
%  Vqrt = volume fraction of quartz 
%  Kv = Voigt bulk modulus of rock matrix
%  Kr = Reuss bulk modulus of rock matrix
%  Ks = Voigt-Reuss-Hill average bulk modulus of rock matrix
%  alf = alfa is Biot coefficient
%  r_w = density of water at given in-situ conditions (T&P)
%  r_b = density of brine (after including saline content) at given T & P
%  r_br = % brine density in Kg/m3
%  vp_w = P wave velcity in water at stock tank conditions
%  vp_b = P wave velcity in water at stock tank conditions
%  Kw = bulk modulus of brine (fluid) 
%  r_o = reference density at 15.6 C
%  r_p = density at pressure P
%  rh_o = in-situ density at T was developed by Dodson and Standing (1945)
%  Vp_o = oil veocity in terms of T & P
%  Vpo = oil veocity in terms of API
%  RG = GOR Gas-oil ratio in terms of T & P
%  Rg = GOR Gas-oil ratio in terms of API
%  Ko = bulk modulus of oil
%  So = saturation of oil (1-Sw)
%  Kfl = Bulk modulus of fluid (oil/brine)
%  rh_fl = density of fluid (oil/brine)
%  rh_mat = density of matrix/rock forming minerals
%  rh_eff = effective density of composite rock unit
%  mu = shear modulus 
%  Ksat = saturated rock bulk modulus via Gassmann's Eq 1951
%  VP = Gassmann VP, after FSM
%  VS = Gassmann VS, after FSM

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
rh_o = (r_p/(0.972+3.81*10^-4*(T+17.78)^1.175))*10^3                       % Batzle & Wang 1992 Eq. 19 (density in Kg/cm3)
Vp_o = 2096*(r_o/(2.6-r_o))^0.5-3.7*T+4.64*P+0.0115*(4.12*...
    (1.08*r_o^-1-1)^0.5-1)*T*P;                                             % Batzle & Wang 1992 Eq. 20a (VP based on r_o and T&P)
Vpo = 15450*(77.1+API)^-0.5-3.7*T+4.64*P+0.0115*(0.36*API^0.5-1)*T*P;       % Batzle & Wang 1992 Eq. 20b (VP based on API and T&P)
RG = 0.02123*SG*(P*exp((4.072/r_o)-0.00377*T))^1.205;                       % Batzle & Wang 1992 Eq. 21a (GOR based on r_o and T&P)
Rg = 2.03*SG*(P*exp(0.02878*API-0.00377*T))^1.205;                          % Batzle & Wang 1992 Eq. 21b (GOR based on API and T&P)
Ko = rh_o*Vpo^2*10^-9;
So = 1-Sw; 
Kfl = 1/((Sw/Kw)+(So/Ko));                                                  % Wood's Equation fluid modulus
rh_fl = r_br*Sw+rh_o*So;
rh_mat = Vcl*rh_cl+Vqrt*rh_qrt;
rh_eff = rh_mat*(1-PhiE)+PhiE*rh_fl;
Vs2 =Vs2*1000; Rho2 = Rho2*1000;
mu = (Vs2^2*Rho2)/10^9;

% Full Gassmann Equation

K_Gass = Kd + ((1-Kd/Ks)^2/((PhiE/Kfl)+(1-PhiE/Ks)+(Kd/Ks^2)));

%% approximated Gassmann eq. version
M = ((PhiE/Kfl)+((alf-PhiE)/Ks))^(-1);
Ksat = Kd+(alf^2)*M;                                                % Gassmann Equation for Saturated Rock Bulk Modulus
VP = sqrt((Ksat*10^9+(4/3)*mu*10^9)/rh_eff); VP = VP/1000;
VS = sqrt(mu*10^9/rh_eff); VS = VS/1000; rh_eff = rh_eff/1000;

% For Patchy saturation case with approximated Gassmann eq K_GH

M1 = (((alf-PhiE)./Ks)+(PhiE./Kw)).^(-1);
M2 = (((alf-PhiE)./Ks)+(PhiE./Ko)).^(-1);
Kwet = Kd+(alf.^2.*M1);
Khyd = Kd+(alf.^2.*M2);
K_GH = ((Sw./(Kwet+(4/3).*mu))+((1-Sw)./(Khyd+(4/3).*mu))).^(-1)-(4/3).*mu;

% patchy seismic velocity VP with approximated Gassmann equation
Rho_FSM1 = rh_eff.*1000;
VP_GH = sqrt((K_GH.*10^9+(4/3).*mu.*10^9)./Rho_FSM1); 
VP_GH = VP_GH./1000;

%% test GPR 2018, 66 Wollner and Dvorkinn CASE2 PATCHY

K_wet = Ks*((PhiE*Kd-((1-PhiE)*Kw*Kd/Ks)+Kw)/((1-PhiE)*Kw+PhiE*Ks-(Kw*Kd/Ks)));
K_hyd = Ks*((PhiE*Kd-((1-PhiE)*Ko*Kd/Ks)+Ko)/((1-PhiE)*Ko+PhiE*Ks-(Ko*Kd/Ks)));

K_GH2 = ((Sw./(K_wet+(4/3).*mu))+((1-Sw)./(K_hyd+(4/3).*mu))).^(-1)-(4/3).*mu;

VP_GH2 = sqrt((K_GH2.*10^9+(4/3).*mu.*10^9)./Rho_FSM1); 
VP_GH2 = VP_GH2./1000;

end 


