%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tad Riley
% ASEN 5053 Rocket Propulsion
% Dr. Lakshmi Kantha
% HW07
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
close all

g=9.806;
%Problem 1

m_N2=10;        %[kg]
m_H2=10;        %[kg]
m_He=15;        %[kg]
P_He=6.7e6;     %[Pa]
T_He=300;       %[k]

MH2=2.02;
MN2=28.01;
MHe=4.00;

etaH2=m_H2/MH2;
etaN2=m_N2/MN2;
etaHe=m_He/MHe;

eta=etaH2+etaN2+etaHe;

chi_He=etaHe/eta;
chi_N2=etaN2/eta;
chi_H2=etaH2/eta;

Mbar=MH2*chi_H2+MN2*chi_N2+MHe*chi_He
%Find M_Bar of the mixture and gamma and c_p of each species

%Problem 2

%moles of each constituent
nH2=.317;
nO2=.101;
nH2O=1.512;
nH=.109;
nO=.054;
nOH=.233;

%cp of each constituent [kJ/kmol*K]
cpH2=.0329;
cpO2=.0363;
cpH2O=.0483;
cpH=.0209;
cpO=.0209;
cpOH=.0335;

%heat of formation of each constituent [kJ/kmol]
hfH2=0;
hfO2=0;
hfH2O=-241.83;
hfH=218.0;
hfO=294.17;
hfOH=39.46;

T_rxn=250;
cpH2_rxn=.0296;
cpO2_rxn=.0323;

%Heat of Reaction
H_rxn=(nH2*hfH2+nO2*hfO2+nH2O*hfH2O+nO*hfO+nH*hfH+nOH*hfOH)-(2*hfH2+hfO2)

%Adiabatic Flame Temperature
H_products=(nH2*hfH2+nO2*hfO2+nH2O*hfH2O+nO*hfO+nH*hfH+nOH*hfOH);
H_reactants=(2*hfH2+hfO2);
T_f=(-H_products+(nO2*cpO2+nH2*cpH2+nH2O*cpH2O+nO*cpO+nH*cpH+nOH*cpOH)*T_rxn)/(nO2*cpO2+nH2*cpH2+nH2O*cpH2O+nO*cpO+nH*cpH+nOH*cpOH)

%Problem 3
%MMH/N2O4 hypergolic engine

ObyF=1.75;       %oxidizer to fuel ratio
L_star=.75;         %characteristic length

m_payl=4914;        %[kg]
DelV_desire=1721;   %[m/s]
TbyW=.3;            %thrust to weight
m_i=12000;    %[kg] minimum initial mass
L=3;                %[m] max length
D=3;                %[m] max diameter
P_c=.7e6;           %[Pa] combustion pressure

DelV_design=1.1*DelV_desire    %[m/s] designed with 10% additional delV for security
T_design=TbyW*g*m_i
%From page 3 appendix B
gamma=1.245;
Isp_vac=305;
T_c=3000;
M_bar=20;

%Engine characteristics
m_engine=(T_design/g)/(0.0006098*T_design+13.44)
l_engine=(.0054*T_design+31.92)/100
D_engine=(.00357*T_design+14.48)/100
TbyW_engine=T_design/(m_engine*g)

eps=100;                %Nozzle expansion ratio
P_e=.00058*P_c          %Exit pressure
P_a=0;                  %Ambient pressure
c_star=1790.4;
C_F=1.892;              %optimal
Isp=(C_F*c_star)/g
C_F2=(Isp_vac*g)/c_star     %given Isp from NASA

rho_MMH=874;
rho_NTO=1431;
u=7;                    %[m/s] flow velocity
%propellant tank
P_prop_tank=1.2*P_c+5e5+.5*rho_MMH*u^2
%oxidizer tank
P_ox_tank=1.2*P_c+5e5+.5*rho_NTO*u^2

%Pressurize to 1.5 MPa
%Pressurant helium tank will be 21 MPa

v_e=Isp_vac*g
m_f=m_i*exp(-DelV_design/v_e)
m_inert=m_f-m_payl
m_prop=m_i-m_f
m_fuel=m_prop/(ObyF+1)
m_ox=m_prop-m_fuel

f_inert=m_inert/(m_inert+m_prop)
V_fuel=1.1*(m_fuel/rho_MMH)
V_ox=1.1*(m_ox/rho_NTO)

mdot=T_design/v_e
mdot_f=mdot/(1+ObyF)
mdot_ox=(mdot*ObyF)/(1+ObyF)

t_b=m_prop/mdot

%Nozzle and thrust chamber
A_t=(mdot*c_star)/P_c
D_t=sqrt((A_t*4)/pi)
A_e=A_t*100
D_e=sqrt((A_e*4)/pi)

%15 degree nozzle
NL=(D_e-D_t)/(2*tand(15))
eta_Noz=.5*(1+cosd(15))
%Bell nozzle
NL_bell=.675*NL
m_noz=.002*NL_bell*(24/(.002*2.1))

%thrust chamber
M_c=.2;
A_c=3*A_t
D_c=sqrt((A_c*4)/pi)
L_c=L_star/3
t_c=2*P_c*1.5*(D_c/2)/310e6
m_c=(.5*pi*D_c^2+pi*D_c*L_c)*t_c*8500
m_eng=(m_c+m_noz)/.4
m_tvc=m_engine-m_eng
m_injector=m_eng*.25
m_ins=m_eng*.35

%propellant tanks continued
m_f_tank=10^6*2*V_fuel/(2500*g)
m_ox_tank=10^6*2*V_ox/(2500*g)

%Pressurant of He at 21 MPA
MB_He=4;
gam_He=1.66;
T_a=273;
P_ini=21e6;
P_fin=1.5e6;

T_fin=T_a*(P_fin/P_ini)^((gam_He-1)/gam_He)
V_He=(T_fin/T_a)*(P_ini/P_fin)
V_He_ini=(V_fuel+V_ox)/(V_He-1)
m_He_loaded=(MB_He*P_ini*V_He_ini)/(8314.41*T_a)
m_He_tank=(P_ini*V_He_ini)/(g*6350)

%toroidal tanks
Rr_f=V_fuel/(2*pi^2)
Rr_ox=V_ox/(2*pi^2)
Rr_He=V_He/(2*pi^2)

rfuel=.3346;
rox=.3480;
rHe=.4979;

Rfuel=Rr_f/rfuel^2
Rox=Rr_ox/rox^2
RHe=Rr_He/rHe^2