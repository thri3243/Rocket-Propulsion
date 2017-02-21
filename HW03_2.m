%%%%%%%%%%%%%%%
% Tad Riley
% HW03
% ASEN 5053 Rocket Propulsion
% Dr. Lakshmi Kantha
%%%%%%%%%%%%%%%

g=9.8067
p_sl=101325;

% Question 1
mi=3000;    	%[kg], Initial mass
mf=300;         %[kg], final mass
tb=45;          %[s], burn time
F_avg=120e3;    %[N], average thrust
p_c=8.00e6;     %[Pa], chamber pressure
p_e=.08e6;      %[Pa], Nozzle exit pressure
D_noz_t=.12;    %[m], Nozzle throat diameter
D_noz_e=.4;     %[m], Nozzle exit diameter

%mass flow rate [kg/s]
mdot=(mi-mf)/tb

%exit velocity [m/s]
ve=F_avg/mdot

cstar=(p_c*(pi*(D_noz_t/2)^2))/mdot
I_sp=F_avg/(mdot*g)
%Effective velocity
c_sl=ve+((p_e-p_sl)*((D_noz_e/2)^2*pi))/mdot
p_30km=1186;
c_30km=ve+((p_e-p_30km)*((D_noz_e/2)^2*pi))/mdot

F_30km=F_avg+(p_e-p_30km)*((D_noz_e/2)^2*pi)
C_F_sl=F_avg/(p_c*pi*(D_noz_t/2)^2)
C_F_30km=F_30km/(p_c*pi*(D_noz_t/2)^2)

%%
%Question 2
I_sp_2=420;     %[s] Specific impulse of all stages

mi1s=100000;     %[kg] initial mass
mpay1s=1000;     %[kg] Payload mass
mprop1s=75000;   %[kg] propellant mass

mi2s=18000;    %[kg] initial mass of 2nd stage
mprop2s=12000; %[kg] Mass of 2nd stage propellant

mi3s=4000;     %[kg] 3rd stage initial mass
mprop3s=2000;  %[kg] 3rd stage porpellant mass

%part a
mf1s=mi1s-mprop1s;
minert1s=mf1s-mpay1s;
DelV1s=I_sp_2*g*log(mi1s/mf1s)

mf2s=mi2s-mprop2s;
minert2s=mf2s-mpay1s;
DelV2s=I_sp_2*g*log(mi2s/mf2s)

mf3s=mi3s-mprop3s;
minert3s=mf3s-mpay1s;
DelV3s=I_sp_2*g*log(mi3s/mf3s)

%part b
f_pay_b=mpay1s/(mpay1s+minert1s)
f_inert_b=minert1s/(mprop1s+minert1s)

MR=(1+f_pay_b)/(f_inert_b+f_pay_b)
DelV_tot_b=3*I_sp_2*g*log(MR)
DelV_tot_a=DelV1s+DelV2s+DelV3s
diff_DelV=DelV_tot_b-DelV_tot_a

%part c
%%%%%%%
%%%%
%
%
%
%
%
%
%
%
%
%
%%
%Question 3

m_pay3=100000;      %[kg] Payload mass
DelV_Mina=13200;     %[m/s] DelV plus loss required
f_inert_3=.06;      %Inert mass fraction
I_sp3avg=400;       %[s] average Isp
Na=3;                %Number of stages

%Part a, 3 stage optimum rocket

mi3a=m_pay3*((exp(DelV_Mina/(I_sp3avg*g*Na))*(1-f_inert_3))/(1-f_inert_3*exp(DelV_Mina/(I_sp3avg*g*Na))))^Na

%Part b, refuel in LEO => 2 Stage to orbit
DelV_Minb=9600;     %[m/s] DelV plus loss required
Nb=2;

mi3b=m_pay3*((exp(DelV_Minb/(I_sp3avg*g*Nb))*(1-f_inert_3))/(1-f_inert_3*exp(DelV_Minb/(I_sp3avg*g*Nb))))^Nb

DelV_Minc=3600;     %[m/s] DelV plus loss required
Nc=1;
mi3c=m_pay3*((exp(DelV_Minc/(I_sp3avg*g*Nc))*(1-f_inert_3))/(1-f_inert_3*exp(DelV_Minc/(I_sp3avg*g*Nc))))^Nc

%Mass reduction to LEO using LEO refueling
advantage=100-((mi3b/mi3a)*100)

mLEO=mi3c-m_pay3;
minertLEO=mLEO*f_inert_3;

%Mass propellant required in LEO for additional DelV
mpropLEO=mLEO-minertLEO






