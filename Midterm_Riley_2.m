
% -------------------------------------------------------------------------
%              ULA Delta II 7920 Performance Estimation
%
%       Data from spaceflight101.com supplemented and validated by
%       Data from SpaceX website, Wikipedia,astronautix.com
%
%                     Copyright Mr. T Riley
%
% Note that the analysis is based on imprecise and often unreliable public
% domain data and may also invove some "guestimates" for important numbers
% Also, a different point in the parameter space may yield a similar result
% So this should be regardes as an illustration of the analysis technique
% Any results must be verified by actual data from the Launchers, which is
% most often not available. Copyright is simply to avoid any legal issues.
% -------------------------------------------------------------------------

clc; 
clear all; close all;

global F A Cd tturn hturn Vturn Ang0 mdot RE g0 nstage steer V Ang
global Vrot ii  t X h gls dls

g0 = 9.806; d2r = pi/180; RE = 6378e3; 

% ------------------------------
% Delta II 7920 Specifications
% ------------------------------

% Type	Delta II 7920	
% Height	38.9m	
% Diameter	2.44m	
% Launch Pad Mass	228,000kg
% Stages	2	
% Boosters	9	
% Mass to LEO	4.850kg	
% Mass to SSO	3,000kg
% Takes __ hrs to fuel
 
% First Stage

% Type:	Extra Long extended Tank Thor	
% Length:	26.1 m	
% Diameter:	2.44 m
% Fuel:	Rocket Propellant 1	(RP-1) 
% Oxidizer:	Liquid Oxygen	 (LOX)
% Inert Mass:	5,680 kg	
% Propellant Mass:	96,120 kg	
% Launch Mass:    101800 kg
% Propellant Tanks:	Aluminum Isogrid
% Tank Pressurization:	Helium
% Propulsion:	(1) RS-27A	  
% Engine Type:	Gas Generator, Open Cycle
% RS-27A Thrust:	Sea Level- 890kN , Vac- 1054kN
% Specific Impulse:	Sea Level- 254s , Vac- 302s
% Throttle Range:	100% Only (No throttling)
% Engine Diameter:	1.07m
% Engine Length:    3.78 m
% Engine Dry Mass:	1147 kg	
% Thrust to Weight:   102.5
% Chamber Pressure:	4.8MPa (696 psi)
% Expansion Ratio:	12	 
% Mixture Ratio: 2.25
% Burn Time:	261s	 
% Guidance:	From 2nd Stage	
% Attitude Control:	Gimbaled Engines (pitch, yaw) 
%                   (2) LR-101 Vernier Engines (roll)
% Thrust: 4.4 kN (pump-fed) 3.7 kN (tank-fed)
% Type: Sustainer Engines
% Mixture Ratio: 1.8
% Area Ratio: 5.6
% Specific Impulse: 207s (pump-fed) 197s(Tank-fed)
% Chamber Pressure: 2.47 MPa (pump-fed) 2.08 MPa (Tank-fed)
% Propellant Pressure: 4.34 MPa (Pump-fed) 3.52 MPa (Tank-fed)
% RP-1 Flowrate: .76 kg/s (Pump) .68 kg/s (Tank)
% LOX Flowrate: 1.40 kg/s (Pump 1.22 kg/s (Tank)
% Engine Mass:  25 kg

% Solid Rocket Motors

% Type: GEM-40
% # Boosters: 9
% Ground-lit: 6
% Air-lit: 3
% Length: 11.05 m
% Diameter: 1.03 m
% Propellant: HTPB
% Inert Mass: 1102 kg
% Propellant Mass: 11766 kg
% Launch Mass: 12962 kg
% Propulsion: GEM-40
% Engine Type: Solid Rocket Motor
% Maximum Thrust: 643.8 kN
% Average Thrust: 499.1 kN
% Specific Impulse: 245s (SL) 274s (Vacuum)
% Nozzle Diameter .82 m
% Chamber Pressure: 5.6 MPa
% Area Ratio: 11
% Burn Time: 63.3 s
% Vehicle Control: None

% Second Stage

% Type:	Delta-K
% Length:	5.97 m	
% Diameter: 2.44m	
% Fuel: Aerozine-50	
% Oxidizer: Nitrogen Tertroxide
% Inert Mass:   950 kg
% Propellant Mass:  6000 kg
% Launch Mass:  6950 kg 
% Propellant Tanks: Aluminum Isogrid	 
% Tank Struture:    Common Blukhead 
% Propulsion:	AJ10-118K
% Engine Type:	Pressure Fed	
% Thrust:   43.4 kN (Vacuum)
% Specific Impulse:   319s
% Engine Diameter:	1.7 m	
% Engine Mass:	98 kg	
% Burn Time:	<500 s	
% Chamber Pressure:	8.9 MPa (1290 psi)	
% Restart Capability:	Yes	
% Engine Cooling: Ablative
% Fuel Flow Rate: 4.75 kg/s
% Oxidizer Flow Rate:    9.10 kg/s
% Mixture Ratio:  1.9
% Expansion Ratio:    65 
% Attitude Control:	Reaction Control System	 

% Payload Fairing
 
% Type:   Composite Fairing
% Diameter:   3 m
% Length: 8.88 m
% Mass: 3524 kg 
% 
% Mach 1 at 33 s, maxQ at 45 s, 2 air lit solids 1m 12s, Mach 5 at 1m 55s, 
% Mach 10 at 2m 55s, MECO at 265 s
% ------------------------------------------------------------------------
%               DeltaV estimate using Ideal Rocket Equation
% ------------------------------------------------------------------------
% launched from Vandenburg AFB
phi = 34.7; alt = 10; ii=98.7;   % Launch Site Specs

%86164 s in a sidereal day
theta = acosd(cosd(ii)/sind(phi));                  % Orbit Inclination in deg
Vrot = ((2*pi*RE)/86164)*cosd(theta)*sind(phi);     % Earth's rotation in m/s

% ------------------------------------------------------------------------
%         PART 1: October 28, 2011 Suommi NPP Mission
% ------------------------------------------------------------------------

m_payl = 2128;  % 1400 kg dry, 464 kg payload 

%  Second Stage

Isp2_vac = 319;                  % Vacuum
T2_vac = 43.4e3;                 % AJ10-118K Engine

mdot2 = T2_vac/(Isp2_vac*g0);    % kg/s
m_prop2 = 6000;                  % Propellant mass (2nd stage)
tb2 = m_prop2/mdot2;
m_inrt2 = 950;                   % Inert mass (2nd stage) (Includes Engine??? - 98 kg)   

mi2 = m_payl+m_inrt2+m_prop2;
mf2 = m_payl+m_inrt2;
MR2 = mi2/mf2;
DelV2 = Isp2_vac*g0*log(MR2);  % Second stage DeltaV
m_tot2 = mi2;
disp([' MR2 m_tot2 = ',num2str(MR2),'   ',num2str(m_tot2)])

%  First Stage Core

t_IG_C = -2.0;              % Ignition
t_LO =  0;                % Lift-off time

m_propC = 96120 ;         % First stage prop mass total
m_inrtC = 5680;           % First Stage inert mass (9 engines - 4500 kg) 
m_fairing = 800;        %unknown ath this point!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IspC_sl = 254;                      
TC_sl = 890.0e3;                     
IspC_vac = 302;
TC_vac = 1054.0e3;

mdotC = TC_sl/(IspC_sl*g0);
% mdotC = TC_vac/(IspC_vac*g0);

% First Stage Boosters, 6 grount lit, 3 air lit

t_IG_GB=t_LO;       %ignition of groung lit motors
t_IG_AB=63.3;       %ignition of air lit motors

t_BO_GB=63.3;       %Burnout of ground lit boosters
t_BO_AB=126.6;      %Burnout of airlit motors

t_J_GB=86;          %Jettison of ground lit boosters
t_J_AB=128;         %Jettison of airlit boosters

m_propB=11766;      %Propellant mass per booster
m_inrtB=1102;       %Inert mass of booster

IspB_sl=245;        %Specific impulse of booster at sea level
TB_sl=499.1e3;      %Average sea level thrust for boosters
IspB_vac=274;

mdotB=TB_sl/(IspB_sl*g0);


m_tot1 = m_inrtC + m_propC + 9*(m_propB+m_inrtB);

% Performance (assume constant thrust) 

m_LP = m_tot1 + m_tot2 + m_fairing;
m_burnLO = mdotC*(t_LO-t_IG_C);
m_LO = m_LP - m_burnLO;

TbyW = (TC_sl+6*TB_sl)/(m_LO*g0);

% At Jettison of ground lit boosters
Isp_J_GB=0.5*(IspC_sl+IspC_vac)+0.5*(IspB_sl+IspB_vac);
m_J_GB=m_LP-(6*(m_propB+m_inrtB))-t_J_GB*mdotC-(t_J_GB-t_IG_AB)*3*mdotB;
MR_J_GB=m_LO/m_J_GB;
DelV_J_GB=Isp_J_GB*g0*log(MR_J_GB);

% At Jettison of air lit boosters
Isp_J_AB=Isp_J_GB;
m_J_AB=m_LP-(9*(m_propB+m_inrtB))-t_J_AB*mdotC;
MR_J_AB=m_J_GB/m_J_AB;
DelV_J_AB=Isp_J_AB*g0*log(MR_J_AB);

% At main engine cutoff

Isp1_MECO = 0.5*(IspC_sl+IspC_vac);
m_MECO = m_LP- m_propC-(9*(m_propB+m_inrtB)); 
MR_MECO = m_J_AB/m_MECO;
DelV1 = Isp1_MECO*g0*log(MR_MECO);

%How do i adjust this properly?????????????????????????????????????????????

disp([' MR1 m_tot1 m_LP m_LO = ',num2str(MR_MECO),'   ',num2str(m_tot1), ...
      '   ',num2str(m_LP),'   ',num2str(m_LO)])
fprintf('\n');

DelV = DelV1 + DelV2 +Vrot + DelV_J_GB + DelV_J_AB;

disp([' T/W = ',num2str(TbyW)])
disp([' DelV_J_GB = ',num2str(DelV_J_GB),' m/s'])
disp([' DelV_J_AB = ',num2str(DelV_J_AB),' m/s'])
disp([' DelV1 = ',num2str(DelV1),' m/s'])
disp([' DelV2 = ',num2str(DelV2),' m/s'])
disp([' DelV_rot = ',num2str(fix(Vrot)),' m/s'])
disp([' DelV  = ',num2str(DelV),' m/s'])

mu=398600;
hp=208; ha=365;                    % ii = 98.7 deg;
Inc=ii;
Rp=hp+RE/1000; Ra=ha+RE/1000; a=0.5*(Rp+Ra); % all in km
Vp=sqrt(mu)*sqrt((2/Rp)-(1/a));
Va=Vp*Rp/Ra;
np=sqrt(mu/a^3);
Tp=(2*pi/np)/60;     % Orbit period in minutes
DelV_orbit=( sqrt( (Vp^2) + 2*hp*(g0/1000)*(0.001*RE/(0.001*RE+hp))^2 ) )*1000;

DelV_loss=DelV-DelV_orbit;

disp([' V orbit = ',num2str(fix(0.5*(Vp+Va)*1000)),' m/s'])
disp([' DelV orbit = ',num2str(DelV_orbit),' m/s'])
disp([' DelV loss (actual?) = ',num2str(DelV_loss),' m/s'])

% ------------------------------------------------------------------------
%      Part 2: Calculate the Flight Trajectory of the Rocket
% ------------------------------------------------------------------------

D = 3.66; A = pi*D^2/4; Cd = 1.0;  % Assume constant Cd

t_MECO = 165.0;  % MECO
t_SEP1 = 168.0;  % First stage separation
t_SIG = 170.0;   % Second stage ignition
t_FJ = 220.0;    % Fairing Jettison
t_SECO = 560.0;  % Second stage engine cutoff

t_turn = 23;     % gravity turn  13     
Ang0 = 88.8;     % initial angle of rocket
steer1 = +0.000;    % First stage steering
steer2 = 0.00125;    % Second stage steering

mi_corrn=0; ntry=1;
% mi1=505537-mi_corrn;  ????

                  for n=1:ntry

            % First stage

nstage=1;
steer = steer1;

mi1=m_LP-mi_corrn;
F1_sl = TC_sl;
Isp1_sl = IspC_sl;
F1_vac = TC_vac;
Isp1_vac = IspC_vac;
mdot1 = F1_sl/(Isp1_sl*g0);
tdelay=t_LO-t_IG_C;
mi1  = mi1-tdelay*mdot1;         
minrt1 = m_inrtC;
mprop1=m_propC;
mf1 = mi1 - mprop1;
MR1 = mi1/mf1;
Isp1_MECO = 0.5*(Isp1_sl+Isp1_vac);
F1_MECO = 0.5*(F1_sl+F1_vac);
DelV1 = Isp1_MECO*g0*log(MR1);
tb1 = mprop1/mdot1;
% tb1 = 164;

ti1 = 0; tf1 = ti1 + t_turn;
t_range1 = [ti1,tf1];

F = F1_MECO; mdot=mdot1;
mi=mi1;

u1=0; v1=0; X1=0; h1=alt; gls1=0; dls1=0;
statei1 = [u1, v1, X1, h1, gls1, dls1, mi];
[t,state] = ode45(@Rocket2,t_range1,statei1);

u     = state(:,1);
v     = state(:,2);        
X     = state(:,3);      
h     = state(:,4);      
gls   = state(:,5);          
dls   = state(:,6);          
m     = state(:,7);

           if (n==ntry)
plotout
fprintf('\n Beginning of Gravity Turn')
fprintf('\n Final speed               = %5.1f (m/s)',V(end))
fprintf('\n Final angle               = %5.1f (deg)',Ang(end)*180/pi)
fprintf('\n Altitude                  = %5.1f (km)',h(end)/1000)
fprintf('\n Downrange                 = %5.1f (km)',X(end)/1000)
fprintf('\n Gravity loss              = %5.1f (m/s)',-gls(end))
fprintf('\n Drag loss                 = %5.1f (m/s)',-dls(end))
mf=m(end);
fprintf('\n First  Stage mi mf        = %5.0f %5.0f(kg)',fix(mi),fix(mf))
fprintf('\n');
           end
           
tig=tf1; tfg=tb1+t_IG_C;
t_rangeg = [tig,tfg];

F = F1_MECO; mdot=mdot1;
mi=m(end);

Vig=V(end);Angig=Ang0*pi/180;Xig=X(end);hig=h(end);glsig=gls(end);dlsig=dls(end);
stateig = [Vig, Angig, Xig, hig, glsig, dlsig, mi];
[t,state] = ode45(@Rocket,t_rangeg,stateig);

V     = state(:,1);      
Ang   = state(:,2);           
X     = state(:,3);      
h     = state(:,4);
gls   = state(:,5);          
dls   = state(:,6);          
m     = state(:,7);  

           if (n==ntry)
plotout
fprintf('\n End of First Stage Burn')
fprintf('\n Final speed               = %5.1f (m/s)',V(end))
fprintf('\n Final angle               = %5.1f (deg)',Ang(end)*180/pi)
fprintf('\n Altitude                  = %5.1f (km)',h(end)/1000)
fprintf('\n Downrange                 = %5.1f (km)',X(end)/1000)
fprintf('\n Gravity loss              = %5.1f (m/s)',-gls(end))
fprintf('\n Drag loss                 = %5.1f (m/s)',-dls(end))
mf=m(end);
fprintf('\n First  Stage mi mf        = %5.0f %5.0f(kg)',fix(mi),fix(mf))
fprintf('\n DelV1                     = %5.1f (m/s)',DelV1)
fprintf('\n');

[DR1,HA1]=Downrange(V(end),Ang(end),h(end));
fprintf('\n First Stage Impact DR     = %5.0f (km)',fix(DR1))
fprintf('\n First Stage Apogee        = %5.0f (km)',fix(HA1))
fprintf('\n');
           end

                % Second stage
nstage=2;
steer = steer2;

mprop2 = m_prop2;   
minrt2 = m_inrt2;
mi2  = m_LP-m_tot1; 
mf2 = mi2 - mprop2;
MR2 = mi2/mf2;
F2 = T2_vac;
Isp2 = Isp2_vac;
DelV2 = Isp2*g0*log(MR2);
mdot2 = F2/(Isp2*g0);
tb2 = mprop2/mdot2;          
% tb2 = 350;

tij = tfg; tfj = t_FJ;
t_rangej = [tij,tfj];

F = F2; mdot=mdot2;
mi=mi2;

Vij=V(end);Angij=Ang(end);Xij=X(end);hij=h(end);glsij=gls(end);dlsij=dls(end);
stateij = [Vij, Angij, Xij, hij, glsij, dlsij, mi];
[t,state] = ode45(@Rocket,t_rangej,stateij);

V     = state(:,1);      
Ang   = state(:,2);           
X     = state(:,3);      
h     = state(:,4);
gls   = state(:,5);          
dls   = state(:,6);          
m     = state(:,7);

          if (n==ntry)
plotout
fprintf('\n Fairing Jettison')
fprintf('\n Final speed               = %5.1f (m/s)',V(end))
fprintf('\n Final angle               = %5.1f (deg)',Ang(end)*180/pi)
fprintf('\n Altitude                  = %5.1f (km)',h(end)/1000)
fprintf('\n Downrange                 = %5.1f (km)',X(end)/1000)
fprintf('\n Gravity loss              = %5.1f (m/s)',-gls(end))
fprintf('\n Drag loss                 = %5.1f (m/s)',-dls(end))
mf=m(end);
fprintf('\n Second  Stage mi mf       = %5.0f %5.0f(kg)',fix(mi),fix(mf))

[DRf, HAf]=Downrange(V(end),Ang(end),h(end));
fprintf('\n');
fprintf('\n Fairing Impact DR         = %5.0f (km)',fix(DRf))
fprintf('\n Fairing Apogee            = %5.0f (km)',fix(HAf))
fprintf('\n');
          end

ti2 = tfj; tf2 = tfg + tb2;
t_range2 = [ti2,tf2];

F = F2; mdot=mdot2;
mi=m(end)-m_fairing;

Vi2=V(end);Angi2=Ang(end);Xi2=X(end);hi2=h(end);glsi2=gls(end);dlsi2=dls(end);
statei2 = [Vi2, Angi2, Xi2, hi2, glsi2, dlsi2, mi];
[t,state] = ode45(@Rocket,t_range2,statei2);

V     = state(:,1);      
Ang   = state(:,2);           
X     = state(:,3);      
h     = state(:,4);
gls   = state(:,5);          
dls   = state(:,6);          
m     = state(:,7);  

          if (n==ntry)
plotout
fprintf('\n End of Second Stage Burn')
fprintf('\n Final speed               = %5.1f (m/s)',V(end))
fprintf('\n Final angle               = %5.1f (deg)',Ang(end)*180/pi)
fprintf('\n Altitude                  = %5.1f (km)',h(end)/1000)
fprintf('\n Downrange                 = %5.1f (km)',X(end)/1000)
fprintf('\n Gravity loss              = %5.1f (m/s)',-gls(end))
fprintf('\n Drag loss                 = %5.1f (m/s)',-dls(end))
mf=m(end);
fprintf('\n Second  Stage mi mf       = %5.0f %5.0f(kg)',fix(mi),fix(mf))
fprintf('\n DelV2                     = %5.1f (m/s)',DelV2)
fprintf('\n');

[DR2,HA2]=Downrange(V(end),Ang(end),h(end));
fprintf('\n Second Stage Impact DR    = %5.0f (km)',fix(DR2))
fprintf('\n Second Stage Apogee       = %5.0f (km)',fix(HA2))
fprintf('\n');

DelV = DelV1 + DelV2 + Vrot;
fprintf('\n DelV                      = %5.1f (km/s)',DelV)
          
          end
          
mfinal=mf2-minrt2;
mi_corrn = mfinal-m_payl;

         if (n==ntry)
fprintf('\n');
fprintf('\n     Altitude            = %5.1f (km)',h(end)/1000)
fprintf('\n     Spacecraft speed    = %5.1f (m/s)',V(end))
fprintf('\n     V orbit             = %5.1f (m/s)',fix(0.5*(Vp+Va)*1000))
fprintf('\n     Total gravity loss  = %5.1f (m/s)',-gls(end))
fprintf('\n');
fprintf('\n     Total drag loss     = %5.1f (m/s)',-dls(end))
fprintf('\n     mpayl               = %5.1f (kg)',mfinal)
fprintf('\n     Launchpad Mass      = %5.1f (kg)',m_LP)
fprintf('\n');
         end       

                  end
                  
% ------------------------------------------------------------------------
%     PART 3: Compute drag and gravity losses using mission profile
% ------------------------------------------------------------------------

D = 3.66; A = pi*D^2/4; Cd = 1.0;           % Constant Cd assumption

t=[0  60  120 150 180 240 300 435 540 560]; % sec  % SpaceX MECO at 163s, Sig at 168 s SECO at 570s?
                             % time is from ignition, not lift-off; 
h=[0 5.5  35  58  90  134 166 205 208 210]; % km
h=h*1000;
x=[0 0.9 14.5 33  67  135 225 495 850 900]; % km
x=x*1000 + Vrot*cosd(ii)*t;
v=[0 0.25 1.0 1.8 1.9 2.1 2.6 4.2 6.9 7.8]; % km/s     Incorrect
v=[0 0.25 1.0 1.8 1.9 2.3 3.0 4.8 7.6 7.8]; % km/s     Corrected
v=v*1000;

% MECO-165; St1_sep-168

mm=length(t);
maxt=max(t);
ind=zeros(size(t));
gls=zeros(size(t));
dls=zeros(size(t));
M=zeros(size(t));
veq=zeros(size(t));
for m=1:mm
    veq(m)=sqrt( (v(m)^2) + 2*h(m)*g0*(RE/(RE+0.001*h(m)))^2 );
end

tt=[0:1:maxt];
for m=1:mm
    ind(m)=find(tt==t(m));
end

hh=zeros(size(tt));
xx=zeros(size(tt));
ang=zeros(size(tt));

hh=interp1(t,h,tt,'pchip');
xx=interp1(t,x,tt,'pchip');
vv=interp1(t,v,tt,'pchip');

                       figure(2)
lwd=2; font=20;
Xmin=0; Xdel=200; Xmax=1000;
hmin=0; hdel=20; hmax=220;
tmin=0; tdel=100; tmax=600;

subplot(1,2,1)
plot(xx/1000,hh/1000,'r','Linewidth',lwd); hold on;
set (gca,'xlim',[Xmin Xmax],'xtick',[Xmin:Xdel:Xmax],'fontsize',font);
xlabel('Range (km)  ','fontsize',font);set(gca,'xcolor','black','fontsize',font);
set (gca,'ylim',[hmin hmax],'ytick',[hmin:hdel:hmax],'fontsize',font);
ylabel('Altitude (km)  ','fontsize',font);set(gca,'xcolor','black','fontsize',font);
grid on
subplot(1,2,2)
plot(tt,hh/1000,'r','Linewidth',lwd); hold on;
set (gca,'xlim',[tmin tmax],'xtick',[tmin:tdel:tmax],'fontsize',font);
xlabel('Time (s)  ','fontsize',font);set(gca,'xcolor','black','fontsize',font);
set (gca,'ylim',[hmin hmax],'ytick',[hmin:hdel:hmax],'fontsize',font);
ylabel('Altitude (km)  ','fontsize',font);set(gca,'xcolor','black','fontsize',font);
grid on

nn=length(tt);

for n=2:nn-1
    ang(n)=atan((hh(n+1)-hh(n-1))/(xx(n+1)-xx(n-1)));
%   ang(n)=90*pi/180;
end
ang(1)=ang(2);ang(nn)=ang(nn-1);

gloss=0;
for n=2:nn-1
    gloss=gloss+g0*sin(ang(n))*(tt(n+1)-tt(n-1))/2;
    for m=1:mm
        if(n==ind(m));gls(m)=gloss;end
    end
end
gls(mm)=gloss;

M(1)=m_LO;
M(2)=m_LO-t(2)*mdotC;
M(3)=m_LO-t(3)*mdotC;
M(4)=m_LO-t(4)*mdotC;
t_ig2=t_SIG;
M(5)=m_tot2-mdot2*(t(5)-t_ig2);
M(6)=m_tot2-mdot2*(t(6)-t_ig2);
M(7)=m_tot2-mdot2*(t(7)-t_ig2);
M(8)=m_tot2-mdot2*(t(8)-t_ig2);
M(9)=m_tot2-mdot2*(t(9)-t_ig2);
MM=interp1(t,M,tt,'pchip');

dloss=0;
for n=2:nn-1
    rho=6.9579e-6;
    if (hh(n) < 84852); [dum dum dum rho]=atmoscoesa(hh(n)); end
    q=0.5*rho*(vv(n)^2);
    dloss=dloss+ ((q*Cd*A)/MM(n))*(tt(n+1)-tt(n-1))/2;
    for m=1:mm
        if(n==ind(m));dls(m)=dloss;end
    end
end
dls(mm)=dloss;

%[t' h' x' v' gls' dls']
fprintf('\n     time(s)   Alt(km)  Range(km) DelV (m/s) gloss (m/s) dloss (m/s)\n')
for m=1:mm
    fprintf('%10.0f %10.1f %10.1f %10.0f %10.0f %10.0f\n',t(m),h(m)/1000,x(m)/1000,veq(m),gls(m),dls(m))
end

fprintf('\n     Total gravity loss  = %5.1f (m/s)',gloss)
fprintf('\n     Total drag loss     = %5.1f (m/s)',dloss)
fprintf('\n     mpayl               = %5.1f (kg)',m_payl)
fprintf('\n Spacecraft speed         = %5.1f (m/s)',DelV-gloss-dloss)
fprintf('\n');

% ------------------------  Output -----------------------------

%  MR2 m_tot2 = 8.0395   102785
%  MR1 m_tot1 m_LP m_LO = 4.0961   402752   507287   500901.4176
% 
%  T/W = 1.1983
%  DelV1 = 4099.6575 m/s
%  DelV2 = 6949.3624 m/s
%  DelV_rot = 288 m/s
%  DelV  = 11337.5668 m/s
%  V orbit = 7734 m/s
%  DelV orbit = 8066.0338 m/s
%  DelV loss (actual?) = 3271.533 m/s
% 
%  Beginning of Gravity Turn
%  Final speed               =  73.0 (m/s)
%  Final angle               =  90.0 (deg)
%  Altitude                  =   0.8 (km)
%  Downrange                 =   0.0 (km)
%  Gravity loss              = 225.5 (m/s)
%  Drag loss                 =   0.5 (m/s)
%  First  Stage mi mf        = 500901 451945(kg)
% 
%  End of First Stage Burn
%  Final speed               = 2712.5 (m/s)
%  Final angle               =  22.9 (deg)
%  Altitude                  =  79.7 (km)
%  Downrange                 = 107.9 (km)
%  Gravity loss              = 1337.4 (m/s)
%  Drag loss                 =  49.1 (m/s)
%  First  Stage mi mf        = 451945 122286(kg)
%  DelV1                     = 4255.6 (m/s)
% 
%  First Stage Impact DR     =   889 (km)
%  First Stage Apogee        =   136 (km)
% 
%  Fairing Jettison
%  Final speed               = 2916.9 (m/s)
%  Final angle               =  16.6 (deg)
%  Altitude                  = 119.4 (km)
%  Downrange                 = 219.2 (km)
%  Gravity loss              = 1472.4 (m/s)
%  Drag loss                 =  49.1 (m/s)
%  Second  Stage mi mf       = 104535 94414(kg)
% 
%  Fairing Impact DR         =  1077 (km)
%  Fairing Apogee            =   154 (km)
% 
%  End of Second Stage Burn
%  Final speed               = 9334.0 (m/s)
%  Final angle               =  -0.0 (deg)
%  Altitude                  = 185.3 (km)
%  Downrange                 = 1887.4 (km)
%  Gravity loss              = 1659.1 (m/s)
%  Drag loss                 =  49.1 (m/s)
%  Second  Stage mi mf       = 92664 12785(kg)
%  DelV2                     = 6577.9 (m/s)
% 
%  Second Stage Impact DR    =   NaN (km)
%  Second Stage Apogee       =   NaN (km)
% 
%  DelV                      = 11122.1 (km/s)
% 
%      Altitude            = 185.3 (km)
%      Spacecraft speed    = 9334.0 (m/s)
%      V orbit             = 7734.0 (m/s)
%      Total gravity loss  = 1659.1 (m/s)
% 
%      Total drag loss     =  49.1 (m/s)
%      mpayl               = 9635.0 (kg)
%      Launchpad Mass      = 507287.0 (kg)
% 
%      time(s)   Alt(km)  Range(km) DelV (m/s) gloss (m/s) dloss (m/s)
%          0        0.0        0.0          0          0          0
%         60        5.5       11.6        413        232         10
%        120       35.0       36.0       1299        680         46
%        150       58.0       59.8       2092        884         49
%        180       90.0       99.2       2318       1069         49
%        240      134.0      178.0       2814       1355         50
%        300      166.0      278.7       3501       1533         50
%        435      205.0      572.9       5202       1715         51
%        540      208.0      946.6       7864       1723         58
%        560      210.0     1000.2       8060       1730         68
% 
%      Total gravity loss  = 1730.4 (m/s)
%      Total drag loss     =  67.6 (m/s)
%      mpayl               = 7885.0 (kg)
%      Spacecraft speed    = 9324.0 (m/s)
