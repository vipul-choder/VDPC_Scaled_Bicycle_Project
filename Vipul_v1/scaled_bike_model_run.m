clear all, clc

%Total 25 Parameters would need to be defined

global Mphiphi Mdelphi Mdeldel Mphidel

global Kophiphi Kophidel Kodelphi Kodeldel K2phiphi K2phidel K2delphi K2deldel

global Cophiphi Cophidel Codelphi Codeldel 

global v gravacc

% Rear Wheel Parameters Values;
mR = 2; %Kgs
iRXX = 0.0603 ; %kgm^2
iRYY = 0.12; %kgm^2
rR = 0.3; %m

%Rear Frame Parameter Values;
mB = 85; %kgs
xB = 0.3; %m
zB = -0.9; %m
iBXX = 9.2;
iBXZ = 2.4;
iBZZ = 2.8;
iBYY = 11;

%Front Frame Parameter Values;
mH = 4.5; %kgs
xH = 0.9; %m
zH = -0.7; %m
iHXX = 0.0589;
iHXZ = -0.00756;
iHZZ = 0.00708;
iHYY = 0.06;

% Front Wheel Parameters Values;
mF = 3; %Kgs
iFXX = 0.1405 ; %kgm^2
iFYY = 0.28; %kgm^2
rF = 0.35; %m

%Vehicle Parameter values
w = 1.02; %Wheelbase, m
lambda = pi/10;
c = 0.08; %Pneumatic Trail


%Getting coefficient matrices
mT = mR + mB + mH + mF;
xT = (1/mT)*(xB*mB + xH*mH + w*mF);
zT = (1/mT)*(-rR*mR + zB*mB + zH*mH -rF*mF);

ITXX = iRXX + iBXX + iHXX + iFXX + mR*rR^2 + mB*zB^2 + mH*zH^2 + mF*rF^2;

ITXZ = iBXZ + iHXZ - mB*xB*zB - mH*xH*zH + mF*w*rF;

iRZZ = iRXX;
iFZZ = iFXX;

ITZZ = iRZZ + iBZZ + iHZZ + iFZZ + mB*xB^2 + mH*xH^2 + mF*w^2;


%Front Assembly 

mA = mH + mF;
xA = (1/mA)*(xH*mH + w*mF);
zA = (1/mA)*(zH*mH - rF*mF);

IAXX = iHXX + iFXX + mH*(zH - zA)^2 + mF*(rF + zA)^2;
IAXZ = iHXZ - mH*(xH - xA)*(zH - zA) + mF*(w - xA)*(rF + zA);
IAZZ = iHZZ + iFZZ + mH*(xH - xA)^2 + mF*(w - xA)^2;

uA = (xA - w - c )*cos(lambda) - zA*sin(lambda);


%Moment of Inertias about the Steering Axis
IALL = mA*uA^2 + IAXX*(sin(lambda))^2 + 2*IAXZ*sin(lambda)*cos(lambda) + IAZZ*(cos(lambda)^2);

IALX = -mA*uA*zA + IAXX*sin(lambda) + IAXZ*cos(lambda);

IALZ = mA*uA*xA + IAXZ*sin(lambda) + IAZZ*cos(lambda);

%Definition of mu

mu = (c/w)*cos(lambda);

%Definition of Gyrostatic Coefficients
SR = iRYY/rR;
SF = iFYY/rF;
ST = SF + SR;

SA = mA*uA + mu*mT*xT;


% Definition of coefficients

%Inertia coefficients
Mphiphi = ITXX;
Mphidel = IALX + mu*ITXZ;
Mdelphi = Mphidel;
Mdeldel = IALL + 2*mu*IALZ + mu*mu*ITZZ;

%Stiffness Coefficients
Kophiphi = mT*zT;
Kophidel = -SA;
Kodelphi = Kophidel;
Kodeldel = -SA*sin(lambda);

K2phiphi = 0;
K2phidel = (1/w)*(ST - mT*zT*cos(lambda));
K2delphi = 0;
K2deldel = (1/w)*(SA + SF*sin(lambda)*cos(lambda));


%Damping Coefficeints
Cophiphi = 0;
Cophidel = mu*ST + SF*cos(lambda) + (ITXZ/w)*cos(lambda) - mu*mT*zT;
Codelphi = -(mu*ST + SF*cos(lambda));
Codeldel = (IALZ/w)*cos(lambda) + mu*(SA + (IAZZ/w)*cos(lambda));

%%% Variable Definition Completed

%Test param values
v = 5; %mpersec
gravacc = 9.81;
x0=[0 0.5 0 0];

tstop=2.5;             % stop-time
dt=0.0002;            % timestep
time=0:dt:tstop;    % time vector (fixed time step)



[timeout, xout]=ode15s(@scaled_bicycle_model,time,x0);
%% 

% Identify variables from simulation results
%------------------------------------------------------
roll_rate=xout(:,2);
steer_rate=xout(:,4);

roll_angle=xout(:,1);
steer_angle = xout(:,3);

figure;
plot(timeout, steer_rate, 'Color','r',LineWidth=2, DisplayName='SteerRate');
hold on;
yyaxis right;
plot(timeout, roll_rate, 'Color','b',LineWidth=2, DisplayName='RollRate');
legend();
