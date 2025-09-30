clear all, clc

%Total 25 Parameters would need to be defined
%
% Mphiphi Mdelphi Mdeldel Mphidel
% 
% Kophiphi Kophidel Kodelphi Kodeldel K2phiphi K2phidel K2delphi K2deldel
% 
% Cophiphi Cophidel Codelphi Codeldel 
% 
% Tphi Tdel
% 
% v gravacc

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
K2phidel = (1/w)*(ST - mT*zT)*cos(lambda);
K2delphi = 0;
K2deldel = (1/w)*(SA + SF*sin(lambda))*cos(lambda);


%Damping Coefficeints
Cophiphi = 0;
Cophidel = mu*ST + SF*cos(lambda) + (ITXZ/w)*cos(lambda) - mu*mT*zT;
Codelphi = -(mu*ST + SF*cos(lambda));
Codeldel = (IALZ/w)*cos(lambda) + mu*(SA + (ITZZ/w)*cos(lambda));


M = [
    Mphiphi Mphidel;
    Mdelphi Mdeldel];

Ko = [
    Kophiphi Kophidel;
    Kodelphi Kodeldel];

K2 = [
     K2phiphi K2delphi;
     K2delphi K2deldel
     ];


Co = [
    Cophiphi Cophidel;
    Codelphi Codeldel
    ];

%Test Parameters
Tphi = 0;
Tdel = 0;


v = 5; %mpersec
gravacc = 9.81;

A = [ 
        Mphiphi Mphidel;
        Mdelphi Mdeldel
    ];


B = [
        v*Cophiphi v*Cophidel;
        v*Codelphi v*Codeldel
        ];

C = [
        gravacc*Kophiphi + v^2*K2phiphi gravacc*Kophidel + v^2*K2phidel;
        gravacc*Kodelphi + v^2*K2delphi gravacc*Kodeldel + v^2*K2deldel
        ];

D = [
        Tphi;
        Tdel
    ];

%Making State Space Matrices

Ass = [
        zeros(2,2) eye(2);
        -A\C -A\B
        ];

Bss = [
        zeros(2,2);
        A\eye(2)
        ];


Css = [
        eye(4)
        ];


Dss = zeros(4, 1);

%%% Variable Definition Completed

%Since we will be working with only one input
Bss = Bss(:,2);

% create Continuous state-space object (requires Control System Toolbox)
sys = ss(Ass, Bss, Css, Dss);
x0 = [0.15 ;0 ;1.5 ;0];

%simulate the state space system
t = 0:0.01:10;
u = zeros(length(t),1);   
[y,t,x] = lsim(sys, u, t,x0);

% Creating Discrete State Space
Ts = 0.01;            

% discrete plant
sysc = ss(Ass, Bss, Css, Dss);  
sysd = c2d(sysc, Ts, 'zoh');
Ad = sysd.A; Bd = sysd.B; Cd = sysd.C; Dd = sysd.D;

% Making the Discrete LQR Controller
% LQR design (discrete)
Q = diag([0.1, 0.1, 1.5, 0.1]);   % tune these %Roll Angle, Steer Angle, Roll Rate, Steer Rate
R = 1;               % One actuation: Steer torque
Kd = dlqr(Ad, Bd, Q, R);        % discrete state-feedback Kd (2x4) -> u[k] = -Kd*x[k]

% Discrete Kalman (dlqe): model process and measurement noise covariances
n = size(Ad,1); p = size(Cd,1);
Gd = eye(n);                    % assume process noise on each state (tweak if you have better model)
Qw = 1e-4 * eye(n);             % process noise covariance (tune)
Rv = 1e-2 * eye(p);             % measurement noise covariance (tune)
[Ld, ~, ~] = dlqe(Ad, Gd, Cd, Qw, Rv);  % Ld is the discrete estimator gain (n x p)

%% Discete Kalman State Space
% x[k+1] = Ad*x[k] + Bd*u[k] + Ld*(y[k] - Cd*x[k])
% x[k+1] = (Ad - Ld*Cd)x[k] + Bd*u[k] + Ld*y[k]

Aobs = Ad - Ld*Cd;
Bobs = [Bd, Ld];
Cobs = eye(4);
Dobs = zeros(4, size(Bobs,2));
xhat0 = [0.15 ;0 ;1.5 ;0];



%Q = diag([0.1, 0.1, 1.5, 0.1]); 