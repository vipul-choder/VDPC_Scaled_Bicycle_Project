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


v = 6; %mpersec
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


Dss = zeros(4, 2);

%%% Variable Definition Completed

% create state-space object (requires Control System Toolbox)
sys = ss(Ass, Bss, Css, Dss);
x0 = [0 ;0 ;0.5 ;0];

%simulate the state space system
t = 0:0.01:10;
u = zeros(length(t),2);   
[y,t,x] = lsim(sys, u, t,x0);
% plot(t, y);
% legend();

phi = y(:,1);
delta = y(:,2);
phi_dot = y(:,3);
delta_dot = y(:,4);

delta_ddot = [0; diff(delta_dot)];
phi_ddot = [0; diff(phi_dot)];

figure;
plot(t,delta_dot, Linewidth = 2, LineStyle="-", DisplayName = "DeltaRate");
hold on;
plot(t,phi_dot, Linewidth = 2, LineStyle=":", DisplayName = "RollRate");
legend();

figure;
plot(t,delta, Linewidth = 2, LineStyle="-", DisplayName = "DeltaAngle");
hold on;
plot(t,phi, Linewidth = 2, LineStyle=":", DisplayName = "RollAngle");
legend();


%% Lets try to make a controller for this which would try to get roll to 0 in 1.5 secs
% Will try to use a LQR Controller

%
% if x is your state = [roll, steer, roll rate, steer rate]
% LQR Controller = integral (x'Q x + u'R u)dt
% u = -K * x 
% Q: 4x4 Penalty matrix on the state vector
% R: Penalty scalar coefficient on the use of control input

%Definiing Q - diagonal matrix
Q = diag([100, 10, 1, 1]);

%Most penalty 100  to Roll Angle difference
%Second Most penalty 10 to Steer Angle difference

%Defining R as a scalar
R = 1;

%Getting the Control Variable, K
K = lqr(Ass, Bss, Q, R);  % 1x4 gain vector

%Using the K for my feedback law
%Creating Closed Loop State Matrices
Acl = Ass - Bss*K;
Bcl = [];
Ccl = eye(4);
Dcl = [];

% create state-space object (requires Control System Toolbox)
syscl = ss(Acl, Bcl, Ccl, Dcl);
x0 = [0 ;0 ;0.5 ;0];

t = 0:0.01:10;
u = zeros(length(t),0);  
[y_cl,t,x] = lsim(syscl, u, t,x0);

phi_cl = y_cl(:,1);
delta_cl = y_cl(:,2);
phi_dot_cl = y_cl(:,3);
delta_dot_cl = y_cl(:,4);

figure(1);
plot(t, phi_cl,Linewidth = 2, LineStyle="-", DisplayName="Roll Controlled"); 
hold on;
grid on;
plot(t, phi, Linewidth = 2, LineStyle="--", DisplayName="Roll Uncontrolled");
xlabel('Time [s]'); ylabel('Roll angle [rad]');
title('Bicycle roll under LQR control');
legend();

%Trying to get the Input which was needed to get this control

u = zeros(length(t),2);  % preallocate

for i = 1:length(t)
    u(i,:) = -K * x(i,:)';    % torque at each time step
end

figure(2);
plot(t, u(:,1), 'LineWidth', 2, DisplayName = "Roll Torque");
hold on;
plot(t, u(:,2), 'LineWidth', 2, DisplayName = "Steer Torque");
xlabel('Time [s]');
ylabel('Torque [Nm]');
title('Control input applied by LQR');
grid on;
legend();