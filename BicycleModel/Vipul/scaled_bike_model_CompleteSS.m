clear all, clc

%Total 25 Parameters would need to be defined

global Mphiphi Mdelphi Mdeldel Mphidel

global Kophiphi Kophidel Kodelphi Kodeldel K2phiphi K2phidel K2delphi K2deldel

global Cophiphi Cophidel Codelphi Codeldel 

global Tphi Tdel

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


v = 20; %mpersec
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


%% Lets try deriving the Lateral positions derivative, yp, psi and yq

% psi_dot = (((v*delta) + c*(delta_dot))/w)*cos(lambda);
% psi_ddot = [0; diff(psi_dot)];
% 
% 
% yp_dot = v*psi_dot;
% yp_ddot = [0;diff(yp_dot)];
% 
% figure;
% plot(t,yp_dot, Linewidth = 2, LineStyle="-", DisplayName = "yp_dot");
% hold on;
% plot(t,yp_ddot, Linewidth = 2, LineStyle=":", DisplayName = "yp_ddot");
% 
% legend();

%% Lets try to get some forces

% TBphi = -mT.*yp_ddot.*zT + ITXX.*phi_ddot + ITXZ.*psi_ddot + IALZ.*delta_ddot + psi_dot.*v*ST + delta_dot.*v*SF*cos(lambda) + gravacc*mT*zT.*phi - gravacc*SA*delta;
% 
% figure;
% plot(t,TBphi, Linewidth = 2, LineStyle="-", DisplayName = "Lean Torque [Nm]");
% hold on;
% yyaxis right;
% plot(t, phi, Linewidth = 2, LineStyle="-", DisplayName = "Roll Angle [rad]");
% legend();
% 
% 
% FFy = (1/w)*(mT.*yp_ddot.*xT + ITXZ.*phi_ddot + ITZZ.*psi_ddot + IALZ.*delta_ddot - phi_dot.*v*ST - delta_dot.*v*SF*sin(lambda));
% 
% figure;
% plot(t,FFy, Linewidth = 2, LineStyle="-", DisplayName = "Lateral Force [N]");
% hold on;
% yyaxis right;
% plot(t, phi, Linewidth = 2, LineStyle="-", DisplayName = "Roll Angle [rad]");
% legend();

%% Lets try to get a Velocity sweep for Eigenvalues

vRange = linspace(0,10,100);

eigVals = [];

for vi = vRange
    A_ins = [ 
            Mphiphi Mphidel;
            Mdelphi Mdeldel
    ];


    B_ins = [
            vi*Cophiphi vi*Cophidel;
            vi*Codelphi vi*Codeldel
        ];

    C_ins = [
            gravacc*Kophiphi + vi^2*K2phiphi gravacc*Kophidel + vi^2*K2phidel;
            gravacc*Kodelphi + vi^2*K2delphi gravacc*Kodeldel + vi^2*K2deldel
        ];

    A_eig = [
            zeros(2,2) eye(2);
            -A_ins\C_ins -A_ins\B_ins
            ];
    
    lambda_ins = eig(A_eig);

    eigVals = [eigVals; lambda_ins.'];

end

% Plot eigenvalues vs speed
figure; hold on; grid on;
for mode = 1:size(eigVals,2)
    plot(vRange, real(eigVals(:,mode)),'-','LineWidth',1.5);
end
xlabel('Forward speed v [m/s]');
ylabel('Real part of eigenvalue');
title('Bicycle Eigenvalue Speed Sweep');
legend('Mode 1','Mode 2','Mode 3','Mode 4');