%% IGM-3D (3Phase)

clear;close all;clc
format short g
set(0,'DefaultLineLineWidth',2);

p.GM   = 3.986004418e14;
p.RE   = 6378137;
p.g0   = 9.80665;

p.Thr3a = 225000.0;
p.Thrc  = 0.0;
p.Thr3b = 15000.0;

p.m03a = 12000.0;
p.m0c  = 3700.0;
p.m03b = p.m0c;

p.mdot3a = 75.0;
p.mdotc  = 0.0;
p.mdot3b = 5.0;

p.Isp3a = p.Thr3a/p.mdot3a/p.g0;
p.vex3a = p.Isp3a*p.g0;

p.Isp3b = p.Thr3b/p.mdot3b/p.g0;
p.vex3b = p.Isp3b*p.g0;

p.tb3a = 100.0;
p.tc   = 393.59;
p.tb3b = 320.0;

p.Tsamp = 20e-3;

% Initial Conditions
x0 = 247541.9;
y0 = 6531426;
z0 = -74628.75;
vx0 = 3429.384;
vy0 = 1445.38;
vz0 = -401.6901;

X0 = [x0;y0;z0;vx0;vy0;vz0];

p.h      = 900e3;
p.rf     = p.h+p.RE;
p.gammaf = 0.0;
p.vf     = sqrt(p.GM/p.rf);

p.zaiT = 0;
p.etaT = p.rf;
p.psiT = 0;

p.zaidotT  = p.vf*cos(p.gammaf);
p.etadotT  = p.vf*sin(p.gammaf);
p.zetadotT = 0.0;

p.lat0   = 25.347755;
p.azim0  = 192-360;
p.incf   = 99.03;
p.omegaf = -(180.0-176.89);

rotxx = @(a) [1,0,0;0,cosd(a),sind(a);0,-sind(a),cosd(a)];
rotyy = @(a) [cosd(a),0,-sind(a);0,1,0;sind(a),0,cosd(a)];
rotzz = @(a) [cosd(a),sind(a),0;-sind(a),cosd(a),0;0,0,1];

p.A = rotyy(p.azim0-90);
p.B = rotxx(p.lat0);
p.C = rotzz(-p.omegaf);
p.D = rotyy(-p.incf);

p.G = p.D*p.C*p.B*p.A; % Constant Rotation Matrix

tt1 = [0 p.tb3a];
tt2 = [p.tb3a p.tb3a+p.tc];
tt3 = [p.tb3a+p.tc p.tb3a+p.tc+p.tb3b];

kk1 = [1 1.8];
kk2 = [2.25 1];
kk3 = [1 2.25];

p.K1 = polyfit(tt1,kk1,1);
p.K2 = polyfit(tt2,kk2,1);
p.K3 = polyfit(tt3,kk3,1);

tic;out = sim('IGM_3D_MP_LC');toc;

outputs = out.simout;
mag = @(v) (sqrt(sum(v.^2,2)));

figure(1);hold on;grid on;box on
plot(outputs.time,(mag(outputs.data(:,1:3))-p.RE)/1e3,'r.-');
xlabel('t [sec]');ylabel('h [km]');title('Altitude');

figure(2);hold on;grid on;box on
plot(outputs.time,mag(outputs.data(:,4:6)),'r.-');
xlabel('t [sec]');ylabel('v [m/s]');title('Velocity');

figure(3);hold on;grid on;box on
plot(outputs.time(20:end),outputs.data(20:end,7),'r.-');
xlabel('t [sec]');ylabel('\theta [deg]');title('Pitch Angle');

figure(4);hold on;grid on;box on
plot(outputs.time(20:end),outputs.data(20:end,8),'r.-');
xlabel('t [sec]');ylabel('\psi [deg]');title('Yaw Angle');