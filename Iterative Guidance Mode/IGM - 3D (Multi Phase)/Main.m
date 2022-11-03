%% IGM-3D (Multi-Phase)

clear;close all;clc
format short g
set(0,'DefaultLineLineWidth',2);

p.GM   = 3.986004418e14;
p.RE   = 6378137;
p.g0   = 9.80665;

p.Thr3a = 7303;
p.Thrc  = 0.0;
p.Thr3b = 7303;

p.m03a = 1332.18;
p.m0c  = 911.98;
p.m03b = p.m0c;

p.mdot3a = 2.4823;
p.mdotc  = 0.0;
p.mdot3b = 2.4823;

p.Isp3a = p.Thr3a/p.mdot3a/p.g0;
p.vex3a = p.Isp3a*p.g0;

p.Isp3b = p.Thr3b/p.mdot3b/p.g0;
p.vex3b = p.Isp3b*p.g0;

p.tb3a = 169.28;
p.tc   = 145.446;
p.tb3b = 160.0;

p.Tsamp = 20e-3;

% Initial Conditions
x0 = 429274.2166;
y0 = 6556331.5265;
z0 = -92502.297615;
vx0 = 4971.958297;
vy0 = 1127.1250285;
vz0 = -398.2442616;

X0 = [x0;y0;z0;vx0;vy0;vz0];

p.h      = 500e3;
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
p.azim0  = 190.73;
p.incf   = 97.414;
p.omegaf = -(180.0-177.46);

rotxx = @(a) [1,0,0;0,cosd(a),sind(a);0,-sind(a),cosd(a)];
rotyy = @(a) [cosd(a),0,-sind(a);0,1,0;sind(a),0,cosd(a)];
rotzz = @(a) [cosd(a),sind(a),0;-sind(a),cosd(a),0;0,0,1];

p.A = rotyy(p.azim0-90);
p.B = rotxx(p.lat0);
p.C = rotzz(-p.omegaf);
p.D = rotyy(-p.incf);

p.G = p.D*p.C*p.B*p.A; % Constant Rotation Matrix

tic;out = sim('IGM_3D_MP');toc;

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