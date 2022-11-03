%% IGM-3D (Single Phase)

clear;close all;clc
format short g
set(0,'DefaultLineLineWidth',2);

p.GM   = 3.986004418e14;
p.RE   = 6371e03;
p.g0   = 9.80665;
p.Thr  = 14.70995e3;
p.m0   = 3630.0;
p.mdot = 4.765;
p.Isp  = p.Thr/p.mdot/p.g0;
p.vex  = p.Isp*p.g0;

p.Tsamp = 20e-3;

% Initial Conditions
x0  = 0.0;
y0  = 7085e03;
z0  = 0.0;
vx0 = 5575.0;
vy0 = 843.264;
vz0 = 0.0;

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

p.lat0   = 13.0;
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

tic;out = sim('IGM_3D_SP');toc;

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