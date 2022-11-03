%% IGM-2D (Single Phase)

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
vx0 = 5575.0;
vy0 = 843.264;

X0 = [x0;y0;vx0;vy0];

p.h      = 900e3;
p.rf     = p.h+p.RE;
p.gammaf = 0.0;
p.vf     = sqrt(p.GM/p.rf);

p.zaiT = 0;
p.etaT = p.rf;

p.zaidotT  = p.vf*cos(p.gammaf);
p.etadotT  = p.vf*sin(p.gammaf);

tic;out = sim('IGM_2D_SP');toc;

outputs = out.simout;
mag = @(v) (sqrt(sum(v.^2,2)));

figure(1);hold on;grid on;box on
plot(outputs.time,(mag(outputs.data(:,1:2))-p.RE)/1e3,'r.-');
xlabel('t [sec]');ylabel('h [km]');title('Altitude');

figure(2);hold on;grid on;box on
plot(outputs.time,mag(outputs.data(:,3:4)),'r.-');
xlabel('t [sec]');ylabel('v [m/s]');title('Velocity');

figure(3);hold on;grid on;box on
plot(outputs.time(20:end),outputs.data(20:end,5),'r.-');
xlabel('t [sec]');ylabel('\theta [deg]');title('Pitch Angle');