%% Optimal Guidance (2D Single Phase - Spherical Earth)

clear;close all;clc
format long g;warning off
set(0,'DefaultLineLineWidth',2);

% Constants
p.G  = 6.67e-11;
p.M  = 5.972e24;
p.GM = p.G*p.M;
p.RE = 6371e3;

% SLV Specifications
p.Thr  = 14.709975e3;
p.m0   = 3630;
p.mdot = 4.765;

% Initial Conditions
p.x0 = 0;
p.y0 = 7085000;
p.vx0 = 5575.0;
p.vy0 = 843.264;

% Terminal Conditions
p.h      = 900e3;
p.rf     = p.h+p.RE;
p.gammaf = 0.0;

% Sim
n = 120;
x_guess = linspace(0,1,n);

y_guess = [ones(1,n)*p.x0;    %x
           ones(1,n)*p.y0;    %y
           ones(1,n)*p.vx0;   %vx
           ones(1,n)*p.vy0;   %vy
           ones(1,n)*4e-4;  %lambda_x
           -ones(1,n)*6e-4; %lambda_y
           -ones(1,n)*1.0;  %lambda_vx
           ones(1,n)*0.2];  %lambda_vy

Tf_guess = 700;

solinit.x = x_guess;
solinit.y = y_guess;
solinit.parameters = Tf_guess;

tol = 1E-10;
Options = bvpset('RelTol',tol,'AbsTol',ones(1,8)*tol,'Nmax',2000);

sol = bvp4c(@odes_2d_sp,@bcs_2d_sp,solinit,Options,p);

Time_of_flight = sol.parameters;
res = bcs_2d_sp(sol.y(:,1),sol.y(:,end),sol.parameters,p);

t = sol.x*Time_of_flight;
x      = sol.y(1,:);
y      = sol.y(2,:);
vx     = sol.y(3,:);
vy     = sol.y(4,:);
lam_x  = sol.y(5,:);
lam_y  = sol.y(6,:);
lam_vx = sol.y(7,:);
lam_vy = sol.y(8,:);

%%
figure(1);hold on;grid on;box on
theta = atan2d(-lam_vy,-lam_vx);
plot(t,theta,'r');
xlabel('t [sec]');
ylabel('\theta [deg]');
title('Pitch Angle');

mag = @(v) (sqrt(sum(v.^2,2)));

figure(2);hold on;grid on;box on
plot(t,(mag([x(:),y(:)])-p.RE)/1e3,'r.-');
xlabel('t [sec]');ylabel('h [km]');title('Altitude');

figure(3);hold on;grid on;box on
plot(t,mag([vx(:),vy(:)]),'r.-');
xlabel('t [sec]');ylabel('v [m/s]');title('Velocity');
