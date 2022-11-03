%% Optimal Guidance (3D Single Phase - Spherical Earth)

clear;close all;clc
format long g;warning off
set(0,'DefaultLineLineWidth',2);

p.GM   = 3.986004418e14;
p.RE   = 6371e03;
p.g0   = 9.80665;
p.Thr  = 14.70995e3;
p.m0   = 3630.0;
p.mdot = 4.765;

p.x0  = 0.0;
p.y0  = 7085e03;
p.z0  = 0.0;
p.vx0 = 5575.0;
p.vy0 = 843.264;
p.vz0 = 0.0;

p.h      = 900e3;
p.rf     = p.h+p.RE;
p.gammaf = 0.0;
p.lat0   = 13.0;
p.azim0  = 192-360;
p.incf   = 99.03;
p.omegaf = -176.89;

%%

% Initial Guesses for Solution
n = 100; % Number of Grid Points
x_guess = linspace(0,1,n);

y_guess = [ ones(1,n)*1117896;    % x
            ones(1,n)*6854812;    % y
            ones(1,n)*-163321.3;    % z
            ones(1,n)*4320.543;   % vx
            ones(1,n)*629.1032;   % vy
            ones(1,n)*-366.6159;   % vz
            ones(1,n)*-0.0001;  % lambda_x
           -ones(1,n)*0.0001;  % lambda_y
            ones(1,n)*-0.001;  % lambda_z
           -ones(1,n)*0.91;   % lambda_vx
            ones(1,n)*0.036;   % lambda_vy
            ones(1,n)*-0.16]; % lambda_vz
        
Tf_guess = 380; % Guess for Time of Flight

solinit.x = x_guess;
solinit.y = y_guess;  
solinit.parameters = Tf_guess;

% Solution Options
tol = 1E-10; % Tolerance
Options = bvpset('RelTol',tol,'AbsTol',tol,'Nmax',2000);

% Solution
sol = bvp4c(@odes_3d_sp,@bcs_3d_sp,solinit,Options,p);

Time_of_flight = sol.parameters;
res = bcs_3d_sp(sol.y(:,1),sol.y(:,end),sol.parameters,p);

%%
t = sol.x*Time_of_flight;
x = sol.y(1,:);
y = sol.y(2,:);
z = sol.y(3,:);
vx = sol.y(4,:);
vy = sol.y(5,:);
vz = sol.y(6,:);
lam_x  = sol.y(7,:);
lam_y  = sol.y(8,:);
lam_z  = sol.y(9,:);
lam_vx = sol.y(10,:);
lam_vy = sol.y(11,:);
lam_vz = sol.y(12,:);

%%
mag = @(v) (sqrt(sum(v.^2,2)));

theta = atan2d(-lam_vy,-lam_vx);
psi   = asind(-lam_vz(:)./mag([lam_vx(:),lam_vy(:),lam_vz(:)]));

figure(1);hold on;grid on;box on
plot(t,(mag([x(:),y(:),z(:)])-p.RE)/1e3,'r.-');
xlabel('t [sec]');ylabel('h [km]');title('Altitude');

figure(2);hold on;grid on;box on
plot(t,mag([vx(:),vy(:),vz(:)]),'r.-');
xlabel('t [sec]');ylabel('v [m/s]');title('Velocity');

figure(3);hold on;grid on;box on
plot(t,theta,'r.-');
xlabel('t [sec]');ylabel('\theta [deg]');title('Pitch Angle');

figure(4);hold on;grid on;box on
plot(t,psi,'r.-');
xlabel('t [sec]');ylabel('\psi [deg]');title('Yaw Angle');