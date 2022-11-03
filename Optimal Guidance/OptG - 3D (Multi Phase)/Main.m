%% Optimal Guidance (3D Multi Phase - Spherical Earth)

clear;close all;clc
format long g;warning off
set(0,'DefaultLineLineWidth',2);

p.GM   = 3.986004418e14;
p.RE   = 6371e03;
p.g0   = 9.80665;

p.Thr3  = 225e3; % Thrust (Stage 3)
p.m03   = 12e3;  % Propellant Mass (Stage 3)
p.mdot3 = 75.0;  % Propellant Consumption Rate (Stage 3)

p.ThrC  = 0.0;   % Thrust (Coasting Phase)
p.m0C   = 3700;  % Propellant Mass (Coasting Phase)
p.mdotC = 0.0;   % Propellant Consumption Rate (Coasting Phase)

p.Thr4  = 15e3;  % Thrust (Stage 4)
p.m04   = 3700;  % Propellant Mass (Stage 4)
p.mdot4 = 5.0;   % Propellant Consumption Rate (Stage 4)

% Initial Conditions
tb3 = 100; % Burnout Time Stage 3
tc  = 350; % Coasting Time
tb4 = 100; % Burnout Time Stage 4

T3  = tb3;
TC  = tc + T3;
T4  = 750;

% Initial Conditions
p.x03  = 247541.9;
p.y03  = 6531426;
p.z03  = -74628.75;
p.vx03 = 3429.384;
p.vy03 = 1445.38;
p.vz03 = -401.6901;

p.h      = 900e3;
p.rf     = p.h+p.RE;
p.gammaf = 0.0;
p.lat0   = 13.0;
p.azim0  = 192-360;
p.incf   = 99.03;
p.omegaf = -176.89;

%% Simulation

% Initial Guesses for Solution
n = 100; % Number of Grid Points
x_guess = linspace(0,1,n);
        
y_guess = [ p.x03;    %x          (Stage 03)
            p.y03;    %y          (Stage 03)
            p.z03;    %z          (Stage 03)
            p.vx03;   %vx         (Stage 03)
            p.vy03;   %vy         (Stage 03)
            p.vz03;   %vz         (Stage 03)
           -1e-8;    %lambda_x   (Stage 03)
           -1e-8;    %lambda_y   (Stage 03)
           -1e-8;    %lambda_z   (Stage 03)
           -0.82;    %lambda_vx  (Stage 03)
           -0.18;    %lambda_vy  (Stage 03)
            0.0074   %lambda_vz  (Stage 03)
            723246.3;     %x          (Coasting Stage)
            6675327;     %y          (Coasting Stage)
            -114196.7;     %z          (Coasting Stage)
            5709.542;    %vx         (Coasting Stage)
            1282.741;    %vy         (Coasting Stage)
            -388.9282;    %vz         (Coasting Stage)
           -0.0001;  %lambda_x   (Coasting Stage)
           -0.001;   %lambda_y   (Coasting Stage)
           -1e-06;   %lambda_z   (Coasting Stage)
           -0.89;    %lambda_vx  (Coasting Stage)
            0.3;     %lambda_vy  (Coasting Stage)
            0.001    %lambda_vz  (Coasting Stage)
            2631771;    %x          (Stage 04)
            6618294;    %y          (Stage 04)
            -238704.1;    %z          (Stage 04)
            5343.232;   %vx         (Stage 04)
            -1586.591;   %vy         (Stage 04)
            -315.9531;   %vz         (Stage 04)
           -0.0001;  %lambda_x   (Stage 04)
           -0.001;   %lambda_y   (Stage 04)
           -1e-06;   %lambda_z   (Stage 04)
           -0.89;    %lambda_vx  (Stage 04)
            0.3;     %lambda_vy  (Stage 04)
            0.001
            ];  %lambda_vz  (Stage 04)        

Tf_guess = [TC;T4]; % Guess for Time of Flight
solinit = bvpinit(x_guess,y_guess);
solinit.parameters = Tf_guess;

% Solution Options
tol = 1E-10; % Tolerance
Options = bvpset('RelTol',tol,'AbsTol',tol,'Nmax',2000);

% Solution
sol = bvp4c(@odes_3d_mp,@bcs_3d_mp,solinit,Options,T3,p);

Time_of_flight = sol.parameters(2);
Coasting_Time  = sol.parameters(1)-T3;

[res] = bcs_3d_mp(sol.y(:,1),sol.y(:,end),sol.parameters,T3,p);

%%

t1 = sol.x*T3;
t2 = sol.x*(sol.parameters(1)-T3)+T3;
t3 = sol.x*(sol.parameters(2)-sol.parameters(1))+sol.parameters(1);

t      = [t1,t2,t3];
x      = [sol.y(1,:)  sol.y(13,:)  sol.y(25,:)];
y      = [sol.y(2,:)  sol.y(14,:)  sol.y(26,:)];
z      = [sol.y(3,:)  sol.y(15,:)  sol.y(27,:)];
vx     = [sol.y(4,:)  sol.y(16,:)  sol.y(28,:)];
vy     = [sol.y(5,:)  sol.y(17,:)  sol.y(29,:)];
vz     = [sol.y(6,:)  sol.y(18,:)  sol.y(30,:)];
lam_x  = [sol.y(7,:)  sol.y(19,:)  sol.y(31,:)];
lam_y  = [sol.y(8,:)  sol.y(20,:)  sol.y(32,:)];
lam_z  = [sol.y(9,:)  sol.y(21,:)  sol.y(33,:)];
lam_vx = [sol.y(10,:) sol.y(22,:)  sol.y(34,:)];
lam_vy = [sol.y(11,:) sol.y(23,:)  sol.y(35,:)];
lam_vz = [sol.y(12,:) sol.y(24,:)  sol.y(36,:)];

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