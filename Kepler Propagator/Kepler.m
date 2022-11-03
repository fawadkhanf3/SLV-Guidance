clear;close all;clc
format short g

p.mu = 3.986004418e14;
p.RE = 6371e3;

r1 = [192.9094e3;
      0.0e3;
      140.5462e3+p.RE];
     
v1 = [-2.9593e3;
      -0.0083e3;
      -0.9686e3];

t0 = 0;
tf = 36.1028;

dotT = @(t,x,p) [x(4);
                 x(5);
                 x(6);
                 -p.mu/(sqrt(x(1)^2+x(2)^2+x(3)^2))^3*x(1);
                 -p.mu/(sqrt(x(1)^2+x(2)^2+x(3)^2))^3*x(2);
                 -p.mu/(sqrt(x(1)^2+x(2)^2+x(3)^2))^3*x(3)];
                 
[~,yf] = ode45(@(t,x) dotT(t,x,p),[t0,tf],[r1;v1]);

r2a = yf(end,1:3);
v2a = yf(end,4:6);

[r2b,v2b] = keplers_propagator([t0,tf],r1,v1,p);

%([r2a(:) r2b(:)]-[0,0;0,0;p.RE,p.RE])/1e3
r_error = r2a(:)-r2b(:);

disp('Position Error (m)');
disp(r_error);

% ([v2a(:) v2b(:)])
v_error = v2a(:)-v2b(:);

disp('Velocity Error (m/s)');
disp(v_error);

