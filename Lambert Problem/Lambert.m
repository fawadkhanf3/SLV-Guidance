%% Lambert (Required Velocity)

clear;close all;clc
format long g

p.G  = 6.67e-11; % Universal Gravitional Constant
p.M  = 5.972e24; % Mass of Earth (kg)
p.RE = 6371e3; % Radius of Earth (m)

r1 = [125044;  
      6540831];

r2 = [2333733;
      5928180]; 

delta_t = 960;

[V,gamma] = lambert_velocity(delta_t,r1,r2,p);

