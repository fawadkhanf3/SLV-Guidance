%*************************************************************************
%
%   3-DOF Trajectory simulation program. Computes the trajectory of a
%   2 stage launch vehicle in a cartesian 3 DOF FOR. This script defines     
%   the vehicle parameters and initial position (lat, long) and 
%   makes the call to integrate numerically.


%   The equations of motion are defined within the EoM.m file, and
%   transformation matrices are defined within this too. Vehicle force
%   and mass computation is done in the ForceMass.m file, in which the
%   atmospheric properties are also computed through the Atmosphere.m
%   function.
    
%           V1.0 01/02/2019 - Euler integrator for verification of correct
%           operation and computation.
%
%           21/04/2019 - Redefinition of angles - defined the thrust vector
%           locally with pitch and yaw, which can be transformed inertially
%           to solve the equations of motion
%
%           13/08/2020 - Addition of UPFG routine based on space shuttle 
%           GNC documents. Guidance converges during prethrust cycling,
%           need to now enter main program and issue steering commands.
%
%
%**************************************************************************

clear
clear EqOfMotion theta psi
close all

global M0 mdotProp mdotTot ue S Re tb tb2 M02 ue2 mdotTot2 gamma...
mdotProp2 mu dt

%   EARTH/ATMOSPHERIC PARAMETERS:

Re = 6378E3;
mu = 5.972E24*0.06673E-9;
gamma = 1.4;


%========================VEHICLE VARIABLES:===============================
%STAGE 2 VARIABLES:
Mprop2 = 111.5E3;                                                           %Stage 2 propellant mass [kg]
MS2 = 4.5E3;                                                                %Stage 2 structural mass [kg]
MP2 = 0;                                                                    %Stage 2 payload mass (=M03)
M02 = Mprop2 + MS2 + MP2;                                                   %Stage 2 total mass [kg]
PRatio2 = MP2/M02;                                                          %Stage 2 payload ratio
Sig2 = MS2/Mprop2;                                                          %Stage 2 propellant tankage efficiency
mdotProp2 = Mprop2/397;                                                     %Total Propellant mass flow rate expended by ONE engine [kg/s]
Isp2 = 348;                                                                 %MVAC ISP [s]
ue2 = Isp2*9.81;                                                            %Stage 2 nozzle exit velocity [m/s]
mdotTot2 = (mdotProp2)*1;                                                   %Stage 2 total mass flow rate [kg/s]

%STAGE 1 VARIABLES:
MProp = 418.7E3;                                                            %Stage 1 propellant Mass [kg]
MS1 = 27.2E3;                                                               %Stage 1 structural mass [kg]
MP1 = M02;                                                                  %Stage 1 payload mass (=M02)
M0 = MP1 + MS1 + MProp;                                                     %Stage 1 total mass [kg]
PRatio1 = MP1/M0;                                                           %Stage 1 payload ratio
Sig1 = MS1/MProp;                                                           %Stage 1 propellant tankage efficiency
mdotProp = MProp/161;
mdotTot = (mdotProp);                                                       %Total mass loss of vehicle each second [kg/s]
Isp = 283;
ue = Isp*9.81;                                                              %Propellant exit velocity [m/s]
S = pi * (3.66)^2/4;                                                        %Frontal rocket area for a rocket diameter D = 1.7 m [m^2]
Ae = 0;
Pe = 0;


%==========================INITIAL CONDITIONS=============================

%Initial latitude and longitude (in degrees):
    %CAPE CANAVERAL AIR FORCE STATION, FL, USA: 28N 280W
    %A'MHOINE PENINSULA SCOTLAND: 58.5N -4.5W
    %KOUROU: 5N -52W

lat0 = 28;
long0 = 280;

x0 = Re * cosd(long0)*cosd(lat0);                                           %Initial x pos
y0 = Re * sind(long0)*cosd(lat0);                                           %Initial y pos
z0 = Re * sind(lat0);                                                       %Initial z pos
vx0 = 0;
vy0 = 0;
vz0 = 0;
ax0 = 0;
ay0 = 0;
az0 = 0;
Y0 = [x0;y0;z0;vx0;vy0;vz0;ax0;ay0;az0];

%Calculate burn time:
tb = MProp/mdotTot;                                                         %Stage 1 burn time
tb2 = Mprop2/mdotTot2;                                                      %Stage 2 burn time
IntegrationTime = tb+339.75;                                                   %Integrate over powered flight
dt = 0.1;
ts = [0, IntegrationTime]; 

%=========================SOLVE FOR DYNAMICS===============================

%   Solve using ODE45:
%[t, State] = ode45('EqOfMotion', ts, Y0);

State = EULER(0,x0,y0,z0,vx0,vy0,vz0,ax0,ay0,az0,IntegrationTime,dt);
St = cell2mat(State(1:10));

%Index solution array for state variables:
Time = St(1,:);
X = St(2,:);
Y = St(3,:);
Z = St(4,:);
Vx = St(5,:);
Vy = St(6,:);
Vz = St(7,:);
Ax = St(8,:);
Ay = St(9,:);
Az = St(10,:);

R_S = [X;Y;Z];
V_S = [Vx;Vy;Vz];

%Index for vehicle variables passed through the integrator
Veh = cell2mat(State(11));
Thrust = Veh(:,1);
Drag = Veh(:,2);
Mass = Veh(:,3);
Gamma = Veh(:,4);
Theta = Veh(:,5);
Psi = Veh(:,6);
Heading = Veh(:,7);

%Calculate position and inertial velocity in the Earth frame:
Mod_R_Pos = sqrt(X.^2 + Y.^2 + Z.^2);
Mod_V = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
Mod_A = sqrt(Ax.^2 + Ay.^2 + Az.^2);
Alt = (Mod_R_Pos - Re)/1E3;                                                 %Altitude in km

%Compute radial and tangential velocity:

%V_r = Mod_V .* sin(Gamma);
%V_t = Mod_V .* cos(Gamma);


%Compute the latitude and logitude along each point in the trajectory:

Lat = asin(Z./Mod_R_Pos);
Long = atan2(Y,X);

%Use this to compute the arc length in the longitude and latitude, and the
%square of this will be the actual trajectory downrange from the launch
%site.

Lat_Range = Mod_R_Pos.*(Lat-Lat(1));
Long_Range = Mod_R_Pos.*(Long - Long(1));

%Downrange distance in KM:
DownR = sqrt(Lat_Range.^2 + Long_Range.^2)./1000;

ORB = State2Orb(R_S(:,end), V_S(:,end));

inr = [x0;y0;z0];
%LaunchPadFrame(X,Y,Z,inr,lat0,long0)

%PLOTTING: 

figure;
hold on
earth_sphere('m');
plot3(X,Y,Z,'r','LineWidth',2)
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');

% figure;
% plot(Time,Alt,'k');
% title('Altitude (km)');
% 
% figure;
% plot(Time,Mod_V,'k',Time,V_r,'r');

figure;
%plot(Time(1:end-1),Gamma*180/pi,'k',Time(1:end-1),Heading*180/pi,'k--');
plot(Time(1:end-1),Gamma*180/pi,'k');
%legend('Flight Path','Heading');
title('Flight Angles');

figure;
plot(Time(1:end-1),Theta*180/pi,'k')
title('Theta');

figure;
plot(Time,Alt,'k');
title('Altitude v time (km)');

%Print terminal state to the command window:

fprintf("Terminal Altitude:             %6.2f km \n", Alt(end));
fprintf("Terminal Downrange:            %6.2f km \n", DownR(end));
fprintf("Terminal Velocity:             %6.2f m/s \n", Mod_V(end));
fprintf("Terminal Flight Path Angle:    %6.2f deg \n",Gamma(end)*180/pi);
fprintf("Time of flight:                %6.2f s \n", Time(end));