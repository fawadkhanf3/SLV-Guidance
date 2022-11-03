function S = EqOfMotion(t,Y,UPFG_in)
%Equations of motion in 3 degrees of freedom describing the trajectory of
%the centre of mass of the launch vehicle. Input the cartesian position and
%velocity and outputs the vector dVdt and drdt.
%   dvdt = 1/m * [ T_BE * F_b] + F_grav
%   drdt = V
global mu Re dt tb tb2 ue2
persistent theta psi countUPFG Tgo0


if isempty(theta)
    
    theta = 0;

end

if isempty(countUPFG)
    countUPFG = 0;
end

x = Y(1);
y = Y(2);
z = Y(3);
vx = Y(4);
vy = Y(5);
vz = Y(6);
ax = Y(7);
ay = Y(8);
az = Y(9);
UPFG_vars = UPFG_in;

% Position and velocity state vectors:

R = [x;y;z];
V = [vx;vy;vz];
Acc = [ax;ay;az];

%Calculate the longitude, latitdue (in rad) and radial position to compute the
%position in geographic NED coordinates:

r = sqrt(x^2+y^2+z^2);

v_i = [vx;vy;vz];
V_i = sqrt(vx^2+vy^2+vz^2);

long = atan2(y,x);
lat = asin(z/r);


%Define the transformation matrix for Earth inertial - local:


T_EG = [cos(lat)*cos(long), cos(lat)*sin(long), sin(lat);
    -sin(long), cos(long), 0;
    -sin(lat)*cos(long), -sin(lat)*sin(long), cos(lat)];

invTEG = inv(T_EG);


%Transform the velocity coordinates to local to compute the velocity
%vector:= (Up-East-North)

V_geo =  T_EG * [vx;vy;vz];


%Flight path and heading angles:

gamm = asin(V_geo(1)/V_i);
hdg = atan2(V_geo(2),V_geo(3));
hdg*180/pi;
%[gamm*180/pi,hdg*180/pi]

%Heading angle:

V_v = ForceMass2(t,r,v_i);

M = V_v{1};
T = V_v{2};
D = V_v{3};

    
%GRAVITY TURN:

    %After the vehicle is approx. 200m high, command a constant thrust
    %angle (~1-1.5 deg) in order to begin inducing a lateral velocity,
    %until some value of flight path angle is desired from which the
    %gravity turn can begin. After this, command the thrust vector to point
    %in the direction of the flight path angle to ensure a zero-lift
    %trajectory. 

%theta = 0*pi/180;
psi = 90*pi/180;

if t < tb
    if r - Re > 150

        if gamm > 86*pi/180

            %1 deg/second change in the thrust vector until the flight path
            %angle is at a suitable angle to begin the gravity turn:

            %theta = (theta + dt)*pi/180;
            theta = 1.7*pi/180;

        elseif gamm <= 86*pi/180

            theta = pi/2 - gamm;

        end
    end
end

% CALL UNIFIED POWERED FLIGHT GUIDANCE ROUTINE 

% Declare parameters
R0 = 0;
V0 = 0;
iy = 0;
rd = 0;
vd = 0;
gammad = 0;

%UPFG variable initialisation:
Rd = 0;
Rbias = 0;
Rgrav = 0;
Vd = 0;
Vgo = 0;
Tgo = 0;
spre = 0;
csetc = 0; csexc = 0; dtc = 0; xc = 0;
iF = 0;
Targ = {iy,rd,vd,gammad};

if t > tb
    
    countUPFG = countUPFG + 1;
    rd = Re + 300e3;
    vd = sqrt(mu/rd);
    gammad = 0;                                                        % Measured from LOCAL HORIZONTAL 
    h = rd*vd;
    [R0, V0] = Orb2State(h,0,deg2rad(29),0,0,0);
    %Target normal - opposite to orbital momentum vector
    iy = -cross(R0,V0)/norm(cross(R0,V0));
    Targ = {iy,rd,vd,gammad};

    if countUPFG == 1
    
        % TARGET CONDITIONS: RD AND VD
        
        % Initial guess values to cycle routine 
        
        Rd = rd*Rodrigues((R/norm(R)), -iy, 20);
        Rbias = [0;0;0];
        Vd = vd*(cross(-iy,Rd)/norm(cross(-iy,Rd))) - V;
        Tgo = tb2;
        spre = 1;
        csetc = 0;
        csexc = 0;
        [iF, Rd, Rbias, Rgrav, Vgo, Tgo, dtc, xc] = CycleUPFG(R, V, Acc, spre, Rbias, 0, Rd, Vd, Tgo, T, M, Targ, csetc, csexc);

    else
    
    % After pre-cycling, enter main routine to advance vehicle state and
    % issue steering commands:

        % Output from last UPFG cycle 
        spre = 0;
        Rd = UPFG_in{1}{1};
        Rbias = UPFG_in{1}{2};
        Rgrav = UPFG_in{1}{3};
        Vgo = UPFG_in{1}{4};
        Tgo = UPFG_in{1}{5};
        dtc = UPFG_in{1}{6};
        xc = UPFG_in{1}{7};

        %Tgo
        [iF, Rd, Rbias, Rgrav, Vgo, Tgo, dtc, xc] = UPFG(R, V, Acc, spre, 0, Rbias, Rgrav, Rd, Vgo, Tgo, T, M, Targ, dtc, xc);
    end
    % IF FIRST ITERATION, COMPUTE THE TARGET CONDITIONS AND SET INITIAL GUESS
    % FOR RD/VGO, THEN CYCLE GUIDANCE (PRETHRUST) UNTIL CONVERGENCE 
    
    % Transformation of thrust vector from inertial to local:
    
    iF_local = T_EG*iF;
    theta = acos(iF_local(1));
    psi = acos(iF_local(3)/sin(theta));
    
    % CUTOFF CONDITIONS
    norm(Vgo);
    if Tgo < 0.01
        
        T = 0;
        fprintf("Cutoff confirmed\n");
        TERMSTATE = State2Orb(R, V);
    end
    
end


%Define the unit vector in the local frame by specifying the Euler angles
%of pitch (theta) and yaw (psi) which orient the body axis (hence thrust
%vector), as it is assumed all thrust acts through the nose of the vehicle.
%This can then be transformed to the inertial frame to define the thrust
%direction. 

    %******UPDATE THIS TO MATCH THE DEFINITION OF GAMMA

    %Theta is defined as 0 when pointing up in the local frame and 90 when
    %against the horizon 

[theta*180/pi, psi*180/pi];
utlx = cos(theta);
utly = sin(theta)*sin(psi);
utlz = sin(theta)*cos(psi);

%Unit vector locally:
Utl = [utlx; utly; utlz];

%Unit vector inertially, transformed from local:

Uti = invTEG * Utl;
F_g = -mu/r^3 * [x; y; z];


%Finally, the equations of motion in inertial coordinates:

dVdt = 1/M * (T*Uti - D) + F_g;
drdt = [vx;vy;vz];

%Vector to store the vehicle variables. These will be passed through the
%integrator and output to validate nominal performance. 

V_variables = [T, norm(D), M, gamm, theta, psi, hdg];
UPFG_vars = {Rd, Rbias, Rgrav, Vgo, Tgo, dtc, xc};
S = {drdt;dVdt;V_variables;UPFG_vars};

end

