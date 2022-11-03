function STATE = State2Orb(R,V)
%Given, a state vector of a spacecraft in inertial coordinates, compute the
%7 orbital elements:
%
%   h:      Specific orbital momentum [km2/s]
%   a:      Semi-major axis [km]
%   e:      Eccentricity
%   Om:     Right-ascension [deg]
%   i:      Inclination [deg]
%   w:      Argument of periapse [deg]
%   theta:  True anomaly [deg]
global Re mu
%Initialise values:

h_mag = 0;
a = 0;
e = 0;
Om = 0;
i = 0;
w = 0;
theta = 0;


%mu = 398600;                                                                %Standard gravitational parameter [km3/s2]

%Find magnitudes of the velocity and position
r = norm(R);
v = norm(V);

vr = dot(R,V)/r;

%Specific orbital momentum:
h = cross(R,V);
h_mag = norm(h);

%Inclination:
i = acosd(h(3)/h_mag);

K = [0,0,1];                                                                %This defines the unit vector in the Z direction

N = cross(K,h);
N_mag = norm(N);

%Right-ascension:

Om = acosd(N(1)/N_mag);

if N(2) < 0
    
    Om = 360 - Om;

end

%Eccentricity:

e = 1/mu * (cross(V,h) - mu * R/r);
e_mag = norm(e);

%Argument of periapsis:

w = acosd(dot(N,e)/(N_mag*e_mag));

if e(3) < 0
    
    w = 360 - w;

end

%True anomaly:

theta = acosd(dot(e,R)/(e_mag*r));

if vr < 0
    
    theta = 360 - theta;
    
end

rp = h_mag^2/mu * 1/(1+e_mag*cosd(0));
ra = h_mag^2/mu * 1/(1+e_mag*cosd(180));

a = 0.5*(ra+rp);

fprintf("Terminal orbital state\n");
fprintf("Periapsis:         %6.2f km\n", (rp-Re)/1000);
fprintf("Apoapsis:          %6.2f km\n", (ra-Re)/1000);
fprintf("Eccentricity:      %6.2f\n", e_mag)
fprintf("RAAN:              %6.2f deg\n", Om);
fprintf("Inclination:       %6.2f deg\n", i);


STATE = [h_mag,a,e_mag,Om,i,w,theta];

end

