function A = Atmosphere(h)

%Computes the pressure, density and temperature at a geopotential altitude h (m) above
%the ground using a standard atmospheric model.

%Define constants for each segment of the atmosphere:
re = 6378E3;                                                                %Earth radius [m]
g0 = 9.81;                                                                  %Gravitational acceleration [m/s^2]
R = 287;                                                                    %Gas constant [J/kgK]

%First linear segment 0-11 km
T1 = 288.16;                                                                %Sea level temp
P1 = 101325;                                                                %Sea level pressure
rho1 = 1.225;                                                               %Sea level density
a1 = -6.5E-3;                                                               %Lapse rate [K/m]

%First isothermal segment 11-25 km
T2 = 216.66;
P2 = 2.2616E4;
rho2 = 0.3636;

%Second linear segment 25-47 km
T3 = 216.66;
P3 = 2.484E3;
rho3 = 0.0399;
a3 = 3E-3;

%Second isothermal segment 47-53 km
T4 = 282.66;
P4 = 120.0438;
rho4 = 0.0015;

%Third linear segment 53-79 km
T5 = 282.66;
P5 = 58.1075;
rho5 = 7.2608E-4;
a5 = -4.5E-3;

%Convert geometric altitude to geopotential altitude:

h = re*h/(re + h);

if h < 0
    fprintf("crash");
    
end


if h <= 11E3

    T = T1 + a1*h;
    Pr = P1 * (T/T1)^(-g0/(a1*R));
    rho = rho1 * (T/T1)^-(g0/(a1*R) + 1);

elseif h > 11E3 && h <= 25E3

    T = T2;
    Pr = P2 * exp(-g0*(h-11E3)/(R*T2));
    rho = rho2 * exp(-g0*(h-11E3)/(R*T2));

elseif h > 25E3 && h <= 47E3

    T = T3 + a3*(h - 25E3);
    Pr = P3 * (T/T3)^(-g0/(a3*R));
    rho = rho3 * (T/T3)^-(g0/(a3*R) + 1);

elseif h > 47E3 && h <= 53E3 
    
    T = T4;
    Pr = P4 * exp(-g0*(h-47E3)/(R*T4));
    rho = rho4 * exp(-g0*(h-47E3)/(R*T4));

elseif h > 53E3 && h <= 79E3
    
    T = T5 + a5*(h - 53E3);
    Pr = P5 * (T/T5)^(-g0/(a5*R));
    rho = rho5 * (T/T5)^-(g0/(a5*R) + 1);     
    
elseif h > 79E3
    
    T = 165.66;
    Pr = 0;
    rho = 0;
end
A = [Pr, rho, T];                                                           %Output a matrix of pressure, density and temp at the requested h
end

