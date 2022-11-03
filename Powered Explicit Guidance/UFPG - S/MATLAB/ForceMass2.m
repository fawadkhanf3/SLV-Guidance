function V_Var = ForceMass2(t,r,V)
%Computes the thrust T, mass M and drag D at time t for the launch vehicle

global M0 mdotProp mdotTot ue S Re tb tb2 M02 ue2 mdotTot2 gamma...
mdotProp2 

%Compute the local atmospheric conditions:
Atm = Atmosphere(r-Re);

Pr = Atm(:,1);
rho = Atm(:,2);
Temp = Atm(:,3);


%Calculate the Mach number, and use this to find CD:

A = (gamma*287*Temp)^0.5;                                                   %The speed of sound [m/s]
%V = (p^2 + (r*q)^2)^0.5;                                                    %Absolute velocity of vehicle [m/s]
Mach = norm(V)/A;                                                                 %The Mach number 

%Drag model approximated using an Atlas V Cd v M curve interpolated in 4
%segments:

if Mach <= 1.25
    
    CD = -0.0415*Mach^3 + 0.3892*Mach^2 - 0.2614*Mach + 0.303;
    
elseif 1.25 < Mach && Mach <= 4
    
    CD = -0.049*Mach^4 + 0.5664*Mach^3 - 2.3265*Mach^2 + 3.8512*Mach + ...
        -1.6625;
    
elseif 4 < Mach && Mach <=10
    
    CD = -0.0037*Mach^3 + 0.0695*Mach^2 - 0.4105*Mach + 0.9732;
    
elseif Mach > 10
    
    CD = 0.255;
    
end


% CD = 0.5;
D = 0.5 * rho * (norm(V)) * S * CD * V;                                             %The drag force [N]
M = 0;
T = 0;

if t < tb
    T = (mdotProp*ue);                                                     %1st stage thrust accounting for ALL 5 engines [N]
    M = M0 - mdotTot*t;                                                     %Mass @ time t accounting for ALL 5 engines [kg]

elseif  t > tb && t <= tb+tb2                                                %After stage 1 MECO and up to 2nd stage termination
    T = (mdotProp2*ue2)*1;                                                  %2nd stage thrust [N]
    M = M02 - mdotTot2*(t-tb);                                              %2nd stage Mass @ time t  [kg]

elseif t > (tb+tb2) 
    T = 0;                                                                  %MECO3, end of powered flight
    M = M02 - (mdotTot2*tb2);                                               %The burnout mass at the end of powered flight
end


%Output the thrust, drag, and mass into one matrix:

V_Var = {M,T,D};

end

