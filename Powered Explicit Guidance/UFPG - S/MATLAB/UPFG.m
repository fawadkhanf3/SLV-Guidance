function [iF, Rd, Rbias, Rgrav, Vgo, Tgo, dtc, xc] = UPFG(R, V, Acc, spre, cyclecount, Rbias, Rgrav, Rd, Vgo, Tgo0, T, M, Targ, dtc0, xc0)
%UNIFIED POWERED FLIGHT GUIDANCE MAIN PROGRAM ROUTINE
%
% Input: current position and velocity vectors as well as the last cycle
% variables for UPFG. Decision flag for prethrust cycling.
%
% Returns: Unit thrust vector and UPFG variables such as Rgo, Vgo, Tgo
%

global mu dt ue2


% Desired conditions:

iy = cell2mat(Targ(1));
rd = cell2mat(Targ(2));
vd = cell2mat(Targ(3));
gammaD = cell2mat(Targ(4));

% CSE input: 


% BLOCK 1 (INITIALIZATION):
% ONLY ENTER ON FIRST CALL DURING PRETHRUST CYCLING 

    if spre == 1 && cyclecount == 1
        Rgrav = -0.5 * mu * R / norm(R)^3;
    end
    
% BLOCK 2 (UPDATES):

    if spre == 0
    
        Vgo = Vgo - (Acc*dt);                                               % WRONG: needs acceleration to decrement Vgo: DONE 
        Tgo = Tgo0 - dt;
        
    end
    
% BLOCK 3 (TIME TO GO):

    a = T/M;
    tau = ue2/a;
    
    %L = -ue2 * log((tau-Tgo0)/tau);
    L = norm(Vgo);
    tb = tau * (1- exp(-L/ue2));
    Tgo = tb;
    
% BLOCK 4 (THRUST INTEGRALS): 

    J = L*tau - ue2*tb;
    S = -J + L*tb;
    Q = S*tau - 0.5*ue2*tb^2;
    P = Q*tau - (1/6)*ue2*tb^3;
    
    %J = J + L*Tgo;
    %S = S + L*tb;
    %Q = Q + J*tb;
    H = J*Tgo - Q;

    %P = P + H*tb;
    
    
% BLOCK 5 (TURNING RATE):

    lambda = Vgo/norm(Vgo);
    Rgrav = (Tgo/Tgo0)^2 * Rgrav;
    Rgo = Rd - (R + V*Tgo + Rgrav);
    
    iz = cross(Rd, iy)/norm(cross(Rd, iy));                                 % RD AND IY DECLARATIONS - done
    
    Rgoxy = Rgo - dot(iz, Rgo)*iz;
    Rgoz = (S - dot(lambda, Rgoxy))/dot(lambda, iz);
    Rgo = Rgoxy + Rgoz*iz + Rbias;
    
    lambdadot = (Rgo - S*lambda)/(Q - (S*J/L));
    ulambdadot = lambdadot/norm(lambdadot);
    
    
    iF = (lambda - ((J/L) * lambdadot))/norm(lambda - (J/L * lambdadot));
    
    phi = acos(dot(iF, lambda));
    phidot = -phi*L/J;
    
    Vthr = (L - 0.5*L*phi^2 - J*phi*phidot - 0.5*H*phidot^2)*lambda - (L*phi + J*phidot)*ulambdadot;
    Rthr = (S- 0.5*S*phi^2 - Q*phi*phidot - 0.5*P*phidot^2)*lambda - (S*phi + Q*phidot)*ulambdadot;
    
    Vbias = Vgo - Vthr;
    Rbias = Rgo - Rthr;
    
% BLOCK 7 (GRAVITY PREDICTIONS): 

    dRc = -0.1*Rthr - ((1/30) * Vthr*Tgo);
    dVc = 1.2*Rthr/Tgo - (0.1 * Vthr);
    
    Rc1 = R + dRc;
    Vc1 = V + dVc;
    
    % CALL CSE ROUTINE HERE 
    
    [Rc2, Vc2, dtc, xc] = CSE(Rc1, Vc1, Tgo, 0, dtc0, xc0);                 % NEED LAST CSE CALL INPUTS - done 
    
    Vgrav = Vc2 - Vc1;
    Rgrav = Rc2 - Rc1 - Vc1*Tgo;
    
% BLOCK 8 (VELOCITY TO BE GAINED):

    Rp = R + V*Tgo + Rgrav + Rthr;
    Rd = rd * (Rp)/norm(Rp);                                                % NEED TO ADD INPUT rd - done
    
    ix = Rd/norm(Rd);
    iz = cross(ix, iy);
    
    
    %Vd = vd*[ix'*[sind(gammaD); 0; cosd(gammaD)]; iy'*[sind(gammaD); 0; cosd(gammaD)];iz'*[sind(gammaD); 0; cosd(gammaD)]];
    
    Vd = vd*(ix*sind(gammaD) + iz*cos(gammaD));
    
    Vgop = Vd - V - Vgrav + Vbias;
    Vgo = Vgop;

end
















function [iF, Rd, Rbias, Rgrav, Vgo, Tgo, dtc, xc] = CycleUPFG(R, V, spre, Rbias0, Rd0, Vgo0, Tgo0, T, M, Targ, dtc0, xc0)

% ROUTINE TO CYCLE UPFG ALGORITHM UNTIL CONVERGENCE OF THE VELOCITY TO BE
% GAINED. 
% The algorithm is initialised with a guess and repeatedly called without
% advancing the state vector of the vehicle, until the Vgo parameter
% converges to a constant value

    %Maximum cycles:
    max_iter = 1000;
    cyclecount = 1;

    while cyclecount < max_iter

        [iF, Rd, Rbias, Rgrav, Vgo, Tgo, dtc, xc] = UPFG(R, V, spre, cyclecount, Rbias0, Rd0, Vgo0, Tgo0, T, M, Targ, dtc0, xc0);
        
        if norm(Vgo - Vgo0) < 100
            break;
        end
        norm(Vgo-Vgo0)
        Rbias0 = Rd;
        Rd0 = Rd;
        Vgo0 = Vgo;
        Tgo0 = Tgo;
        dtc0 = dtc;
        xc0 = xc;
        cyclecount = cyclecount+1;
    end
end
