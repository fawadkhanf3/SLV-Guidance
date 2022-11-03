function [R, V, dtc, xc] = CSE(R0,V0,Tgo,x,dtc,xc)
% NASA CONIC STATE EXTRAPOLATION ROUTINE: Given a state vector on a conic
% trajectory and a time interval the routine can propagate the state
% forwards or backwards for the specified time. 
%
%The required inputs are:
% 
%   Initial Position and Velocity
%   Time to go
%   Value of x
%   Converged values of time intervals and x from last iteration
%
% Returns:
%   New position and velocity
%   Converged values for dT and x
%
global mu

    x = xc;
    
    % Set dtc to Tgo if first run
    if dtc == 0
        dtc = Tgo;
    end

    f0 = 1;
    n = 0;
    r0 = norm(R0);
    f1 = f0*sqrt(r0/mu);
    f2 = 1/f1;
    f3 = f2/r0;
    f4 = f1*r0;
    f5 = f0/sqrt(r0);
    f6 = f0*sqrt(r0);
    
    ir0 = R0/norm(R0);
    V0til = f1*V0;
    sig0til = dot(ir0, V0til);
    
    b0 = dot(V0til, V0til) - 1;
    alpha = 1 - b0;
    
    xg = f5*x;
    xlast = f5*xc;
    xmin = 0; 
    xmax = 2*pi/sqrt(abs(alpha));
    dtil = f3*Tgo;
    dtlast = f3*dtc; 
    dtmin = 0; 
    xp=0;
    P=0;
    
    if alpha > 0
    
        dtmax = xmax/alpha;
        xp = xmax;
        P = dtmax;
        
        while dtil >= P
        
            n = n + 1;
            dtil = dtil - P;
            dtlast = dtlast - P;
            xg = xg - xp;
            xlast = xlast - xp;
            
        end
    
    else 
        
        % CALL KEPLER TRANSFER TIME INTERVAL 
        [dtmax, A, D, E] = KeplerTT(xmax, sig0til, alpha);
        %if dtmax < dtil
        while dtmax < dtil

            dtmin = dtmax;
            xmin = xmax;
            xmax = 2*xmax;
            [dtmax, A, D, E] = KeplerTT(xmax, sig0til, alpha);
        end
        %end
    end
    
    if xg <= xmin || xg >= xmax
        
        xg = 0.5*(xmin + xmax);
    end
    
    [dtguess, A, D, E] = KeplerTT(xg, sig0til, alpha);
    
    if dtil < dtguess
       
        if xg < xlast && xlast < xmax && dtguess<dtlast && dtlast<dtmax
            
            xmax = xlast;
            dtmax = dtlast;
        
        end
        
    else 
        
        if xmin<xlast && xlast<xg && dtmin<dtlast && dtlast<dtguess
            
            xmin = xlast;
            dtmin = dtlast;
        end
        
    end
    
    % CALL KEPLER ITERATION LOOP 
    
    [xg,dtguess,A,D,E] = KIL(dtil, dtguess, xg, xmin, dtmin, xmax, dtmax, sig0til, alpha,A,D,E);
    
    rtil = 1 + 2*(b0*A + sig0til*D*E);
    b4 = 1/rtil;
    
    % CONVERGED VALUES OF X AND T
    
    xc = f6 * (xg+n*xp);
    dtc = f4*(dtguess + n*P);
    
    % LASTLY, EXTRAPOLATE THE STATE VECTOR
    
    F = 1-2*A;
    Gtil = 2*(D*E + sig0til*A);
    Ftilt = -2*b4*D*E;
    Gt = 1 - 2*b4*A;
    
    R = r0*(F*ir0 + Gtil*V0til);
    V = f2*(Ftilt*ir0 + Gt*V0til);
    
    
end

% SUBROUTINES USED IN THE MAIN LOOP, CONTAINS:
%
%   KEPLER TRANSFER TIME INTERVAL
%   U1 SERIES SUMMATION
%   Q CONTINUED FRACTION
%   KEPLER ITERATION LOOP
%   SECANT ITERATOR


function [dtarg, A, D, E] = KeplerTT(xarg, sigtil0, alpha)

% KEPLER TRANSFER TIME INTERVAL ROUTINE 

    u1 = U1(xarg, alpha);
    
    ztil = 2*u1;
    E = 1 - 0.5*alpha*ztil^2;
    w = sqrt(0.5*(1+E));
    
    D = w * ztil;
    A = D^2;
    B = 2*(E + sigtil0*D);
    
    Q = Qfrac(w);
    
    dtarg = D*(B + A*Q);

end

function u1 = U1(xarg, alpha)

% U1 SERIES SUMMATION
%
% Iterate and sum u1 until there is no change and return this variable
%
    du1 = xarg/4;
    u1 = du1;
    f7 = -alpha * du1^2;
    kmax = 100;
    k = 3;
    
    while k < kmax
        
        du1 = f7*du1/(k*(k-1));
        u1old = u1;
        u1 = u1 + du1;        
        
        if u1 == u1old
            
            break;
            
        end
        k = k+2;
    end
end


function Q = Qfrac(w)
        
    if w < 1

        xq = 21.04 - 13.04*w;

    elseif w < 4.625

        xq = 5/3 * (2*w+5);

    elseif w < 13.846

        xq = 10/7 *(w+12);

    elseif w < 44 

        xq = 0.5*(w+60);

    elseif w < 100

        xq = 0.25*(w+164);

    else 
        xq = 70;

    end

    b = 0;
    y = (w-1)/(w+1);

    j = floor(xq);

    b = y/(1 + ((j-1)/(j+1))*(1-b));

    while j > 2

        j = j - 1;
        b = y/(1 + ((j-1)/(j+1))*(1-b));

    end

    Q = 1/w^2 * (1 + 2*(1-b/4)/(3*w*(w+1)));

end


function [xguess,dtguess,A,D,E] = KIL(dtil, dtguess, xguess, xmin, dtmin, xmax, dtmax, sig0til, alpha,A,D,E)

% KEPLER ITERATION LOOP ROUTINE

    imax = 1000;
    i = 1;

    while i < imax
    
        dterr = dtil - dtguess;
    
        if abs(dterr) < 1E-6
            break;
        end
        
        % CALL SECANT INTEGRATOR 
        
        [dxtil, xmin, xmax, dtmin, dtmax] = SI(dterr, xguess, dtguess, xmin, xmax, dtmin, dtmax);

        xold = xguess;
        xguess = xguess + dxtil;
        
        if xguess == xold
            break;
        end
        
        dtold = dtguess;
        [dtguess, A, D, E] = KeplerTT(xguess, sig0til, alpha);
        
        if dtguess == dtold
            break;
        end
        
        i = i + 1;
    end

end


function [dxtil, xmin, xmax, dtmin, dtmax] = SI(dterr, xg, dtguess, xmin, xmax, dtmin, dtmax)

% SECANT INTEGRATOR ROUTINE #

    dtminp = dtguess - dtmin;
    dtmaxp = dtguess - dtmax;
    
    if abs(dtminp) < 1E-6 || abs(dtmaxp) < 1E-6
        
       dxtil = 0;
       
    else
        
        if dterr < 0
           
            dxtil = (xg-xmax)*dterr/dtmaxp;
            
            if xg + dxtil <= xmin
                
                dxtil = (xg-xmin)*dterr/dtminp;

            end
                
            xmax = xg;
            dtmax = dtguess; 
     
        else
            
            dxtil = (xg-xmin)*dterr/dtminp;
            
            if xg + dxtil >= xmax
                
                dxtil = (xg-xmax)*dterr/dtmaxp;
                
            end
            
            xmin = xg;
            dtmin = dtguess;                
                
        end 
    end
end


