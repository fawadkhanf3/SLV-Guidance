function [iF, Rd, Rbias, Rgrav, Vgo, Tgo, dtc, xc] = CycleUPFG(R, V, Acc, spre, Rbias0, Rgrav0, Rd0, Vgo0, Tgo0, T, M, Targ, dtc0, xc0)

% ROUTINE TO CYCLE UPFG ALGORITHM UNTIL CONVERGENCE OF THE VELOCITY TO BE
% GAINED. 
% The algorithm is initialised with a guess and repeatedly called without
% advancing the state vector of the vehicle, until the Vgo parameter
% converges to a constant value

    %Maximum cycles:
    max_iter = 100;
    cyclecount = 1;

    while cyclecount <= max_iter

        [iF, Rd, Rbias, Rgrav, Vgo, Tgo, dtc, xc] = UPFG(R, V, Acc, spre, cyclecount, Rbias0, Rgrav0, Rd0, Vgo0, Tgo0, T, M, Targ, dtc0, xc0);
        
        if norm(Vgo - Vgo0) < 5
            fprintf("GUIDANCE CONVERGED \n"); 
            fprintf("Estimated time to cutoff = %6.2f s\n", Tgo); 
            fprintf("Estimated Velocity to go = %6.2f km/s \n", norm(Vgo)/1e3); 
            break;
        end
        norm(Vgo);
        Rbias0 = Rbias;
        Rd0 = Rd;
        Rgrav0 = Rgrav;
        Vgo0 = Vgo;
        Tgo0 = Tgo;
        dtc0 = dtc;
        xc0 = xc;
        cyclecount = cyclecount+1;
    end
    
    if cyclecount == max_iter 
        
        fprintf("GUIDANCE NOT CONVERGED");
        
    end
    
end
