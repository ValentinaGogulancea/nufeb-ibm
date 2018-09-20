function R = integTimeTest(R)
% taking only what is needed out of the structure
Tol = R.Sxy.Tol;
dT = R.Sxy.dT; % time step for diffusion
dT_bac = R.Sxy.dT_bac;  % time step for bacteria mass balance
dT_Print = R.Sxy.dT_Print; % time step for printing
nTDiff =ceil (R.Sxy.maxT/dT + 1); % number of time steps for diffusion
nTBac = ceil (R.Sxy.maxT/dT_bac + 1);
numStVLiq2 = R.St.numStVLiq2;  % number of liquid species
% 1st stage in computing the initial concentrations at each point in the
% domain - producing the matrices required for the difusion calculations


[R, ~, BB1, BB2, BB3, LL, UU] = DiffMatrices(R, 1); % ccommented further inside
Sbc_Dir = 1000*kron(R.Sxy.Sbc_Dir(1:numStVLiq2), ones(R.Sxy.nT,1)); % concentrations in mol/L for the boundary at the top of the biofilm

NRVliq = R.rm.NRVliq;
U = R.Sxy.StVLiq2*1000; % concentrations of liquid species - mol/L
V = R.Sxy.StVGas; % concentrations of gas species - bar
W = R.bac.bac_m;  % mass of bacteria - moles
U0 = U; % starting values


i = 1;
j=1;

numSt = numStVLiq2*R.Sxy.nT;
% make directory for the run
out_integTime(0, R, 0); % see comments there

SS = 2*Tol;
% SS= 0;
while i< nTBac
    j = 1;
    while SS > Tol
        % solving the diffusion reaction here
        % this is the 1 + Laplacian term - BB1*U
        % update for the boundary on the top - BB2*Sbc_Dir
        % the reaction term - dT*BB3*(1000)*NRVliq
        
        B = BB1*U  + BB2*Sbc_Dir + dT*BB3*(1000)*NRVliq;
        y = LL\B;
        U = UU\y;
        U = (U > 0).*U;
        R.Sxy.StVLiq2 = U/1000;
        
        j = j + 1;
        
        if mod(j,1000) == 0
            [NRVbac, NRV, ~, ~, ~, ~] = my_kinetics(R);  % update the rates
            NRVliq = NRV(1:numSt); R.rm.NRVliq = NRVliq;
            NRVgas = NRV(numSt+1:end); R.rm.NRVgas = NRVgas;
            R.rm.NRVbac = NRVbac;
            
            V = V + dT*NRVgas;  R.Sxy.StVGas = V;
            W = W + dT*NRVbac; R.bac.bac_m = W; % updating biomass concentration
            R.Sxy.StVLiq = [U/1000; V];
        end

        SS = 100*max(abs(U - U0)./U0);
        
        U0 = U;
        if j > dT_bac/dT
            fprintf('no steady state diffusion')
            return
        end
        
    end
   
    SS = 2 * Tol;
    dT_bac = dT_bac - j * dT;
    R = massBal( U/1000, V, W, R); % update the reaction matrices - for next diffusion step
    V = V + dT_bac*R.rm.NRVgas;  % update of solubles concentrations after biomass growth
    W = W + dT_bac*R.rm.NRVbac; % update of bacteria mass after biomass growth
    NRVliq = R.rm.NRVliq;
    dT_bac = R.Sxy.dT_bac; % adjusting time for division
    R = boundary(R, dT_bac, 1);  % updating the Dirichlet boundary condition for the diffusion-reaction
    
    fprintf('>>>> Current simulation time is %.1f hours.',dT_bac*i)
    if mod(i,dT_Print/dT_bac) == 0
        [R, Fdiv] = bacteria( U/1000, V, W, R); % this calls the shoving algorithm
        [R, Fb_layer, BB1, BB2, BB3, LL, UU] = DiffMatrices(R, Fdiv); % computing the matrices again
        fprintf('Number of bacteria: %.0f', R.bac.bac_n)
         fprintf('\n')
         fprintf('Mass of bacteria: %.6e', sum(R.bac.bac_m)*24.6)
        save('R.mat', 'R', '-v7.3');
        out_integTime(dT_bac*i, R, j == nTDiff);
%        
        if Fdiv == 1
            if R.bac.bac_n > R.bac.bac_nmax
                out_integTime(Time, R, 1);
                error('\n\n Maximum number of cells allowed exceeded !!\n\n Number of cells calculated: %d\n Maximum allowed %d\n\n\n >>> END of SIMULATION\n\n\n', R.bac.bac_n,R.bac.bac_nmax)
            end
            W = R.bac.bac_m;
            NRVliq = R.rm.NRVliq;
        end
        if Fb_layer == 1
            U = R.Sxy.StVLiq2*1000;
            U0 = U;
            V = R.Sxy.StVGas;
            numSt = numStVLiq2*R.Sxy.nT;
        end
        
    end
    
    
    Sbc_Dir = 1000*kron(R.Sxy.Sbc_Dir(1:numStVLiq2), ones(R.Sxy.nT,1));
    fprintf('\n')
    i = i + 1;
end
end





