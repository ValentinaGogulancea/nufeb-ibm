function R = noDiffinteg(R)
% taking only what is needed out of the structure
Tol = R.Sxy.Tol;
dT = R.Sxy.dT; % time step for diffusion
dT_bac = R.Sxy.dT_bac;  % time step for bacteria mass balance
dT_Print = R.Sxy.dT_Print; % time step for printing
nTi = ceil (R.Sxy.maxT/dT_bac + 1); % number of time steps for diffusion
numStVLiq2 = R.St.numStVLiq2;  % number of liquid species
% 1st stage in computing the initial concentrations at each point in the
% domain - producing the matrices required for the difusion calculations


[R, ~, BB1, BB2, BB3, LL, UU] = DiffMatrices(R, 1); % ccommented further inside
R = boundary(R, dT, 0); % solves the reactor mass balance and assigns new boundary conditions
Sbc_Dir = 1000*kron(R.Sxy.Sbc_Dir(1:numStVLiq2), ones(R.Sxy.nT,1)); % concentrations in mol/L for the boundary at the top of the biofilm

NRVliq = R.rm.NRVliq;
U = R.Sxy.StVLiq2*1000; % concentrations of liquid species - mol/m3
V = R.Sxy.StVGas; % concentrations of gas species - bar
W = R.bac.bac_m;  % mass of bacteria - moles
U0 = U; % starting values

mm = 1000;
nn = 100;

j=2; kk=2; i = 2;  l = i + mm;
TimeBac = dT_bac*(j-1); TimePrint = dT_Print*(kk-1); Time = dT*(i-1); TimeSS = dT*(l-1); %TimePrev = Time;
numSt = numStVLiq2*R.Sxy.nT;
% make directory for the run
out_integTime(0, R, 0); % see comments there


while i < nTi
    % solving the diffusion reaction here
    % this is the 1 + Laplacian term - BB1*U
    % update for the boundary on the top - BB2*Sbc_Dir
    % the reaction term - dT*BB3*(1000)*NRVliq
    
    
    
    if Time >= TimeBac
        R = massBal( U/1000, V, W, R); % update the reaction matrices - for next diffusion step
        V = V + dT_bac*R.rm.NRVgas;  % update of solubles concentrations after biomass growth
        W = W + dT_bac*R.rm.NRVbac; % update of bacteria mass after biomass growth
        NRVliq = R.rm.NRVliq;
        U = U + dT_bac*NRVliq;
        dT_bac = R.Sxy.dT_bac; % adjusting time for division
        R = boundary(R, dT_bac, 1);  % updating the Dirichlet boundary condition for the diffusion-reaction
        
        fprintf('>>>> Current simulation time is %.1f hours.',Time)
        fprintf('Number of bacteria: %.0f', R.bac.bac_n)
        save('R.mat', 'R', '-v7.3');
        out_integTime(Time, R, i == nTi);
%         if Fdiv == 1
%             if R.bac.bac_n > R.bac.bac_nmax
%                 out_integTime(Time, R, 1);
%                 error('\n\n Maximum number of cells allowed exceeded !!\n\n Number of cells calculated: %d\n Maximum allowed %d\n\n\n >>> END of SIMULATION\n\n\n', R.bac.bac_n,R.bac.bac_nmax)
%             end
%             W = R.bac.bac_m;
%             NRVliq = R.rm.NRVliq;
%         end
        
    end
    j = j+1; TimeBac = dT_bac*(j-1);
    l = i + mm; TimeSS = dT*(l-1);
    
    fprintf('\n')
    
    
    
    i = i + 1; Time = dT_bac*(i-1);
end
end
function  out_integTime(t, R, fin)
nTSys = R.Sxy.nTSys; % number of grid cells for whole domain
nT = double(R.Sxy.nT); % number of grid cells for the biofilm + b-layer domain
indnT = R.Sxy.pos_xySys==1; % number of grid cell in which we have biofilm + b-layer

if t == 0
    All_StatesVar = [];
    B = [];
    
    if R.flagGas == 2
        % if gas-liquid mass transfer
        % retain also the gas flow rate
        
        All_StatesVar(end+1, :) = [0; R.St.StV; R.Qgas; (1/R.pOp.invHRT); sum(R.Sxy.Sbc_Dir(1:3))]'; %#ok<*NASGU>
    else
        % hold the concentrations on the boundary, HRT,
        All_StatesVar(end+1, :) = [0; R.St.StV; (1/R.pOp.invHRT); sum(R.Sxy.Sbc_Dir(1:3))]';  % holds the concentrations of all species over the entire computational domain
    end
    
    S = zeros(nTSys,R.St.numStVLiq2);
    ind = [nT*((1:R.St.numStVLiq2)'-1)+1,nT*(1:R.St.numStVLiq2)'];
    for k = 1:R.St.numStVLiq2
        A = R.Sxy.Sbc_Dir(k)*ones(nTSys,1); % values in the bulk liquid
        A(indnT) = R.Sxy.StVLiq2(ind(k,1):ind(k,2)); % overwritten with the concentrations in the BF+BL
        S(:, k) = A; % for one liquid species at a time
    end
    Sxy(:, :, 1) = S; % this holds the concentrations for each species at each grid cell at each time step recorded
    
    A = R.Sxy.pHnT;
    A(indnT) = R.Sxy.pH;
    pH(:, 1) = A; % pH values for time step
    
    DCat = zeros(nTSys,R.St.numX); DAn = DCat;
    ind = [nT*((1:R.St.numX)'-1)+1,nT*(1:R.St.numX)'];
    for k = 1:R.St.numX
        A = R.DGCat_bulk(k)*ones(nTSys,1);
        A(indnT) = R.Sxy.DGcat(ind(k,1):ind(k,2));
        DCat(:, k) = A;
        A =R.DGAn_bulk(k)*ones(nTSys,1);
        A(indnT) = R.Sxy.DGan(ind(k,1):ind(k,2));
        DAn(:, k) = A;
    end
    DGRCat(:, :, 1) = DCat; % computed values for catabolism free energy at  each grid cell at each time step
    DGRAn(:, :, 1) = DAn; % computed values for anabolism free energy at  each grid cell at each time step
else
    load ResultsSim.mat;
    bac = R.bac;
    B1 = t(end)*ones(bac.bac_nmax,1);
    B2 = [bac.bac_x;zeros(bac.bac_nmax-length(bac.bac_x),1)]; % coordinates of the existing bacteria -Ox
    B3 = [bac.bac_y;zeros(bac.bac_nmax-length(bac.bac_y),1)]; % see above - Oy
    B4 = [bac.bac_s;zeros(bac.bac_nmax-length(bac.bac_s),1)]; % types of cells
    B5 = [bac.bac_r;zeros(bac.bac_nmax-length(bac.bac_r),1)]; % radius of cells
    B(:, :, end+1) = [B1 B2 B3 B4 B5];
    if R.flagGas == 2
        All_StatesVar(end+1, :) = [t(end); R.St.StV; R.Qgas; (1/R.pOp.invHRT); sum(R.Sxy.Sbc_Dir(1:3))]';
    else
        All_StatesVar(end+1, :) = [t(end); R.St.StV; (1/R.pOp.invHRT); sum(R.Sxy.Sbc_Dir(1:3))]';
    end
    
    S = zeros(nTSys,R.St.numStVLiq2);
    ind = [nT*((1:R.St.numStVLiq2)'-1)+1,nT*(1:R.St.numStVLiq2)'];
    for k = 1:R.St.numStVLiq2
        A = R.Sxy.Sbc_Dir(k)*ones(nTSys,1);
        A(indnT) = R.Sxy.StVLiq2(ind(k,1):ind(k,2));
        S(:, k) = A;
    end
    
    Sxy(:, :, end+1) = S;
    
    A = R.Sxy.pHnT;
    A(indnT) = R.Sxy.pH;
    pH(:, end+1) = A;
    
    DCat = zeros(nTSys,R.St.numX); DAn = DCat;
    ind = [nT*((1:R.St.numX)'-1)+1,nT*(1:R.St.numX)'];
    for k = 1:R.St.numX
        A = R.DGCat_bulk(k)*ones(nTSys,1);
        A(indnT) = R.Sxy.DGcat(ind(k,1):ind(k,2));
        DCat(:, k) = A;
        A =R.DGAn_bulk(k)*ones(nTSys,1);
        A(indnT) = R.Sxy.DGan(ind(k,1):ind(k,2));
        DAn(:, k) = A;
    end
    DGRCat(:, :, end+1) = DCat;
    DGRAn(:, :, end+1) = DAn;
    
    if fin == 1
        toc
    end
end
save('ResultsSim.mat', 'All_StatesVar', 'B', 'DGRAn', 'DGRCat', 'pH', 'Sxy', '-v7.3');
end