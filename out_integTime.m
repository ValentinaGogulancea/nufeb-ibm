function  out_integTime(t, R, fin,PathOutput)
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
        S(:,k) = A; % for one liquid species at a time
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
    B6 = [bac.bac_m;zeros(bac.bac_nmax-length(bac.bac_m),1)]; % mass of cells
    B(:, :, end+1) = [B1 B2 B3 B4 B5 B6];
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
        A = R.DGAn_bulk(k)*ones(nTSys,1);
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

[SUCCESS,MESSAGE,MESSAGEID] = copyfile( 'ResultsSim.mat', PathOutput,'f');
end