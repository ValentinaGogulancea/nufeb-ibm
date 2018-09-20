function [NRVbac, NRV, pH, DGCatM, DGAnM, Df,spcM] = my_kinetics(R)
% computes oh so many things
% NRVbac = rate of bacterial growth
% NRVxybac = rate of bacterial growth in each grid cell
% NRV = rates of soluble species consumption/formation
% R.Sxy.pH = pH values in each grid cell, same goes of DGcat & DGan
% Df - diffusion correction coefficients at each grid cell
% X - molar concentration of biomass (all bacteria)

Keq = R.pTh.Keq;
chrM = R.pTh.chrM;
flagpH = R.flagpH; % if pH is constant or not
flagDG = R.flagDG; % if yield is constant or not
MatrixCat = R.rm.MatrixCat_mod; % catabolic matrix that contains only the solubles and no H2O or H+
MatrixAn = R.rm.MatrixAn_mod; % anabolic ----''-------
MatrixDecay = R.rm.MatrixDecay_mod; % decay small matrix
rMatrixFull = R.rm.rMatrixFull'; % full species, without the decay tho
Vg = R.Sxy.Vg; % volume of one grid cell
nT = R.Sxy.nT; % number of grid cells
Sh_ini_all = 10.^(-R.Sxy.pH); % proton concentrations in comp. domain
SLiqxy_all = R.Sxy.StVLiq; % concentrations of all solubles in all grids
RthT = R.pOp.Rth*R.pOp.T; % the product of R and temperature - in kJ/mole
VRgT = R.pOp.Vgas/(R.pOp.Rg*R.pOp.T); % basically pV = nuRT---nu = p* V/R*T
KhV = R.pTh.Kh.KhV; % Henry's constants
KhV = 0;
kLa = R.kTr.kLa;
DG0 = R.pTh.DG0; % free energy of formation for all species
DGr0 = R.pTh.DGr0; % free energy of the cat+ ana
DGdis = R.pTh.DGdis; % dissipation energy
spcR = R.pTh.spcR; % matrix that specifies the types of charges the species posses in the reactions written
lCat_eD = R.pTh.lCat_eD; % stoichiometric coeff of eDonors

numStVLiq = R.St.numStVLiq; % all soluble
numStVLiq2 = R.St.numStVLiq2; % all liquid
numStVGas2Liq = R.St.numStVGas - 1; % all gas
numX = R.St.numX; % bacteria
numStFull = R.St.numStFull; % all components
Liq2Gas = R.St.Liq2Gas; % those from gas to liquid

Glu_pos = strcmp(R.St.StNames(1:numStVLiq), 'Glu');
AcH_pos = strcmp(R.St.StNames(1:numStVLiq), 'AcH');
H2_pos = strcmp(R.St.StNames(1:numStVLiq), 'H2');
% O2_pos = strcmp(R.St.StNames(1:numStVLiq), 'O2');
% NH3_pos = strcmp(R.St.StNames(1:numStVLiq), 'NH3');
% NO2_pos = strcmp(R.St.StNames(1:numStVLiq), 'NO2');
rv = R.pTh.react_v;

bac_s = R.bac.bac_s;
bac_m = R.bac.bac_m;
bac_mu_max = R.bac.bac_mu_max;
bac_yield = R.bac.bac_yield;
bac_Ks = R.bac.bac_Ks;
bac_maint = R.bac.bac_maint;
CC = R.Sxy.c; % matrix that tells if bacteria i is in grid cell j
% initialization
NRVbac = zeros(R.bac.bac_n, 1); NRV = zeros(numStVLiq*nT,1); Ngas = NRV; X = zeros(nT,1); pH = X;
DGCatM = zeros(nT*numX,1); DGAnM = DGCatM; DG_mat = ones(numStFull,1);

for i=1:nT
    c = CC(:,i); % all bacteria vs grid cell i - length of bacterial number
    
    indLiq = i+((1:numStVLiq)-1)*nT; % indices for where the i-th soluble species concentration begins
    indCat = i+((1:numX)-1)*nT; % same but for bacteria
    
    SLiqxy = SLiqxy_all(indLiq);  % values for the concentrations of all solubles  in grid i
    [spcM, Sh] = solve_pH(Sh_ini_all(indLiq(1)), [SLiqxy;1;0], Keq, chrM, flagpH); % commented inside
    pH(indLiq(1)) = -log10(Sh); % pH at every grid cell
    
    
    if R.flagGas == 2
        % if the reactor is closed and mass transfer is computed
        for jGas = 1:numStVGas2Liq
            kGas = Liq2Gas == jGas;
            jGas_aux = numStVLiq2 + jGas;
            if DG0(kGas,1)== 1e10
                GasT = -kLa(jGas)*(spcM(jGas_aux,2) - spcM(kGas,2)/KhV(jGas));
            else
                GasT = -kLa(jGas)*(spcM(jGas_aux,2) - spcM(kGas,1)/KhV(jGas));
            end
            Ngas(indLiq(jGas_aux)) = GasT;
            Ngas(indLiq(kGas)) = - GasT*VRgT; % partial pressures
        end
    end
    spcMc = log((spcM == 0)*1e-20 + spcM);
    DG_mat(1:(numStVLiq+2)) = sum(spcMc.*spcR, 2);
    DGr_mat = rMatrixFull*DG_mat;
    DGr = DGr0 + RthT*DGr_mat; % after pH calculations - compute the free energy corrections -
    % deviations from standard
    DGCatM(indCat) =  DGr(1:2:end);
    DGAnM(indCat) = DGr(2:2:end);
    bioX = sum(bac_m(c == 1));     % the moles of bacteria in grid cell i
    X(i) = bioX/Vg;  % molar concentration of biomass
    
    if ne(bioX, 0) % if the biomass concentration is not 0 compute the rates at every grid cell
        for k = 1: numX
            ex = find(ne(c.*(bac_s==k),0));  % check what type of bacteria we have in grid cell i and retain their position - i.e which of the cells it is
            if (1-isempty(ex)) == 1 % if there actually is a bacteria
                for e = 1:length(ex)  % there might be more than one bacteria of the same type
                    if flagDG == 2    % if the yields are not assumed constant
                        DGCataux = -DGr(1:2:end);
                        lCat = (DGr(2:2:end) + DGdis)./DGCataux + lCat_eD; % compute the yields as gS-gX
                        lCat(DGCataux < 0) = 0;
                        
                        InvYield = repmat(lCat(k), numStVLiq,1); % make it the length of the soluble species vector
                        if lCat(k) == 0
                            bac_yield(ex(e)) = 0;
                        else
                            bac_yield(ex(e)) = 1/lCat(k); % invert them as gX-gS
                        end
                        R.bac.bac_yield(ex(e)) = bac_yield(ex(e));
                    else
                        InvYield = repmat(1/bac_yield(ex(e)), numStVLiq,1); % if the yields are constant
                    end
                    MatrixMet = MatrixCat(:,k).*InvYield + MatrixAn(:,k); % assemble the metabolic matrix
                    
                    M = 1;
%                     % compute the rates for each bacteria now
%                     if ne(bac_Ks(ex(e),1),0)
%                         M = M*spcM(NH3_pos,rv(NH3_pos))/(spcM(NH3_pos,rv(NH3_pos)) + bac_Ks(ex(e),1));
% %                         use the proper concentration, corresponding to the proper state
%                     end
%                     if ne(bac_Ks(ex(e),2),0)
%                         M = M*(spcM(NO2_pos,rv(NO2_pos))/(spcM(NO2_pos,rv(NO2_pos)) + bac_Ks(ex(e),2)));
%                     end
%                     if ne(bac_Ks(ex(e),3),0)
%                         M = M*(spcM(O2_pos,rv(O2_pos))/(spcM(O2_pos,rv(O2_pos)) + bac_Ks(ex(e),3)));
%                     end
%                     
%                     M = M*spcM(Glu_pos,rv(Glu_pos))/(spcM(Glu_pos,rv(Glu_pos)) + bac_Ks(ex(e),1));

                    if ne(bac_Ks(ex(e),1),0)
                        M = M*spcM(Glu_pos,rv(Glu_pos))/(spcM(Glu_pos,rv(Glu_pos)) + bac_Ks(ex(e),1));
%                         use the proper concentration, corresponding to the proper state
                    end
                    if ne(bac_Ks(ex(e),2),0)
                        M = M*(spcM(H2_pos,rv(H2_pos))/(spcM(H2_pos,rv(H2_pos)) + bac_Ks(ex(e),2)));
                    end
                    if ne(bac_Ks(ex(e),3),0)
                        M = M*(spcM(AcH_pos,rv(AcH_pos))/(spcM(AcH_pos,rv(AcH_pos)) + bac_Ks(ex(e),3)));
                    end
                    mu = bac_mu_max(ex(e))*M - bac_maint(ex(e));
                    
                    if mu > 0  % the cell grows
                        rg = mu*(bac_m(ex(e)));
                        NRVxy = MatrixMet*rg; % rates of all soluble species
                        
                        %                     elseif mu = 0
                        %                         rg = 0;
                        %                         NRVxy = MatrixCat(:,k)*bac_mu_max(ex(e))*bac_m(ex(e));
                    else
                        rg = mu*(bac_m(ex(e))); % the cell shrinks
                        NRVxy = MatrixDecay(:,k)*(-rg); % rates of all soluble species
                    end
                    NRVbac(ex(e)) = rg; % rates for bacteria
                    NRVxy(1:numStVLiq2) = NRVxy(1:numStVLiq2)/Vg; % rates for liquid components - divided by volume - concentrations
                    NRV(indLiq) = NRV(indLiq)+NRVxy; % all rates in all grid cells

                end
            end
        end
    end
end
NRV = NRV + Ngas; % update the gas reaction rates

X = X*R.bac.bac_MW; % compute biomass mass concentration
if R.flagGas == 2
    nrv_gas = sum(-(Ngas(1:numStVLiq2*nT)< 0).*Ngas(1:numStVLiq2*nT)*Vg);
    R.Qgas = nrv_gas*R.pOp.Rg*R.pOp.T/R.pOp.P; %L
end
% compute the corr ections for diffusion coefficients depending on biomass
% concentration in each grid cell
Df = 1 - (0.43.*X.^0.92)./(11.19 + 0.27.*X.^0.99);
% Df = ones(size(X));
end


function [spcM, Sh] = solve_pH(Sh_ini, StV, Keq, chrM, flagpH)
w = 1;
% the pH solver
if flagpH == 1
    % if pH is assumed constant - spcM - computed as before - all that
    % changes is the initial concentrations
    
    spcM = zeros(size(chrM));
    Sh = Sh_ini;
    Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);
    
    spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;
    spcM(:,2) = (StV * Sh^3)                                    ./Denm;
    spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
    spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
    spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;
    
    
else
    % if pH is allowed to vary - must solve the charge balance
    % Checking the existence of a zero pool in the function between pH 1 and 14
    a=1e-14;
    b=1;
    Sh_v = [a; b]; F = zeros(2,1);
    for nn = 1:2
        Sh=Sh_v(nn);
        spcM = zeros(size(chrM));
        Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);
        
        spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;
        spcM(:,2) = (StV * Sh^3)                                    ./Denm;
        spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
        spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
        spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;
        
        % Evaluation of the charge balance for the current Sh value, F(Sh)
        F(nn) = Sh + sum(sum(spcM.*chrM)); % compute the charge balance at pH = 1 and pH = 14
    end
    FF = prod(F);
    if FF > 0 || isnan(FF)
        error('ERROR.- The sum of charges returns a wrong value')
    end
    % if the product > 0 - they have the same sign - so no 0 can exist for
    % the charge balance to close
    
    fa = F(1);
    fb = F(2);
    
    % Newton-Raphson method.-
    Sh = Sh_ini;
    % Counter of convergences
    ipH=1;    Tol = 5.e-15;    maxIter = 20;
    % Inicialization of matrix of species
    spcM = zeros(size(chrM)); dspcM = spcM;
    while ipH <= maxIter
        Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);
        spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;
        spcM(:,2) = (StV .* Sh^3)                                   ./Denm;
        spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
        spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
        spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;
        
        % Evaluation of the charge balance for the current Sh value, F(Sh)
        F = Sh + sum(sum(spcM.*chrM));
        
        % Calculation of all derivated functions
        dDenm = Denm.^2;
        aux = 3*Sh^2*(Keq(:,1)/w + 1) + 2*Sh*Keq(:,2) + Keq(:,2).*Keq(:,3);
        
        dspcM(:,1) =  (3*Sh^2*Keq(:,1).*StV)./(w*Denm) - ((Keq(:,1).*StV*Sh^3).*aux) ./(w*dDenm);
        dspcM(:,2) =  (3*Sh^2*StV)./Denm - (StV*Sh^3.*aux) ./dDenm;
        dspcM(:,3) = (2*Sh*Keq(:,2).*StV)./Denm - ((Keq(:,2).*StV*Sh^2).*aux) ./dDenm;
        dspcM(:,4) = (Keq(:,2).*Keq(:,3).*StV)./Denm - ((Keq(:,2).*Keq(:,3).*StV*Sh).*aux)./dDenm;
        dspcM(:,5) = -(Keq(:,2).*Keq(:,3).*Keq(:,4).*StV.*aux) ./dDenm;
        
        % Evaluation of the charge balance for the current Sh value, dF(Sh)
        dF = 1 + sum(sum(dspcM.*chrM));
        %Error
        err = F/dF;
        % Newton-Raphson algorithm
        Sh = Sh - err;
        
        if (abs(err) < 1e-14) && (abs(F) < Tol)
            % Checking if a valid pH was obtained
            if (Sh > 1e-14) && (Sh < 1)
                ipH = maxIter;
            else
                % Counter of convergence
                ipH = 1; maxIter = 50;
                n1 = 0; n2 = 0;
                while (ipH < maxIter)
                    Sh = (fb*a-fa*b)/(fb-fa);
                    spcM = zeros(size(chrM));
                    Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);
                    
                    spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;
                    spcM(:,2) = (StV .* Sh^3)                                   ./Denm;
                    spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
                    spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
                    spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;
                    
                    fc = Sh + sum(sum(spcM.*chrM));
                    if fa*fc > 0
                        n1 = n1+1;
                        if n1 == 2
                            fb = (fc/(fc+fa))*fb;
                            n1 = 0;
                        end
                        a = Sh; fa = fc;
                    elseif fb*fc > 0 % To avoid problems when fc == 0
                        n2 = n2+1;
                        if n2 == 2
                            fa = (fc/(fc+fb))*fa;
                            n2 = 0;
                        end
                        b = Sh; fb = fc;
                    end
                    err1 = abs(fc);
                    err2 = abs(Sh-(fb*a-fa*b)/(fb-fa));
                    if (err1 < Tol) && (err2 < 1e-14)
                        ipH = maxIter;
                    end
                    ipH = ipH+1;
                end
            end
        end
        ipH = ipH+1;
    end
end
end