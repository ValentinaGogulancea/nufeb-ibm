function [NRV,Df,NRVbac] = react_term(AllConc,R)

% take stuff out of the structure that we need

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
SLiqxy_all = AllConc; % concentrations of all solubles in all grids
RthT = R.pOp.Rth*R.pOp.T; % the product of R and temperature - in kJ/mole
DGr0 = R.pTh.DGr0; % free energy of the cat+ ana
DGdis = R.pTh.DGdis; % dissipation energy
spcR = R.pTh.spcR; % matrix that specifies the types of charges the species posses in the reactions written
lCat_eD = R.pTh.lCat_eD; % stoichiometric coeff of eDonors

numStVLiq = R.St.numStVLiq; % all soluble
numStVLiq2 = R.St.numStVLiq2; % all liquid
% numStVGas2Liq = R.St.numStVGas - 1; % all gas
numX = R.St.numX; % bacteria
numStFull = R.St.numStFull; % all components
% Liq2Gas = R.St.Liq2Gas; % those from gas to liquid

% Glu_pos = strcmp(R.St.StNames(1:numStVLiq), 'Glu');
% AcH_pos = strcmp(R.St.StNames(1:numStVLiq), 'AcH');
% H2_pos = strcmp(R.St.StNames(1:numStVLiq), 'H2');
O2_pos = strcmp(R.St.StNames(1:numStVLiq), 'O2');
NH3_pos = strcmp(R.St.StNames(1:numStVLiq), 'NH3');
NO2_pos = strcmp(R.St.StNames(1:numStVLiq), 'NO2');
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
                    % compute the rates for each bacteria now
                    if ne(bac_Ks(ex(e),1),0)
                        M = M*spcM(NH3_pos,rv(NH3_pos))/(spcM(NH3_pos,rv(NH3_pos)) + bac_Ks(ex(e),1));
%                         use the proper concentration, corresponding to the proper state
                    end
                    if ne(bac_Ks(ex(e),2),0)
                        M = M*(spcM(NO2_pos,rv(NO2_pos))/(spcM(NO2_pos,rv(NO2_pos)) + bac_Ks(ex(e),2)));
                    end
                    if ne(bac_Ks(ex(e),3),0)
                        M = M*(spcM(O2_pos,rv(O2_pos))/(spcM(O2_pos,rv(O2_pos)) + bac_Ks(ex(e),3)));
                    end
                    
%                     M = M*spcM(Glu_pos,rv(Glu_pos))/(spcM(Glu_pos,rv(Glu_pos)) + bac_Ks(ex(e),1));

%                     if ne(bac_Ks(ex(e),1),0)
%                         M = M*spcM(Glu_pos,rv(Glu_pos))/(spcM(Glu_pos,rv(Glu_pos)) + bac_Ks(ex(e),1));
%                         use the proper concentration, corresponding to the proper state
%                     end
%                     if ne(bac_Ks(ex(e),2),0)
%                         M = M*(spcM(H2_pos,rv(H2_pos))/(spcM(H2_pos,rv(H2_pos)) + bac_Ks(ex(e),2)));
%                     end
%                     if ne(bac_Ks(ex(e),3),0)
%                         M = M*(spcM(AcH_pos,rv(AcH_pos))/(spcM(AcH_pos,rv(AcH_pos)) + bac_Ks(ex(e),3)));
%                     end
                    mu = bac_mu_max(ex(e))*M - bac_maint(ex(e));
                    
                    if mu > 0  % the cell grows
                        rg = mu*(bac_m(ex(e)));
                        NRVxy = MatrixMet*rg; % rates of all soluble species
                    
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
% compute the corrections for diffusion coefficients depending on biomass
% concentration in each grid cell
% Df = 1 - (0.43.*X.^0.92)./(11.19 + 0.27.*X.^0.99);


 Df = 0.8 * ones(size(X));
end


