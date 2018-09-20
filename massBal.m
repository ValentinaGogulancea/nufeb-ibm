function R = massBal( L, G, X, R)
% updates the structure - computes the rates of growth of bacteria,
% consumption/formation of soluble species in each grid cell by calling
% mykinetics

St = R.St; % concentrations of all species - averaged
Sxy = R.Sxy; % values in each grid cell
numStVLiq2 = St.numStVLiq2; % num = number - only liquid
numStVLiq = St.numStVLiq; % liquid + gas
Sxy.StVLiq = [L; G]; % concentrations in each grid cell
Sxy.StVLiq2 = L; 
Sxy.StVGas = G; 
R.Sxy = Sxy;% just updated the concentration field
R.bac.bac_m = X; % and the bacteria - molar 
% NRVbac = rate of bacterial  particles growth - length == number of bacterial cells
% NRV = rates of soluble species consumption/formation; - length == number of grid cells* number of soluble species
% % R.Sxy.pH = pH values in each grid cell, same goes of DGcat & DGan
% Df - diffusion correction coefficients at each grid cell
% X - molar concentration of biomass (all bacteria)


[NRVbac, NRV, R.Sxy.pH, R.Sxy.DGcat, R.Sxy.DGan, Df,spcM] = my_kinetics(R);
rm = R.rm;
% update the difusion coefficients using the correction computed above
Diff = R.kTr.Diffwater.*kron(Df, ones(St.numStVLiq2,1));
% divide the diffusion coefficients by two - centered finite differences
% equation uses them as such
R.kTr.a = repmat(Sxy.dT*Diff/2, 1, Sxy.nT*numStVLiq2);
% update the values in the structure
R.kTr.Diff = Diff;
% NH3 = sum( U1(1:25))/25;
% NO2 = sum( U1(26:50))/25;
% O2 = sum( U1(51:75))/25;
% CO2 = sum( U1(76:100))/25;
% ammonium = spcM(1,2); % NH3
% 
% ammonia = spcM(1,3); % NH4+
% hno2 = spcM(2,2);
% no2 = spcM(2,3);
% co2 = spcM(4,2);
% carbonate = spcM(4,3);

if R.flagGas == 1
    %%%%% Gas correction
    ind = [Sxy.nT*((numStVLiq2+1:numStVLiq)'-1)+1,Sxy.nT*(numStVLiq2+1:numStVLiq)'];
    for k = 1:(numStVLiq - numStVLiq2)
        if R.flagGas == 2
            nrv_gas = sum(NRV(ind(k,1):ind(k,2)))/Sxy.nT - Sxy.StVGas(1+(k-1)*Sxy.nT:k*Sxy.nT)*R.Qgas/R.pOp.Vgas;
            NRV(ind(k,1):ind(k,2)) = nrv_gas;
        else
            NRV(ind(k,1):ind(k,2)) = 0;
        end
    end
end
aux = numStVLiq2*Sxy.nT;
% updates the reaction part of the structure
rm.NRVgas = NRV(aux+1:end);
rm.NRVliq = NRV(1:aux);
rm.NRVbac = NRVbac; % rates for bacteria 
R.rm = rm;
end