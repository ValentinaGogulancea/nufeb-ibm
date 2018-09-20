% %%%% loadModelXls.m -Load all the parameters from Excel sheet and
% generates the initial structure

function R = loadModelXlsx

R = []; % R initialization

route = char('AOBNOB-run1.xlsx');% % specifying the route to the excel file for xlsread function

% route = char('AOBNOB-run1.xlsx');
% GENERAL MODEL PARAMETERS from the Excel file
 fprintf('> LOADING AND CREATING MODEL STRUCTURE AND PARAMETERS...\n')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General parameters of the simulation (Temperature, pH, ...)
% Operating parameters of the simulation (Temperature, pH, ...)
% inlet concentrations of C & N - source
% Temperature - K 
% Values for ideal gas constant 
% - Rth - in kJ/moleK for free energy computations
% - Rgas - atmL/molK for the gas liquid mass transfer
% Pressure - bar
% Volume for gas  Vgas == V computational domain m^3
% Volume for reactor Vr == V computational domain m^3
% pH - initial value
% HRT - hydraulic residence time - h
% Qliq - volume flow rate - L/h
[OperatParam, pOpNames]   = xlsread(route, strcat('Parameters'));
aux = '';
for i=1:length(OperatParam)
    % Building the string command for the structure, last without ', '
    if i<length(OperatParam),   
        aux = strcat(aux, char(39), pOpNames(i), char(39), ', ',num2str(OperatParam(i)),', '); 
    else
        aux = strcat(aux, char(39), pOpNames(i), char(39), ', ',num2str(OperatParam(i))); 
    end
end
eval(strcat('pOp = struct(', char(aux), ');'));
pOp.Vgas = pOp.Vgas*1000; %#ok<NODEF> %L
pOp.Vr = pOp.Vr*1000; %%L - used for computing concentrations
pOp.invHRT = 1/pOp.HRT;% used in the reactor mass balance - bulk liquid,h^-1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE VARIABLES Structure and initial values

% initial concentrations and the aggregation states of the main components
% of the system
% names that are to be used further
% L - solubles, G - gas components, S - bacteria & EPS

[States, StNames] = xlsread(route, strcat('States'));
St.StV = States(:,1); % concentrations of all - if initial concentration is zero - use 1e-20 to avoid dividing by 0

St.Phase = StNames(:,end); % the state of the components

St.StNames = StNames(:,1); % the name of the component

St.StVLiq = St.StV(1:(find(strcmp(St.Phase,'S'),1)-1)); % soluble species, i.e. not particles

St.StVLiq2 = St.StV(1:(find(strcmp(St.Phase,'G'),1)-1)); % liquid species - i.e. soluble - gas species

St.StVIni = St.StVLiq; % initial concentrations for gas + liquid

St.StVX = St.StV(find(strcmp(St.Phase,'S'),1):end); % initial bacterial concentration

St.numX = length(St.StVX);% the number of bacterial types

St.numStVLiq = length(St.StVLiq); % number of soluble species considered
St.numStVLiq2 = length(St.StVLiq2); % number of liquid species considered
St.numStVGas = St.numStVLiq - St.numStVLiq2;% number of gaseous species considered
St.numSt = length(St.StNames); % number of total species considered: gas, liquid, bacteria

name = char(St.StNames(St.numStVLiq2 + 1:St.numStVLiq)); % auxiliarry variable storing the names of gas components as char
name = name(:,2:end);% slices the first letter (g) of the name of gas components

% this is to see which of the model species can be transfered to or from
% the gas phase

Liq2Gas = zeros(St.numSt, 1);
for k = 1 : St.numStVGas
    Liq2Gas(strcmp(strcat(name(k, :)),St.StNames)) = k; % check if the gas species is found as a liquid species as well
end
St.Liq2Gas = Liq2Gas; % vector that assigns the indices to the liquid species corresponding to the position of the gas species

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STOICHIOMETRY/REACTION MATRIX 

[rMatrixFull, rNames] = xlsread(route, strcat('ReactionMatrix'));
rm.rMatrixFull = rMatrixFull(:,[1:3:end;2:3:end]); % the entire excel writing - aka
% all species of the model + the stoichiometric coefficients for all the metabolisms encountered
rm.rNamesX = rNames(1,2:3:end)';  % names of the bacteria/eps for which the metabolism is provided
rm.rNamesS = rNames(3:end,1); % names of the components of the model
MatrixCat = rMatrixFull(:,1:3:end); % catabolic matrix in full for all 
MatrixAn = rMatrixFull(:,2:3:end); % anabolic matrix in full
MatrixDecay = rMatrixFull(:,3:3:end); % decay matrix in full 
St.numStFull = length(MatrixCat); % all the species in the model

MatrixCat_mod = [MatrixCat(1:(find(strcmp(rm.rNamesS, 'H2O'))-1),:); MatrixCat((find(strcmp(rm.rNamesS, 'H2O'))+1):end,:)]; 
namesAux = [rm.rNamesS(1:(find(strcmp(rm.rNamesS, 'H2O'))-1),:); rm.rNamesS((find(strcmp(rm.rNamesS, 'H2O'))+1):end,:)]; % removing water
MatrixCat_mod = [MatrixCat_mod(1:(find(strcmp(namesAux, 'H'))-1),:);MatrixCat_mod((find(strcmp(namesAux, 'H'))+1):end,:)]; %removing  hydrogen
MatrixCat_mod = MatrixCat_mod(1:St.numStVLiq,:); %removing bacteria but keeping gas components
MatrixAn_mod = [MatrixAn(1:(find(strcmp(rm.rNamesS, 'H2O'))-1),:); MatrixAn((find(strcmp(rm.rNamesS, 'H2O'))+1):end,:)];% removing water
MatrixAn_mod = [MatrixAn_mod(1:(find(strcmp(namesAux, 'H'))-1),:);MatrixAn_mod((find(strcmp(namesAux, 'H'))+1):end,:)];%removing  hydrogen
MatrixAn_mod = MatrixAn_mod(1:St.numStVLiq,:);%removing bacteria but keeping gas components
MatrixDecay_mod = [MatrixDecay(1:(find(strcmp(rm.rNamesS, 'H2O'))-1),:); MatrixDecay((find(strcmp(rm.rNamesS, 'H2O'))+1):end,:)];% removing water
MatrixDecay_mod = [MatrixDecay_mod(1:(find(strcmp(namesAux, 'H'))-1),:);MatrixDecay_mod((find(strcmp(namesAux, 'H'))+1):end,:)];%removing  hydrogen
MatrixDecay_mod = MatrixDecay_mod(1:St.numStVLiq,:);%removing bacteria but keeping gas components
% and transfering everything to the structure
rm.MatrixCat = MatrixCat; 
rm.MatrixAn = MatrixAn;
rm.MatrixDecay = MatrixDecay;
rm.MatrixCat_mod = MatrixCat_mod;
rm.MatrixAn_mod = MatrixAn_mod;
rm.MatrixDecay_mod = MatrixDecay_mod;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure with the THERMODYNAMIC PARAMETERS

% it holds the values of standard free energies of formation for each state of each component
% it also holds the charges that can be encountered

[pThValues, pThNames]   = xlsread(route, strcat('ThermoParam')); 
pThValues = pThValues(2:end,:); % all values except charges of deltag's and charge matrix
pThNames = pThNames(2:end,:); % all names of the species x 2 - 
DG0 = pThValues(1:(size(pThValues,1)/2), 1:(end-1)); % the values of the standard free energies of formation
nan_locations = isnan(DG0);% where a species is not found in the corresponding charged state and thus has no free energy
DG0(nan_locations) = 1e10; % the nan value is replaced with one large eneough to cancel out in the energy computations
chrM = pThValues((size(pThValues,1)/2 +1):size(pThValues,1), 1:(end-1)); % charge matrix
nan_locations = isnan(chrM); %  where a species is not found in the corresponding charged state
chrM(nan_locations) = 0; % the nan value is replaced with 0 to cancel out in the pH computations

sDG0col = size(DG0,2); 

Keq = zeros(size(DG0,1),sDG0col-2); % initialization for computation of equilibrium constants for dissociations
DG0H2O = -237.18; % Water free energy - used for the hydration equilibrium computations kJ/mol

% the first column specifies if the species undergo hydration before the
% dissociate 
Keq(:,1) = exp((DG0H2O + DG0(:,1) - DG0(:,2))/(-pOp.Rth*pOp.T)); 

for i=2:(sDG0col - 1)
    Keq(:,i) = exp((DG0(:,i+1)+(i*1e10*(DG0(:,i+1) == 1e10))-(DG0(:,i)))/(-pOp.Rth*pOp.T)); % dissociations
end
v1 = [(1:(St.numStFull))', pThValues(1:length(pThNames)/2,end)]; % the components and the number of values to be specified for charge computations    
 
v2 = [St.numStFull, sDG0col];
 
spcR = accumarray(v1, 1 , v2); % matrix that specifies which of the charged states is encountered in the systems
% this shows the charge in the system

pTh.Phase = pThNames(1:length(pThNames)/2,end);  % state of aggregation
pTh.chrM = chrM(1:(find(strcmp(pTh.Phase,'S'),1)-1), :);% charge matrix is restricted to the liquid and gas phase (obvs)
pTh.Keq = Keq(1:(find(strcmp(pTh.Phase,'S'),1)-1), :);   % dissociation constants too
pTh.DGr0 = rm.rMatrixFull'*sum(spcR.*DG0,2);  % computing the standard free energy for the anabolism and catabolism 
pTh.DG0 = DG0; % retaining the free energies of formation

spcR = spcR(1:(find(strcmp(pTh.Phase,'S'),1)-1), :); % keeping only the gas and liquid species

pTh.spcR = spcR == 1; % the same as above but of logical type

pTh.react_v = pThValues(1:length(pThNames)/2,end);  % indicates how many charged species can appear
vaux = find(strcmp(St.Phase,'G'));
vaux = vaux(ne(vaux, find(strcmp(St.StNames,'gN2')))); % no nitrogen mass transfer 

KhV = zeros(length(St.StV),1); % initialization for the values for Henry's constants
for i=1:length(vaux)
    name = char(St.StNames(vaux(i)));
    name = name(2:end);
    w = strcmp(strcat(name),St.StNames);
    if DG0(w,1)== 1e10
        % computing the Henry's constant values - if the species doesn't
        % exist in the charged/hydrated form

        Kh.(name) = exp((DG0(w,2) - DG0(vaux(i),2))  / (-pOp.Rth*pOp.T));% M/atm
    else
        Kh.(name) = exp((DG0(w,1) - DG0(vaux(i),2))  / (-pOp.Rth*pOp.T));% M/atm
    end
    KhV(i) = Kh.(name);   
end

pTh.Kh = Kh; % 0 if the species does not exist in the gas phase - computed otherwise
pTh.Kh.KhV = KhV;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the microbial species (qmax, Ks, Yield)
[SpecParam, SpecNames] = xlsread(route, strcat('SpecParam'));
% pretty self explanatory here
aux = strcmp('mu_max',SpecNames(1,2:end));
pTh.mu_max = SpecParam((isfinite(SpecParam(:,aux))),aux);

aux = strcmp('Ks_O2',SpecNames(1,2:end));
Ks_O2 = SpecParam((isfinite(SpecParam(:,aux))),aux);

aux = strcmp('Ks_NH3',SpecNames(1,2:end));
Ks_NH3 = SpecParam((isfinite(SpecParam(:,aux))),aux);

aux = strcmp('Ks_NO2',SpecNames(1,2:end));
Ks_NO2 = SpecParam((isfinite(SpecParam(:,aux))),aux);

aux = strcmp('Ks_Glu',SpecNames(1,2:end));
Ks_Glu = SpecParam((isfinite(SpecParam(:,aux))),aux);

aux = strcmp('Ks_AcH',SpecNames(1,2:end));
Ks_AcH= SpecParam((isfinite(SpecParam(:,aux))),aux);

aux = strcmp('Ks_H2',SpecNames(1,2:end));
Ks_H2 = SpecParam((isfinite(SpecParam(:,aux))),aux);

pTh.Ks = [Ks_NH3,Ks_NO2,Ks_O2];
% pTh.Ks = [Ks_Glu,Ks_AcH,Ks_H2];
aux = strcmp('Yield',SpecNames(1,2:end));
pTh.Y = SpecParam((isfinite(SpecParam(:,aux))),aux);

aux = strcmp('DGdis',SpecNames(1,2:end));
pTh.DGdis = SpecParam((isfinite(SpecParam(:,aux))),aux);

aux = strcmp('Maintenance',SpecNames(1,2:end));
pTh.maint = SpecParam((isfinite(SpecParam(:,aux))),aux);


eD = SpecNames(2:end,end); % specifies the electron donor
lCat_eD = zeros(St.numX,1);
for i=1:St.numX
   lCat_eD(i) = -rm.MatrixAn_mod(strcmp(rm.rNamesS, eD(i,:)), i);
end
pTh.lCat_eD = lCat_eD;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the IbModel (maximum cell's weight, parameters for shoving...)
% bac_mmax = maximum bacteria mass kg - thereshold for division; 
% bac_mmin = minimum bacteria mass kg - thereshold for being taken out of the system
% bac_rmax = maximum bacteria radius m - alternative thereshold for division
% bac_rho = density of bacteria, kg/m^3
% bac_MW = molecular weight of bacteria, g/mol
% EPS_MW = molecular weight of EPS, g/mol
% frac_EPS = mass ratio of bound EPS to bacteria initially, g/g 
% bac_nmax = maximum number of bacteria that can exist in the domain
% k - multiplied with the radius - distance between the centres of bacterial cells - for initial positioning, 
% overlap - allowed overlap distance
% s - diameter of bacteria

[BacteriaParam, BacteriaNames]   = xlsread(route, strcat('Bacteria'));
aux = '';
for i=1:length(BacteriaParam)
    % Building the string command for the structure, last without ', '
    if i<length(BacteriaParam),   
        aux = strcat(aux, char(39), BacteriaNames(i), char(39), ', ',num2str(BacteriaParam(i)),', '); 
    else
        aux = strcat(aux, char(39), BacteriaNames(i), char(39), ', ',num2str(BacteriaParam(i))); 
    end
end
eval(strcat('bac = struct(', char(aux), ');'));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Space discretization
% nx - number of divisions on the Ox
% ny - number of divisions on the Oy
% maxx & maxy - the limits of the computational domain on Ox & Oy, m
% dx,dy - the length & height  of the discretization cells, m
% dz - the depth of the Oz axis - considered for cell volume computations,
% T_blayer - height of boundary layer, m
% maxT - computational time - hours!
% dT - time step for difusion
% dT_bac - time step for bacteria mass balance
% dT_print - time step for printing
% Tol - tolerance for steady state approximation in the difusion computations

[DiscretizationParam, DiscretizationNames] = xlsread(route, strcat('Discretization'));

aux = '';
for i=1:length(DiscretizationParam)
    % Building the string command for the structure, last without ', '
    if i<length(DiscretizationParam),   
        aux = strcat(aux, char(39), DiscretizationNames(i), char(39), ', ',num2str(DiscretizationParam(i)),', '); 
    else
        aux = strcat(aux, char(39), DiscretizationNames(i), char(39), ', ',num2str(DiscretizationParam(i))); 
    end
end
eval(strcat('Sxy = struct(', char(aux), ');'));
% the actual number of coordinates required for ox & oy
Sxy.nxSys = Sxy.nx +2; %#ok<NODEF> the number of x coordinates
Sxy.nySys = Sxy.ny +2; % the number of y coordinates
Sxy.nTSys=(Sxy.nxSys-2)*(Sxy.nySys-2); % the number of cells for the domain
Sxy.maxxSys = Sxy.maxx;% maximum length of the domain in the x direction
Sxy.maxySys = Sxy.maxy;% maximum length of the domain in the y direction
Sxy.xSys = 0:Sxy.dx:Sxy.maxxSys; % x-coordinates
Sxy.ySys = 0:Sxy.dy:Sxy.maxySys; % y-coordinates 
xp = Sxy.xSys+ Sxy.dx/2; % x-coordinates of the center of grid cells
Sxy.x_pnSys = xp(1:end-1); % the last one is outside the domain - so it goes out
yp = Sxy.ySys+ Sxy.dy/2; % y-coordinates of the center of grid cells
Sxy.y_pnSys = yp(1:end-1); % the last one is outside the domain - so it goes out
Sxy.Vg = (Sxy.dx*Sxy.dy*Sxy.dz)*1000; %L - volume of one grid cell (assuming constant z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Influent - reactor balance --

[Influent, BLayer] = xlsread(route, strcat('Influent'));
R.Inf.St = Influent(:,1); % inlet concentrations for liquid species 

% the type of conditions - if N is specified - the boundary conditions are recalculated using 
% the reactor mass balance implemented in the  function 'boundary'
BLayer = BLayer(:,end); 
% specifies which conditions are Dirichlet - constant concentrations - Oxygen + CO2 in our case
Dir_k = zeros(St.numStVLiq2,1);

for k = 1: St.numStVLiq2
    if strcmp(BLayer(k),'D')
        Dir_k(k) = 1;
    else
        Dir_k(k) = 0;
    end
end
R.Inf.Dir_k = Dir_k;
Sxy.Sbc_Dir = St.StVIni; % initial values for the boundary conditions

Sxy.pHnT = pOp.pH*ones(Sxy.nTSys,1); % values of pH in every grid cell

% % Initial positioning of bacterial cells

bac_x = (0+0*bac.k*bac.bac_rmax:2*bac.k*bac.bac_rmax:Sxy.maxxSys-0*bac.k*bac.bac_rmax)'; %#ok<NODEF>
bac_x = bac_x(1:(floor(length(bac_x)/St.numX)*St.numX));
aux = (Sxy.maxxSys - bac_x(end))/2;
bac.bac_x = bac_x + aux;
auxSize = size(bac.bac_x);
bac.bac_y = bac.bac_rmax*ones(auxSize);

bac_mm = 0.9*bac.bac_mmax*ones(auxSize);
bac.bac_m = bac_mm/bac.bac_MW; %mol
St.StVX = bac.bac_m; 
bac.bac_r = ((bac_mm/bac.bac_rho)*(3/(4*pi))).^(1/3);%radius
bac.bac_n = auxSize(1); % number of cells
bac.bac_s = kron(ones(round(auxSize(1)/St.numX), 1), (1:St.numX)');        % species attribute (types of cells)
% makeshift = [2*ones(1,1); 3* ones(auxSize(1)-1,1)];
% bac.bac_s  = ones(auxSize(1),1);        % species attribute (types of cells)

bac.bac_mu_max = pTh.mu_max(bac.bac_s); 
bac.bac_yield = pTh.Y(bac.bac_s);
bac.bac_Ks = pTh.Ks(bac.bac_s, :);
bac.bac_maint = pTh.maint(bac.bac_s);

bac_ns = zeros(St.numX, 1);
for i  = 1:St.numX
    for j = 1:auxSize(1)
        bac_ns(i) =  (i == bac.bac_s(j)) + bac_ns(i); % number of cells of each
    end
end
bac.bac_ns = bac_ns; % the type of atoms/elements associated - each in its own cluster

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[DiffusionParam, DiffusionNames] = xlsread(route, strcat('Diffusion'));
% the diffusion coefficients for each liquid species (except water and the
% proton - obvs) - m^2/h
% the mass transfer coefficients for the liquid species - kLa - h^-1
 
kTr.Diffn = DiffusionParam(1:St.numStVLiq2,:); % values for liquid species
kTr.DiffNames = DiffusionNames(1:St.numStVLiq2,1); % names assigned
kTr.kLa = DiffusionParam(St.numStVLiq2+1:end,:); % values
kTr.kLaNames = DiffusionNames(St.numStVLiq2+1:end,1); % names
kTr.DiffwaternT = kron(kTr.Diffn, ones(Sxy.nTSys,1)); % the difusion coefficients for all the species in the every grid

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sxy.pos_xySys = 555;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pH calculations and concentration re-calculation
Keq = pTh.Keq; % values of eq. constants for deprotonations
chrM = pTh.chrM;    
spcR = pTh.spcR; w = 1; % activity coefficient
numStVLiq = St.numStVLiq;
spcM = zeros(size(chrM)); % matrix for the species
StV = [St.StVLiq;1;0];  % initial concentrations
Sh = 10.^(-pOp.pH); % proton concentration
Denm =(1+Keq(:,1)/w)*Sh^3 + Keq(:,2)*Sh^2 + Keq(:,3).*Keq(:,2)*Sh + Keq(:,4).*Keq(:,3).*Keq(:,2);
% computing the concentrations of each deprotonaitons 
spcM(:,1) = ((Keq(:,1)/w).*StV*Sh^3)                        ./Denm;        
spcM(:,2) = (StV * Sh^3)                                    ./Denm;
spcM(:,3) = (StV * Sh^2 .* Keq(:,2))                        ./Denm;
spcM(:,4) = (StV * Sh .* Keq(:,2) .* Keq(:,3))              ./Denm;
spcM(:,5) = (StV      .* Keq(:,2) .* Keq(:,3) .* Keq(:,4))  ./Denm;

spcMc = (spcM == 0)*1e-20 + spcM;
% concentrations of charged/uncharged species - where the corresponding charge does not appear
% use 1e-20 - to go 0 in the next line which represent
% computation of the free energy for catabolism and anabolism for each of
% the bacterial species
 % assigning the values   
DGr = pTh.DGr0 + pOp.Rth*pOp.T*(rm.rMatrixFull'*[nansum(spcR.*log(spcMc),2);ones(St.numStFull-(numStVLiq+2),1)]);
R.DGCat_bulk =  DGr(1:2:length(DGr));
R.DGAn_bulk =  DGr(2:2:length(DGr));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Update of the structure R
R.pOp = pOp;
R.pTh = pTh;
R.rm = rm;
R.St = St;
R.bac = bac;
R.Sxy = Sxy;
R.kTr = kTr;
assignin('base', 'R', R) % Writing R in the Workspace
% choosing which modules to use or not
% R.flagGas - gas-liquid mass transfer
% R.flagpH - pH constant or allowed to vary
% R.flagDG - yields constant or not
% R.flagN - fixed or variable ammonia

% prompt = '>> Do you want to keep the reactor open? (default [Y]) [Y = 1/N = 2]\n>> Answer: ';
% x = input(prompt);
% if isempty(x)
%     x = 1; 
% end
% if x == 1
%     fprintf('--> YES\n')
% else
%     fprintf('--> NO\n')
%     R.Qgas = 0;
%     x = 2;
% end
 R.flagGas = 1;
 R.Qgas = 0;
% 
% prompt = '>> Do you want to keep the pH fixed? (default [Y]) [Y = 1/N = 2]\n>> Answer: ';
% x = input(prompt);
% if isempty(x)
%     x = 1; 
% end
% if x == 1
%     fprintf('--> YES\n')
% else
%     fprintf('--> NO\n')
%     x = 2;
% end
 R.flagpH = 1;
% 
% prompt = '>> Do you want to keep the Yield fixed? (default [Y]) [Y = 1/N = 2]\n>> Answer: ';
% x = input(prompt);
% if isempty(x)
%     x = 1; 
% end
% if x == 1
%     fprintf('--> YES\n')
% else
%     fprintf('--> NO\n')
%     x = 2;
% end
 R.flagDG = 1;
% 
% 
% prompt = '>> Do you want to keep the NH3 fixed (variable HRT)? (default [Y]) [Y = 1/N = 2]\n>> Answer: ';
% x = input(prompt);
% if isempty(x)
%     x = 1; 
% end
% if x == 1
%     fprintf('--> YES\n')
%         if R.Inf.St(1) ==  R.St.StVIni(1)
%             error('NH3 influent = NH3 initial - Change the NH3 influent!')
%         end
% else
%     fprintf('--> NO\n')
%     x = 2;
%     if ne(R.Inf.St(1),  R.St.StVIni(1))
%             prompt = '>> NH3 influent is different than NH3 initial. Do you want to continue? (default [Y]) [Y = 1/N = 2]\n>> Answer: ';
%             y = input(prompt);
%             if isempty(y)
%                 y = 1; 
%             end
%             if y == 1
%                 fprintf('--> YES\n')
%             else
%                 fprintf('--> NO\n')
%                 error('NH3 influent is different than NH3 initial')
%             end
%     end
%     
% end
 R.flagN = 1;
% clear x prompt
% 
% fprintf('> ... MODEL STRUCTURE AND PARAMETERS LOADED.')
% fprintf('\n------------ooo\n')
% VN = Sxy.dT*kTr.Diffn/(Sxy.dx*Sxy.dy);
% for i =1:St.numStVLiq2
%     fprintf('Von Neumann coefficient of stability for %s: %.2e\n', char(St.StNames(i)),VN(i))
% end
% fprintf('>>> If Von Neumann coefficient larger than 0.5, consider Backward Euler instead of Crank-Nicolson')
% fprintf('\n------------ooo\n')
% for i = 1:(R.St.numX)
%     fprintf('Number of cells of type %d Name %s : %d\n', i,char(R.rm.rNamesX(i)),R.bac.bac_ns(i))
% end
% fprintf('>>> INITIAL TOTAL Number of cells: %d', R.bac.bac_n)
% fprintf('\n------------ooo\n')
% fprintf('\n>>> Push any key to start the simulation...\n')
% pause()