% %%%% bacteria.m -Generates the IbM: Position of the cells and division
% %%%% R.bac is the structure that keeps all the information related with this
% subscript

function [R, Fdiv] = bacteria1(L, G, X, R)
St = R.St;
Sxy = R.Sxy;
dz = R.Sxy.dz;
bac = R.bac;
numStVLiq = St.numStVLiq;
numStVLiq2 = St.numStVLiq2;

U = zeros(St.numSt,1);
for k = 1: St.numSt
    if k <= numStVLiq2 % values for the liquid species
        U(k) = Sxy.Sbc_Dir(k);
    elseif k <= numStVLiq % values for the gas species
        U(k) = G(Sxy.nT*((k-numStVLiq2)-1)+1);
    elseif k > numStVLiq % bacteria concentrations
        U(k) = sum(X.*(bac.bac_s==(k-numStVLiq)))/R.pOp.Vr;
    end
end
% %%%% Update of the variables
St.StV = U; % update the concentrations on the boundary
St.StVLiq = St.StV(1:numStVLiq);% Liquid and gas variables
St.StVLiq2 = St.StV(1:numStVLiq2);% Liquid variables
St.StVX = St.StV(numStVLiq+1:end);% Biomass
Sxy.StVLiq = [L; G];
Sxy.StVLiq2 = L;
Sxy.StVGas = G;
R.bac.bac_m = X;
R.Sxy = Sxy;
R.St = St;

[R.bac, R.EPS, Fdiv] = bac_division(R.pTh,R.bac,R.EPS, R.Sxy.maxxSys,R.Sxy.dz);
end

function [bac, EPS,Fdiv] = bac_division(pTh,bac, EPS,maxxSys,dz)
% number of bacteria + EPS
bac_n = bac.bac_n;

Fdiv = 0;
bac_m = bac.bac_m * bac.bac_MW; % Using grams from here;


  
j = 0;
for i = 1:bac.bac_n
 
        bac.bac_r(i) = sqrt((bac_m(i)/bac.bac_rho)/(dz*pi)); % radius of bacteria
end
% check for division only on bacteria - no EPS
% aux = bac_m(~EPS.type); % just the bacterial mass
% sum(aux> bac.bac_mmax) >= 1 - if they are big enough to divide
% (sum(bac_m < bac.bac_mmin) >= 1) - if they are small enough to be taken out
% sum(bac_m < bac_bEPS * bac.bac_rho/bac.EPS_rho - if the EPS volume surpases the one of the cell it is attached to
iffdiv = (sum(bac_m> bac.bac_mmax) >= 1)+(sum(bac_m < bac.bac_mmin) >= 1);
if ne(iffdiv,0) % if any of the above take place
    nold = bac_n; % use this as a comparison later on
    Fdiv = 1;
    bac_ns = bac.bac_ns;
    bac_x = bac.bac_x;
    bac_y = bac.bac_y;
    bac_s = bac.bac_s;
    bac_r = bac.bac_r;
    bac_ra = bac.bac_ra;
    bac_mu_max = bac.bac_mu_max;
    bac_yield = bac.bac_yield;
    bac_Ks = bac.bac_Ks;
    bac_maint = bac.bac_maint;
    bac_decay = bac.bac_decay; % 
    
    while sum(bac_m > bac.bac_mmax) >= 1 || sum(bac_m < bac.bac_mmin) >= 1 
        for i=1:length(bac_m)
            % first check for EPS secretion
            
            % second check the division
            if bac_m(i-j) > bac.bac_mmax  % check it is not EPS
                % choose an angle for the division
                fi = rand*pi;
                % increase number of cells with one
                bac_n = bac_n + 1; % a new cell is formed
                bac_ns(bac_s(i-j)) = bac_ns(bac_s(i-j)) + 1; % the cell type of the original cell gets a new addition
                bac_x(bac_n,1)  = bac_x(i-j) + 2*bac_r(i-j)*cos(fi); % 0x coordinates for the daughter
                bac_y(bac_n,1) = bac_y(i-j) + 2*bac_r(i-j)*sin(fi); % 0y coordinates for the daughter
                if bac_y(bac_n,1) < bac.bac_rmax
                    bac_y(bac_n,1) = bac.bac_rmax;
                end
                bac_s(bac_n,1)  = bac_s(i-j);                    % same cell type
                % same kinetic parameters
                bac_mu_max(bac_n,1)  = bac_mu_max(i-j);
                bac_yield(bac_n,1)  = bac_yield(i-j);
                bac_Ks(bac_n,:)  = bac_Ks(i-j, :);
                bac_maint(bac_n,1)  = bac_maint(i-j);
                bac_decay(bac_n,1) =  bac_decay(i-j);
                %
                bac_m(bac_n,1)  = (0.45 + 0.1*rand(1,1))*bac_m(i-j); % mass of the new cell
                %
                bac_bEPS(bac_n,1) = bac_m(bac_n,1) / bac_m(i-j) * bac_bEPS(i-j); % EPS content of the new cell
                %
                bac_m(i-j) = bac_m(i-j) - bac_m(bac_n) ;        % remaining mass with the parent
                %
                bac_bEPS(i-j) = bac_bEPS(i-j) - bac_bEPS(bac_n) ;        % remaining EPS mass with the parent
                bac.producers(bac_n,1) = bac.producers(i-j,1);
               
                EPS.type(bac_n,1) = 0;
                
                if ~bac.producers(i-j) % if the cells do not produce EPS go as usual
                    bac_r(bac_n) = sqrt((bac_m(bac_n)/bac.bac_rho)/(dz*pi)); % radius of new cell
                    
                    bac_ra(bac_n) =  bac_r(bac_n);
                    bac_r(i-j) = sqrt((bac_m(i-j)/bac.bac_rho)/(dz*pi));% new radius of the parent
                    bac_ra(i-j) =  bac_r(i-j);
                else % if cells produce EPS
                    % radius of new cell + EPS
                    bac_r(bac_n) = sqrt((bac_m(bac_n)/bac.bac_rho + bac_bEPS(bac_n)/bac.EPS_rho)/(dz*pi));
                    bac_ra(bac_n) = sqrt(bac_m(bac_n)/bac.bac_rho /(dz*pi)); % active biomass radius
                    % radius of old cell + EPS
                    bac_ra(i-j) = sqrt((bac_m(i-j)/bac.bac_rho )/(dz*pi)); % active biomass
                    bac_r(i-j) = sqrt((bac_bEPS(i-j)/bac.EPS_rho/dz+ pi*bac_ra(i-j)^2)/pi);% new radius of the parent
                    
                end
                
            end
            if bac_m(i-j) < bac.bac_mmin
                % decrease number of cells with one
                bac_n = bac_n - 1;
                bac_ns(bac_s(i-j)) = bac_ns(bac_s(i-j)) - 1;
                % delete of coordinates, species and radius
                bac_m(i-j)  = [];
                bac_r(i-j)  = [];
                bac_ra(i-j)  = [];
                bac_x(i-j)  = [];
                bac_y(i-j)  = [];
                bac_s(i-j)  = [];
                bac_mu_max(i-j) = [];
                bac_yield(i-j) = [];
                bac_Ks(i-j, :) = [];
                bac_maint(i-j) = [];
                bac.producers(i-j) = [];
                EPS.type(i-j) = [];
                j = j + 1;
            end
            
        end
    end
end

if Fdiv == 1
    bac.bac_n = bac_n;
    bac.bac_ns = bac_ns;
    bac.bac_s = bac_s;
    for i = 1:bac_n
        if ~EPS.type
            bac.bac_m(i,1) = bac_m(i)/bac.bac_MW; %%Again converting to mol
        else
            bac.bac_m(i,1) = bac_m(i)/bac.EPS_MW;
        end
    end
    % and update the structure to account for the new addition
    bac.EPS_m = bac_bEPS/bac.EPS_MW; % moles
    bac.bac_bEPS = bac_bEPS; % grams
    bac.bac_r = bac_r;
    bac.bac_ra = bac_ra;
    bac.bac_mu_max = bac_mu_max;
    bac.bac_yield = bac_yield;
    bac.bac_Ks = bac_Ks;
    bac.bac_maint = bac_maint;
    bac.bac_decay = bac_decay;
    % actual shoving algorithm
    [bac_x, bac_y] = bac_shovingloops_mex(bac_x, bac_y, bac_r', bac_m, bac.k, bac.s_dist, bac.overlap, bac.bac_rmax, maxxSys);
    % new coordinates of bacteria
    bac.bac_x = bac_x;
    bac.bac_y = bac_y;
end
end