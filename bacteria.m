% %%%% bacteria.m -Generates the IbM: Position of the cells and division
% %%%% R.bac is the structure that keeps all the information related with this
% subscript

function [R, Fdiv] = bacteria(L, G, X, R)
St = R.St;
Sxy = R.Sxy;
bac = R.bac;
numStVLiq = St.numStVLiq;
numStVLiq2 = St.numStVLiq2;
dz = R.Sxy.dz;
U = zeros(St.numSt,1);
for k = 1: St.numSt
    if k <= numStVLiq2
        U(k) = Sxy.Sbc_Dir(k);
    elseif k <= numStVLiq
        U(k) = G(Sxy.nT*((k-numStVLiq2)-1)+1);
    elseif k > numStVLiq
        U(k) = sum(X.*(bac.bac_s==(k-numStVLiq)))/R.pOp.Vr;
    end
end
% %%%% Update of the variables
St.StV = U;
St.StVLiq = St.StV(1:numStVLiq);% Liquid and gas variables
St.StVLiq2 = St.StV(1:numStVLiq2);% Liquid variables
St.StVX = St.StV(numStVLiq+1:end);% Biomass
Sxy.StVLiq = [L; G];
Sxy.StVLiq2 = L;
Sxy.StVGas = G;
R.bac.bac_m = X;
R.Sxy = Sxy;
R.St = St;
[R.bac, Fdiv] = bac_division(R.bac, R.Sxy.maxxSys,R.Sxy.dz);
end

function [bac, Fdiv] = bac_division(bac, maxxSys,dz)
bac_n = bac.bac_n;

Fdiv = 0;
bac_m = bac.bac_m*bac.bac_MW;%Using grams from here

for i = 1:bac.bac_n
    
    bac.bac_r(i) = sqrt((bac_m(i)/bac.bac_rho)/(dz*pi)); % radius of bacteria
end
j = 0;
iffdiv = (sum(bac_m > bac.bac_mmax) >= 1)+(sum(bac_m < bac.bac_mmin) >= 1);
if ne(iffdiv,0)
    nold = bac_n;
    Fdiv = 1;
    bac_ns = bac.bac_ns;
    bac_x = bac.bac_x;
    bac_y = bac.bac_y;
    bac_s = bac.bac_s;
    bac_r = bac.bac_r;
    bac_rmax = bac.bac_rmax;
    bac_mu_max = bac.bac_mu_max;
    bac_yield = bac.bac_yield;
    bac_Ks = bac.bac_Ks;
    bac_maint = bac.bac_maint;
    while sum(bac_m > bac.bac_mmax) >= 1 || sum(bac_m < bac.bac_mmin) >= 1
        for i=1:nold
            if bac_m(i-j) > bac.bac_mmax
                % choose an angle for the division
                fi = rand*pi;
                % increase number of cells with one
                bac_n = bac_n + 1;
                bac_ns(bac_s(i-j)) = bac_ns(bac_s(i-j)) + 1;
                
                % coordinates, species and radius of the new cell (bac_n)
                bac_x(bac_n,1)  = bac_x(i-j) + 2*bac_r(i-j)*cos(fi); % 0x coordinates for the daughter
                bac_y(bac_n,1) = bac_y(i-j) + 2*bac_r(i-j)*sin(fi); % 0y coordinates for the daughter
                if bac_y(bac_n,1) < bac.bac_rmax
                    bac_y(bac_n,1) = bac.bac_rmax;
                end
                bac_s(bac_n,1)  = bac_s(i-j);                    % same cell type
                bac_mu_max(bac_n,1)  = bac_mu_max(i-j);                    % same mu max
                bac_yield(bac_n,1)  = bac_yield(i-j);
                bac_Ks(bac_n,:)  = bac_Ks(i-j, :);                    % same Ks
                bac_maint(bac_n,1)  = bac_maint(i-j);                    % same maintenance
                bac_m(bac_n,1)  = (0.45 + 0.1*rand(1,1))*bac_m(i-j);
                bac_r(bac_n,1)  = ((bac_m(bac_n)/bac.bac_rho)*(3/(4*pi)))^(1/3);
                % move old parent cell (ng)
                bac_m(i-j) = bac_m(i-j) - bac_m(bac_n) ;        % remaining mass with the parent
                bac_r(i-j) = ((bac_m(i-j)/bac.bac_rho)*(3/(4*pi)))^(1/3);% new radius of the parent
            end
            if bac_m(i-j) < bac.bac_mmin
                % decrease number of cells with one
                bac_n = bac_n - 1;
                bac_ns(bac_s(i-j)) = bac_ns(bac_s(i-j)) - 1;
                % delete of coordinates, species and radius
                bac_m(i-j)  = [];
                bac_r(i-j)  = [];
                bac_x(i-j)  = [];
                bac_y(i-j)  = [];
                bac_s(i-j)  = [];
                bac_mu_max(i-j) = [];
                bac_yield(i-j) = [];
                bac_Ks(i-j, :) = [];
                bac_maint(i-j) = [];
                j = j + 1;
            end
        end
    end
end

if Fdiv == 1
    bac.bac_n = bac_n;
    bac.bac_ns = bac_ns;
    bac.bac_s = bac_s;
    bac.bac_m = bac_m/bac.bac_MW; %%Again converting to mol
    bac.bac_r = bac_r;
    bac.bac_mu_max = bac_mu_max;
    bac.bac_yield = bac_yield;
    bac.bac_Ks = bac_Ks;
    bac.bac_maint = bac_maint;
    
    if length(bac_x)<3e4
        %      [bac_x, bac_y] = bac_shovingloops_mex(bac_x, bac_y, bac_r, bac_m, bac.k, bac.s_dist, bac.overlap, bac.bac_rmax, maxxSys);
        [bac_x, bac_y] = bac_shovingloops(bac_x, bac_y, bac_r, bac_m, bac.k, bac.s_dist, bac.overlap, bac.bac_rmax, maxxSys);
    else
        X = TestOverlap(bac_x, bac_y, bac_r, bac_m,bac.overlap, maxxSys,bac_rmax);
        
        bac.bac_x = X(:,1);
        bac.bac_y = X(:,2);
    end
    bac.bac_x = bac_x;
    bac.bac_y = bac_y;
    
end
end