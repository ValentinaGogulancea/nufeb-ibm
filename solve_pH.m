function [spcM, Sh] = solve_pH(Sh_ini, StV, Keq, chrM)
w = 1;
% the pH solver

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
