function R = boundary(R, dT, Flag)

Sxy = R.Sxy;
Dir_k = R.Inf.Dir_k;
numStVLiq = R.St.numStVLiq;
numStVLiq2 = R.St.numStVLiq2;
ind = [Sxy.nT*((1:numStVLiq2)'-1)+1,Sxy.nT*(1:numStVLiq2)'];    
invHRT = R.pOp.invHRT;

if ne(Flag,0)
G = zeros(numStVLiq2, 1);
for k = 1:numStVLiq2
    G(k) = (sum(R.rm.NRVliq(ind(k,1):ind(k,2)))/Sxy.nT); % compute the mean rates for all liquid species
end

options = odeset('RelTol',1e-8,'AbsTol',1e-10,'NonNegative',ones(numStVLiq2,1));% set options
[~,Y] = ode15s(@(t,y) massbal(t,y,G, R.Inf.St, R.pOp.Pollutant, R.flagN),[0 dT],R.Sxy.Sbc_Dir(1:numStVLiq2), options);% solve the mass balace
Reactor = Y(end, :)';

if ne(sum(Reactor < 0), 0)
    ap = find(Reactor < 0);
  for k = 1:length(ap)
      if ap(k) == 1
            Reactor(1) = 1e-20;
      elseif ap(k) == 2
          Reactor(2) = 1e-20;
      end
  end
      
end
end

aux = zeros(numStVLiq2, 1);
for k = 1:numStVLiq2
    if Dir_k(k) %Dirichlet  
        aux(k) = Sxy.Sbc_Dir(k);  % if the boundary is Dirichlet - use the constant concentration             
    else %Neumann
       
    if ne(Flag,0)
    aux(k) = Reactor(k); % this updates the concetration values using the mass balance
    else
        aux(k) = Sxy.Sbc_Dir(k);   % or not
    end
  
    end
end

Sbc_Dir_new = aux;

R.Sxy.Sbc_Dir(1:numStVLiq2) = Sbc_Dir_new;

% R.pOp.invHRT = invHRT;


function dy = massbal(~,y,G, Xi, Pollutant, flagN)
        % the mass balance over the reactor
dy = zeros(length(y),1);

if flagN == 2 && y(1) > Pollutant
    dy(1) = 0;
    invHRTaux = - G(1)/(Xi(1) - y(1));
    if ne(invHRTaux, 0)
        invHRT = invHRTaux;
    end
else
    dy(1) = invHRT*(Xi(1) - y(1)) + G(1);  
end
dy(2:end) = invHRT*(Xi(2:end) - y(2:end)) + G(2:end);  % 1/HRT * (cin-cbulk) + mean of consumption/formation
dy(4:end) = 0;
end

end