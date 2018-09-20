% check the amount of Na to add

R = loadModelXlsx;
Sh_ini = 10.^(-R.pOp.pH);
Keq = R.pTh.Keq;
chrM = R.pTh.chrM; 

   
StV = [R.St.StVLiq;1;0];  % initial concentrations
StV(6)  = 1e-4;
StV(7)  = 1e-5;
[~, Sh] = solve_pH(Sh_ini,StV, Keq, chrM);

ph = -log10(Sh)