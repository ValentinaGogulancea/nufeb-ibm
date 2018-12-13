function F = DiffReaction(AllConc,R)

[NRVliq,~,~] = react_term(AllConc,R);
[~,~,A,BB2,BB3] = DiffMatricesEuler(R, 1);
Sbc_Dir = 1000*kron(R.Sxy.Sbc_Dir(1:numStVLiq2), ones(R.Sxy.nT,1)); 
F = A*AllConc + BB2 * Sbc_Dir + BB3 *NRVliq;

end
