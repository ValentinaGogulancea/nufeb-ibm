S2 = 2e-4; % oxygen concentration
ks1 = 3.94E-09; % ammonia
ks2 = 1.88E-06; % oxygen
S1 = linspace(1.5e-20,5e-7,100); % ammonia concentration
mu_max = 0.031; % h^-1

mu = mu_max * (S2/(ks2 + S2))* (S1./(ks1 + S1)) ;% * (S1/(ks1 + S1))

A = gradient(mu);

plot(S1,mu)

