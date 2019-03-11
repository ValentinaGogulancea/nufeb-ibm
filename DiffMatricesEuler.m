function [R, Fb_layer, A, BB2, BB3] = DiffMatricesEuler(R, Fdiv)
% the matrices are updated after each diffusion event
persistent Lxy Bc

Fb_layer = 0;

if Fdiv == 1
    [R.Sxy, R.kTr, Fb_layer, Bc, Lxy] = b_layer(Bc, Lxy, R.Sxy, R.kTr, R.St, R.bac.bac_y); % commented bellow
    aux_x= R.Sxy.x_pn-R.Sxy.dx/2; aux_x1 = R.Sxy.x_pn + R.Sxy.dx/2; % coordinates of the beginning and the end of a grid cell
    aux_y= R.Sxy.y_pn-R.Sxy.dy/2; aux_y1 = R.Sxy.y_pn + R.Sxy.dy/2;
    c = zeros(R.bac.bac_n,R.Sxy.nT);
    for i = 1: R.Sxy.nT
        % this is to check if bacteria i is found in grid cell j
        c1 = (R.bac.bac_x >= aux_x(i)).*(R.bac.bac_x <= aux_x1(i)); % x-coordinate of bacteria i is between the edges of the cell
        c2 = (R.bac.bac_y >= aux_y(i)).*(R.bac.bac_y <= aux_y1(i)); % y-coordinate
        c(:,i) = c1.*c2;
    end
    R.Sxy.c = c; %this tells us if the centre of bacteria x is in grid defined by i-i+1 and j-j+1
end

St = R.St; 
Sxy = R.Sxy;

if Fdiv == 1
    % returns the rates of growth for bacteria, reaction rates for solubles
    R = massBal( Sxy.StVLiq2, Sxy.StVGas, R.bac.bac_m, R); % commented in the corresponding function
% updated reaction rates    

end


ind = [Sxy.nT*((1:St.numStVLiq2)'-1)+1,Sxy.nT*(1:St.numStVLiq2)'];
% tells one from where to where are the concentrations of species 1,2, etc
% in each grid cell of the biofilm + b- layer domain

aux = Sxy.nT; % number of grid cells
aux2 = aux*St.numStVLiq2; % number of concentration values 
a = R.kTr.a; % diffusion coeff/2 
% full of zeros
A = sparse(aux2, aux2); BB2 = sparse(aux2, aux2);  BB3 = sparse(aux2, aux2);
[I,J] = ind2sub([aux, aux],1:aux*aux);
for k = 1:St.numStVLiq2
    aa = 2 * a(ind(k,1):ind(k,2),ind(k,1):ind(k,2));
    A1 = aa.*Lxy; % the term at t+1 2 * 
    B2 = aa.*Bc; % the coefficients for boundary conditions of Dirichlet type 2* aa.*
    B3 = sparse(1:aux,1:aux,ones(Sxy.nT,1),aux,aux); % identity sparse matrix
    A = A + sparse((k-1)*aux+I,(k-1)*aux+J,reshape(A1, [],1),aux2,aux2);
    BB2 = BB2 + sparse((k-1)*aux+I,(k-1)*aux+J,reshape(B2, [],1),aux2,aux2);
    BB3 = BB3 + sparse((k-1)*aux+I,(k-1)*aux+J,reshape(B3, [],1),aux2,aux2);
end
end

function [Sxy, kTr, Fb_layer, Bc, Lxy] = b_layer(Bc, Lxy, Sxy, kTr, St, bac_y)
% establishes where the boundary layer + biofilm are
y_biof = [min(bac_y), max(bac_y)]; % corodinates of the biofilm minimum and maximum height

dx = Sxy.dx; dy = Sxy.dy; % discretization size

maxy =  y_biof(2) + dy + Sxy.T_blayer; % maximum height of the biofilm + boundary layer + one discretization cell higher

pos_maxy = find((Sxy.ySys >= maxy),1); % finds the position of the y coordinates corresponding to this point in the discretization

if isempty(pos_maxy)
    pos_maxy = length(Sxy.ySys); % if the maximum is not reached - position of the height of the computational domain
end
x = Sxy.xSys;  % number of points on Ox                 
y = Sxy.ySys(1:pos_maxy);% number of points on Oy    

nx = length(x)-1+2; % number of coordinates on Ox

ny = length(y)-1+2; % number of coordinates on Oy - limited to biofilm + b-layer

pos_xySys_old = Sxy.pos_xySys; % keep the position of the old biofilm + b-layer

ax = ones(Sxy.nxSys-2,1); % auxiliary variables used for the computation of Sxy

ay = zeros(Sxy.nySys-2,1);

ay(1:ny-2) = 1;

% Sxy.pos_xySys - length = number of grid cells in the discretization
% vector of 0 & 1's that shows which grid cells are in the bulk (0) and which are in the biofilm + b-layer (1)
Sxy.pos_xySys = kron(ax, ay); 
%
Fb_layer = sum(ne(Sxy.pos_xySys, pos_xySys_old)) > 0; % checks if the new b-layer is different from the old one

if Fb_layer == 1 % if it is 
    
    nT = (nx-2)*(ny-2); % number of grid cells
    
    xp = x + dx/2; % half intervals
    
    Sxy.x_pn = kron(xp(1:end-1)', ones(ny-2,1)); % x-coordinates of the center of a grid cell
    
    yp = y + dy/2;
    
    Sxy.y_pn = kron(ones(nx-2,1),yp(1:end-1)'); % y-coordinates of the center of a grid cell
    
    Sxy.nx = nx; 
    Sxy.ny = ny; 
    Sxy.nT = nT; 
    Sxy.x = x; 
    Sxy.y = y;
    
    % multiplies & extends the boundary conditions for all the biofilm + b-layer domain
    S = kron([Sxy.Sbc_Dir;St.StVLiq(St.numStVLiq2:St.numStVLiq)], ones(Sxy.nTSys,1));
    
    if ne(pos_xySys_old, 555)% if it's not the first time it is run - use for all the initial concentrations
        S(kron(ones(St.numStVLiq,1),pos_xySys_old)==1) = Sxy.StVLiq; % in the places where it was before - use the existing computed concentrations 
    end
    
    Sxy.StVLiq = S(kron(ones(St.numStVLiq,1),Sxy.pos_xySys)==1); % values for the soluble species concentrations in biofilm + b-layer
    Sxy.StVGas = Sxy.StVLiq(St.numStVLiq2*Sxy.nT+1:end );  % values for the gas species concentrations in biofilm + b-layer
    Sxy.StVLiq2 = Sxy.StVLiq(1:St.numStVLiq2*Sxy.nT); % values for the liquid species concentrations in biofilm + b-layer
    
    S = Sxy.pHnT;
    
    if ne(pos_xySys_old, 555)
        S(pos_xySys_old==1) = Sxy.pH;
    end
    Sxy.pH = S(Sxy.pos_xySys==1); % use the values of pH only for biofilm + b-layer grid cells

    kTr.Diffwater = kron(kTr.Diffn, ones(Sxy.nT,1));% diffusion coefficients in biofilm+BL

    S = kron(kTr.Diffn, ones(Sxy.nTSys,1)); % extend these values to the whole biofilm + b-layer
    
    if ne(pos_xySys_old, 555)
        S(kron(ones(St.numStVLiq2,1),pos_xySys_old)==1) = kTr.Diff;
    end
    kTr.Diff = S(kron(ones(St.numStVLiq2,1),Sxy.pos_xySys)==1);
    
    % this is the part that assembles the laplacian matrices for
    % computations - the left hand side of the diffusion reaction equation

    Ix = speye(nx-2);
    Iy = speye(ny-2);
    
    %Laplacian Matrix
    
    Ex1 = sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
    Ax = Ex1+Ex1'-2*Ix;
    Ey1 = sparse(2:ny-2,1:ny-3,1,ny-2,ny-2);
    Ay = Ey1+Ey1'-2*Iy;
    %Dirichlet boundary conditions
    %%Wall
    Ay(1,1)=-1;
    %%Cyclic boundary conditions in the sides
    Ax(1,nx-2)=1;  Ax(nx-2,1)=1;
    Lxy = kron( Ix, Ay/dy^2) + kron( Ax/dx^2, Iy);% they go at Ay and Ay respectively

    %Boundary conditions Matrix
    bcy = zeros(ny-2, ny-2);
    bcx = zeros((nx-2),(nx-2));
    %%Dirichlet boundary condition
    bcy(ny-2,ny-2) = 1; 
    % the other laplacian that must be added for boundary conditions
    % updating in the diff reaction solving algorithm

    Bc = kron(Ix, bcy/dy^2 ) + kron( bcx/dx^2 , Iy); % tells where the conditions are Dirichlet and where they are not /dy^2 /dx^2
end
end