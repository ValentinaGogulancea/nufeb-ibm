bac_n = 100; % number of cells
k = 0.5; % distance allowed between them (1 = basically linked)
bac_rmax = 5e-7; % radius
maxxSys = 2e-4; % length of the domain
n_max_row = floor(maxxSys/k/bac_rmax/2); % how many fit in a row
n_need_column = ceil (bac_n/n_max_row); % how many columns will one need

% bac_x = zeros(n_max_row * n_need_column,1); 
bac_y = zeros(n_max_row * n_need_column,1);

aux = 0: 2*k*bac_rmax: maxxSys;
% bac_x (1:n_max_row) = aux(1:end-1)' + bac_rmax; % first row
b = aux(1:end-1)' + bac_rmax;

bac_x = repmat(b,n_need_column,1); % all other rows will have the same x-coordinates

bac_y(1:n_max_row) = 2*k*bac_rmax*ones(n_max_row,1); % first row

for i = 1:n_need_column-1
bac_y(i* n_max_row + 1: (i+1)* n_max_row) = (2*i+1)*bac_rmax*ones(n_max_row,1)*k*2;
end
bac_x = bac_x(1:bac_n);
bac_y = bac_y(1:bac_n);


X = [bac_x bac_y]; Y = X;

bac_r = bac_rmax*ones(bac_n,1);
bac_m = bac_rmax*ones(bac_n,1);
s_dist = 2*bac_r;
overlap = 0.1*bac_rmax;
tic
X = TestOverlap(bac_x,bac_y,bac_r,bac_m,overlap,maxxSys,bac_rmax);
toc
bac_x1 = X(:,1);
bac_y1 = X(:,2);

figure
 for k =1:length(bac_x)    
        
            rectangle('Curvature',[1 1],'Position',[bac_x1(k)-bac_rmax bac_y1(k)-bac_rmax 2*bac_rmax 2*bac_rmax],'LineWidth',1);
 end

    axis equal;

   
tic
[bac_x, bac_y] = bac_shovingloops(bac_x, bac_y, bac_r, bac_m, 1, s_dist, overlap, bac_rmax, maxxSys);
toc
figure
 for k =1:length(bac_x)    
        
            rectangle('Curvature',[1 1],'Position',[bac_x(k)-bac_rmax bac_y(k)-bac_rmax 2*bac_rmax 2*bac_rmax],'LineWidth',1);
 end

    axis equal; 

