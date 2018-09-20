clc
clear all
close all
% the first column is the values on
% 3D representation - X = fixed
% for all the files
folder='E:\NUFEB\fph-dbulk-4mg';
filetype='*.h5';  % or xlsx
f = fullfile(folder,filetype);
d = dir(f);

for k = 1:numel(d)
    filename = fullfile(folder,d(k).name);
    % fid = H5F.open(filename);
    % h5disp(filename); % do it once to see how the hdf5 file was saved
    % the values for the free energy of anabolism/catabolism - as value corresponding to coordinate point
    aux_time = num2str(filename);
    aux_time = aux_time(1:end-3);
    aux_time = aux_time(29:end);
    time1(k) = str2num(aux_time)/3600;
end
% because Matlab reads the files alphabetically - some rearrangements are
% required
[time,Index] = sort(time1);

% initialization
aux_Ana_AOB{length(time)} = [];
aux_Ana_NOB{length(time)} = [];
%
aux_Cat_AOB{length(time)} = [];
aux_Cat_NOB{length(time)} = [];
%
% concentration fields
aux_NH3_concfield{length(time)} = [];
aux_NO3_concfield{length(time)} = [];
aux_NO2_concfield{length(time)} = [];
aux_O2_concfield{length(time)} = [];
aux_CO2_concfield{length(time)} = [];
% pH field
aux_Hplus{length(time)} = [];
% and now for the particle plotting

% If we want to track one special particle
aux_ID_tracking{length(time)} = [];

% the type of particle
aux_Type_part{length(time)} = [];
% the coordinates of particles
aux_x_coord{length(time)} = [];
aux_y_coord{length(time)} = [];
aux_z_coord{length(time)} = [];
aux_radius{length(time)} = [];
%
for k = 1:numel(d)
    filename = fullfile(folder,d(k).name);
    aux_Ana_AOB{k} = h5read(filename,'/anabolism/aob');
    aux_Ana_NOB{k} = h5read(filename,'/anabolism/nob');
    %
    aux_Cat_AOB{k} = h5read(filename,'/catabolism/aob');
    aux_Cat_NOB{k} = h5read(filename,'/catabolism/nob');
    %
    % concentration fields
    aux_NH3_concfield{k} = h5read(filename,'/concentration/nh3');
    aux_NO3_concfield{k} = h5read(filename,'/concentration/no3');
    aux_NO2_concfield{k} = h5read(filename,'/concentration/no2');
    aux_O2_concfield{k} = h5read(filename,'/concentration/o2');
    aux_CO2_concfield{k} = h5read(filename,'/concentration/co2');
    % pH field
    aux_Hplus{k} = h5read(filename,'//hydronium');
    % and now for the particle plotting
    
    % If we want to track one special particle
    aux_ID_tracking{k} = h5read(filename,'//id');
    
    % the type of particle
    aux_Type_part{k} = h5read(filename,'//type');
    % the coordinates of particles
    aux_x_coord{k} = h5read(filename,'//x');
    aux_y_coord{k} = h5read(filename,'//y');
    aux_z_coord{k} = h5read(filename,'//z');
    aux_radius{k} = h5read(filename,'//radius');
end
Ana_AOB{length(time)} = [];
Ana_NOB{length(time)} = [];
%
Cat_AOB{length(time)} = [];
Cat_NOB{length(time)} = [];
%
% concentration fields
NH3_concfield{length(time)} = [];
NO3_concfield{length(time)} = [];
NO2_concfield{length(time)} = [];
O2_concfield{length(time)} = [];
CO2_concfield{length(time)} = [];
% pH field
Hplus{length(time)} = [];
% If we want to track one special particle
ID_tracking{length(time)} = [];
% the type of particle
Type_part{length(time)} = [];

% the type of particle
% the coordinates of particles
x_coord{length(Index)} = [];
y_coord{length(Index)} = [] ;
z_coord{length(Index)} = [];
radius{length(Index)} = [];
for k = 1:length(Index)
    % SORTED CONCENTRATIONS
    O2_concfield{k} =  aux_O2_concfield{Index(k)};
    CO2_concfield{k} = aux_CO2_concfield{Index(k)};
    NO2_concfield{k} = aux_NO2_concfield{Index(k)};
    NO3_concfield{k} = aux_NO3_concfield{Index(k)};
    NH3_concfield{k} = aux_NH3_concfield{Index(k)};
    %SORTED ANABOLISM AND CATABOLISM FREE ENERGIES
    Ana_AOB{k} = aux_Ana_AOB{Index(k)};
    Ana_NOB{k} = aux_Ana_NOB{Index(k)};
    %
    Cat_AOB{k} = aux_Ana_AOB{Index(k)};
    Cat_NOB{k} = aux_Ana_NOB{Index(k)};
    %SORTED pH
    Hplus{k} = aux_Hplus{Index(k)};
    %SORTED BACTERIA COORDINATES & ALL
    Type_part{k} = aux_Type_part{Index(k)};
    ID_tracking{k} = aux_ID_tracking{k};
    x_coord{k} = aux_x_coord{Index(k)};
    y_coord{k} = aux_y_coord{Index(k)} ;
    z_coord{k} = aux_z_coord{Index(k)};
    radius{k} = aux_radius{Index(k)};
end
%%%
%%%
lmax_x = 120; lmax_y = 40; lmax_z = 300;
lmin_x = 0; lmin_y = 0; lmin_z = 0;
x_conc = linspace(lmin_x,lmax_x,40);
y_conc = linspace(lmin_y,lmax_y,20);
z_conc = linspace(lmin_z,lmax_z,100);
%%%%%%%% The actual plotting ------------------------------------------
%
% % now for making an average on the biofilm
% % see the biofilm max height
b_layer = 40;
% % The Oxygen Profile vs time
% average_O2 = zeros(1,numel(d));
% index_height = zeros(numel(d),1);
% for k = 1: numel(d)
%     [~,index_height(k)] = min(abs(max(z_coord{k})+ b_layer - z_conc));
%     
%     n_voxels = length(x_conc)*length(y_conc)*index_height(k);
%     average_auxO2 = zeros(index_height(k),1);
%     for i = 1: index_height(k)
%         average_auxO2(i) =  sum(sum(O2_concfield{k}(:,:,i)));
%     end
%     average_O2(k) = sum(average_auxO2)/n_voxels;
% end
% 
% figure
% plot(time, average_O2,'LineWidth',2)
% xlabel('time,h')
% ylabel('O_2 concentration')
% % The NH3 profile over time
% average_NH3 = zeros(1,numel(d));
% for k = 1: numel(d)
%     [~,index_height(k)] = min(abs(max(z_coord{k})+ b_layer - z_conc));
%     
%     n_voxels = length(x_conc)*length(y_conc)*index_height(k);
%     average_auxNH3 = zeros(index_height(k),1);
%     for i = 1: index_height(k)
%         average_auxNH3(i) =  sum(sum(NH3_concfield{k}(:,:,i)));
%     end
%     average_NH3(k) = sum(average_auxNH3)/n_voxels;
% end
% %
% figure
% plot(time, average_NH3,'LineWidth',2)
% xlabel('time,h')
% ylabel('NH_3 concentration')
% % NO2
% average_NO2 = zeros(1,numel(d));
% for k = 1: numel(d)
%     [~,index_height(k)] = min(abs(max(z_coord{k})+ b_layer - z_conc));
%     
%     n_voxels = length(x_conc)*length(y_conc)*index_height(k);
%     average_auxNO2 = zeros(index_height(k),1);
%     for i = 1: index_height(k)
%         average_auxNO2(i) =  sum(sum(NO2_concfield{k}(:,:,i)));
%     end
%     average_NO2(k) = sum(average_auxNO2)/n_voxels;
% end
% figure
% plot(time, average_NO2,'LineWidth',2)
% xlabel('time,h')
% ylabel('NO_2 concentration')
% %NO3
% average_NO3 = zeros(1,numel(d));
% for k = 1: numel(d)
%     [~,index_height(k)] = min(abs(max(z_coord{k})+ b_layer - z_conc));
%     
%     n_voxels = length(x_conc)*length(y_conc)*index_height(k);
%     average_auxNO3 = zeros(index_height(k),1);
%     for i = 1: index_height(k)
%         average_auxNO3(i) =  sum(sum(NO3_concfield{k}(:,:,i)));
%     end
%     average_NO3(k) = sum(average_auxNO3)/n_voxels;
% end
% figure
% plot(time, average_NO3,'Linewidth',2)
% xlabel('time,h')
% ylabel('NO_3 concentration')
% %CO2
% average_CO2 = zeros(1,numel(d));
% for k = 1: numel(d)
%     [~,index_height(k)] = min(abs(max(z_coord{k})+ b_layer - z_conc));
%     
%     n_voxels = length(x_conc)*length(y_conc)*index_height(k);
%     average_auxCO2 = zeros(index_height(k),1);
%     for i = 1: index_height(k)
%         average_auxCO2(i) =  sum(sum(CO2_concfield{k}(:,:,i)));
%     end
%     average_CO2(k) = sum(average_auxCO2)/n_voxels;
% end
% figure
% plot(time, average_CO2,'LineWidth',2)
% xlabel('time,h')
% ylabel('CO_2 concentration')
%%%% the average number of bacteria as well

n_AOB = zeros(1,numel(d)); n_NOB = zeros(1,numel(d)); n_dead = zeros(1,numel(d));
for k = 1: numel(d)
    for i = 1:length(Type_part{k})
        if Type_part{k}(i) == 1
            n_AOB(k) = n_AOB(k) + 1;
        elseif Type_part{k}(i) == 2
            n_NOB(k) = n_NOB(k) + 1;
        else
            n_dead(k) = n_dead(k) + 1;
        end
    end
end

figure
plot(time, n_AOB)
hold on
plot(time, n_NOB)
xlabel('time,h')
ylabel('number of cells')
legend('AOB','NOB')
% average values
% now lets see how the 2D plotting
% O2 = O2_concfield{20};
%
% O2_ox = reshape(O2(:,1,1),[40,1]);
%
% O2_oz = reshape (O2(1,1,:),[100,1]);
%
% %
% O2_oy = reshape (O2(1,:,1),[20,1]);
% figure
% plot(x_conc,O2_ox)
% figure
% plot(y_conc,O2_oy)
% figure
% plot(z_conc,O2_oz)
%%%%%%Plotting and saving the videos
NH3fig = figure('Name','NH3','Color','w');
[X,Z] = ndgrid(x_conc,z_conc);
%at the 5 um distance
v1 = VideoWriter('NH3.avi');
v1.FrameRate = 15;
v1.Quality  = 75;
open(v1);

aux = [];
for i = 1:numel(d)
    aux = [aux NH3_concfield{i}];
end
NH3max = max(aux(:));
NH3min = min(aux(:));
for k = 1:numel(d)
    NH3 = reshape(NH3_concfield{k}(:,1,:),[40,100]);
    figure(NH3fig)
    surf(X,Z,NH3,'EdgeColor','none');%C = 
%     clabel(C,'LabelSpacing',100,'FontSize',8,'LineStyle',':')
    az = 0;
    el = 90;
    view(az, el);
    shading interp
    colormap(parula(256));
    caxis([0.5*NH3min NH3max])
    c = colorbar;
    c.Label.String = 'Concentration (mol/L)';
   
    title({'NH_3 concentration profile';['time (\itt) = ',num2str(time(k))]})
    xlabel('Spatial co-ordinate (x), \mum')
    ylabel(' Spatial co-ordinate (z), \mum')
    axis ([0 lmax_x 0 lmax_z])%NH3min NH3max])
    frame = getframe(NH3fig);
    writeVideo(v1,frame);
    
end
close(v1);

v2 = VideoWriter('NO2.avi');
open(v2);
aux = [];
for i = 1:numel(d)
    aux = [aux NO2_concfield{i}];
end
NO2max = max(aux(:));
NO2min = min(aux(:));
NO2fig = figure('Name','NO2','Color','w');
for k = 1:numel(d)
    NO2 = reshape(NO2_concfield{k}(:,1,:),[40,100]);
    figure(NO2fig)
    surf(X,Z,NO2,'EdgeColor','none')
    shading interp
    az = 0;
    el = 90;
    view(az, el);
    axis ([0 lmax_x 0 lmax_z])
    caxis([NO2min NO2max])
     colormap(parula(512));
    c = colorbar;
    c.Label.String = 'Concentration (mol/L)';
    title({'NO_2 concentration profile';['time (\itt) = ',num2str(time(k))]})
     xlabel('Spatial co-ordinate (x), \mum')
    ylabel(' Spatial co-ordinate (z), \mum')
    frame = getframe(NO2fig);
    writeVideo(v2,frame);
end
colorbar('Limits',[NO2min NO2max]);
close(v2);

%
v3 = VideoWriter('NO3.avi');
open(v3);
aux = [];
for i = 1:numel(d)
    aux = [aux NO3_concfield{i}];
end
NO3max = max(aux(:));
NO3min = min(aux(:));
NO3fig = figure('Name','NO3','Color','w');
for k = 1:numel(d)
    NO3 = reshape(NO3_concfield{k}(:,1,:),[40,100]);
    figure(NO3fig)
    surf(X,Z,NO3,'EdgeColor','none')
    shading interp
    az = 0;
    el = 90;
    view(az, el);
    caxis([NO3min NO3max])
     colormap(parula(512));
    c = colorbar;
    c.Label.String = 'Concentration (mol/L)';
    title({'NO_3 concentration profile';['time (\itt) = ',num2str(time(k))]})
     xlabel('Spatial co-ordinate (x), \mum')
    ylabel(' Spatial co-ordinate (z), \mum')
    axis ([0 lmax_x 0 lmax_z])
    frame = getframe(NO3fig);
    writeVideo(v3,frame);
end
close(v3);

v4 = VideoWriter('O2.avi');
open(v4);
aux = [];
for i = 1:numel(d)
    aux = [aux O2_concfield{i}];
end
O2max = max(aux(:));
O2min = min(aux(:));
O2fig = figure('Name','O2','Color','w');
for k = 1:numel(d)
    O2 = reshape(O2_concfield{k}(:,1,:),[40,100]);
    figure(O2fig)
    surf(X,Z,O2,'EdgeColor','none')
    shading interp
    az = 0;
    el = 90;
    view(az, el);
     colormap(parula(512));
    caxis([O2min O2max])
    c = colorbar;
    c.Label.String = 'Concentration (mol/L)';
    title({'O_2 concentration profile';['time (\itt) = ',num2str(time(k))]})
     xlabel('Spatial co-ordinate (x), \mum')
    ylabel(' Spatial co-ordinate (z), \mum')
    axis ([0 lmax_x 0 lmax_z])
    frame = getframe(O2fig);
    writeVideo(v4,frame);
end
close(v4);

v5 = VideoWriter('CO2.avi');
open(v5);
aux = [];
for i = 1:numel(d)
    aux = [aux CO2_concfield{i}];
end
CO2max = max(aux(:));
CO2min = min(aux(:));
CO2fig = figure('Name','CO2','Color','w');
for k = 1:numel(d)
    CO2 = reshape(CO2_concfield{k}(:,1,:),[40,100]);
    figure(CO2fig)
    surf(X,Z,CO2,'EdgeColor','none')
    shading interp
    az = 0;
    el = 90;
    view(az, el);
    colormap(parula(512));
    caxis([CO2min CO2max])
    c = colorbar;
    c.Label.String = 'Concentration (mol/L)';
    title({'CO_2 concentration profile';['time (\itt) = ',num2str(time(k))]})
     xlabel('Spatial co-ordinate (x), \mum')
    ylabel(' Spatial co-ordinate (z), \mum')
    axis ([0 lmax_x 0 lmax_z])
    frame = getframe(CO2fig);
    writeVideo(v5,frame);
end
close(v5);


