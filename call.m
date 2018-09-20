% %%%% Script call.m: Initializates the Simulation and calls the numerical
% solver

tic
clc
close all
clear all %#ok<CLSCR>
% prompt = '>> Do you want to create a new R? (default [Y]) [Y = 1/N = 2]\n>> Answer: ';
% x = input(prompt);
% if isempty(x)
x = 1;
% end
% if x == 1
%     fprintf('--> YES\n')
% else
%     fprintf('--> NO\n')
% end

if x == 1
% %%%% Reading the model parameters from Excel file: calling the subscript
% %%%% loadModelXlsx.m that generates the initial structures
    R = loadModelXlsx;
    % save R in a folder with its name on it
else
    load('R.mat')
    R.Sxy.pos_xySys = 555;
end
clear x prompt
fprintf('\n> MODEL RUNNING >>>>>\n')

%      R = integTime1(R);
     R = integTimeEuler(R);
fprintf('\n \n SIMULATION FINISHED >>>> \n')

for i = 1:(R.St.numX)
fprintf('\nNumber of cells of type %d Name %s : %d\n', i,char(R.rm.rNamesX(i)),R.bac.bac_ns(i))
end

% fprintf('>>> Push any key to continue with the visualization of results...\n\n')
% pause()


fprintf('\n \n SIMULATION FINISHED >>>> \n  \n ')

for i = 1:(R.St.numX)
fprintf('\nNumber of cells of type %d Name %s : %d\n', i,char(R.rm.rNamesX(i)),R.bac.bac_ns(i))
end
fprintf('\nTOTAL Number of cells: %d\n\n\n', R.bac.bac_n)
clear i

%%%%%%% VISUALIZATION

% draw
