%
% This script processes the results of the targeting optimization problem
% run by the propagationTargetingExample.cpp Tudat/Pagmo2 example.
%

set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

clc
clear all
close all

% Load grid search data for unperturbed case
dataFolder = '../SimulationOutput/';
load(strcat(dataFolder,'propagationTargetingGridSearch.dat'))
load(strcat(dataFolder,'propagationTargetingGridSearch_x_data.dat'))
load(strcat(dataFolder,'propagationTargetingGridSearch_y_data.dat'))

% Plot objective function as contour plot for unperturbed case
figure(1)
contour(propagationTargetingGridSearch_x_data,propagationTargetingGridSearch_y_data,propagationTargetingGridSearch')
colorbar
xlabel('Argument of periapsis [deg]')
ylabel('Longitude of asc. node [deg]')
title('Minimum targeting error [m]')

set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
set(gcf,'PaperPositionMode','auto');

pause(0.1)
saveas(gcf,strcat('unperturbedTargetingGridSearch'),'png');
%%

figure(2)
population = cell(8,1);
fitness = cell(8,1);

% Plot population cloud over contour plot for 8 generations (unperturbed case)
for k=1:8
    
    % Specify current generation
    if( k == 1 )
        indexToUse = 1;
    elseif( k == 2 )
        indexToUse = 2;
    elseif( k == 3 )
        indexToUse = 3;
    elseif( k == 4 )
        indexToUse = 4;
    elseif( k == 5 )
        indexToUse = 5;
    elseif( k == 6 )
        indexToUse = 8;
    elseif( k == 7 )
        indexToUse = 12;
    elseif( k == 8 )
        indexToUse = 25;
    end
    
    % Retrieve population/fitness for requested genertion
    j=indexToUse;
    population{j} = load(strcat(dataFolder,'population_targetingPropagation_',num2str(j-1),'_',num2str(j-1),'.dat'));
    fitness{j} = load(strcat(dataFolder,'fitness_targetingPropagation_',num2str(j-1),'_',num2str(j-1),'.dat'));
    
    % Plot current generation on contour plot
    subplot(2,4,k)
    contour(propagationTargetingGridSearch_x_data,propagationTargetingGridSearch_y_data,propagationTargetingGridSearch')
    hold on
    scatter(population{indexToUse}(:,1),population{indexToUse}(:,2),25,fitness{indexToUse}(:,1),'*')
    xlabel('Argument of periapsis [deg]')
    ylabel('Longitude of asc. node [deg]')
    title(strcat('Gen.=',{' '},num2str(indexToUse),{' '},', Min.=',num2str(min(fitness{indexToUse}(:,1))/1000,3),' km'));
end

set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
set(gcf,'PaperPositionMode','auto');

pause(0.1)
saveas(gcf,strcat('unperturbedTargetingOptimization'),'png');


% Load grid search data for perturbed case
load(strcat(dataFolder,'propagationTargetingGridSearch_pert.dat'))
load(strcat(dataFolder,'propagationTargetingGridSearch_pert_x_data.dat'))
load(strcat(dataFolder,'propagationTargetingGridSearch_pert_y_data.dat'))

%%

% Plot objective function as contour plot for perturbed case
figure(3)
contour(propagationTargetingGridSearch_pert_x_data,propagationTargetingGridSearch_pert_y_data,propagationTargetingGridSearch_pert')
colorbar
xlabel('Argument of periapsis [deg]')
ylabel('Longitude of asc. node [deg]')
title('Minimum targeting error [m]')

set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
set(gcf,'PaperPositionMode','auto');

pause(0.1)
saveas(gcf,strcat('perturbedTargetingGridSearch'),'png');

figure(4)

population_pert = cell(8,1);
fitness_pert = cell(8,1);

% Plot population cloud over contour plot for 4 generations (perturbed case)
for k=1:4
    
    if( k == 1 )
        indexToUse = 1;
    elseif( k == 2 )
        indexToUse = 2;
    elseif( k == 3 )
        indexToUse = 3;
    elseif( k == 4 )
        indexToUse = 4;
    end
    
    %Retrieve population/fitness for requested genertion
    j=indexToUse;
    population_pert{j} = load(strcat(dataFolder,'population_targetingPropagation_pert_',num2str(j-1),'_',num2str(j-1),'.dat'));
    fitness_pert{j} = load(strcat(dataFolder,'fitness_targetingPropagation_pert_',num2str(j-1),'_',num2str(j-1),'.dat'));
    
    % Plot current generation on contour plot
    subplot(1,4,k)
    contour(propagationTargetingGridSearch_pert_x_data,propagationTargetingGridSearch_pert_y_data,propagationTargetingGridSearch_pert')
    hold on
    scatter(population_pert{indexToUse}(:,1),population_pert{indexToUse}(:,2),25,fitness_pert{indexToUse}(:,1),'*')
    xlabel('Argument of periapsis [deg]')
    ylabel('Longitude of asc. node [deg]')
    min(fitness_pert{indexToUse}(:,1))
    title(strcat('Gen.=',{' '},num2str(indexToUse),{' '},', Min.=',num2str(min(fitness_pert{indexToUse}(:,1))/1000,3),' km'));
end

set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
set(gcf,'PaperPositionMode','auto');

pause(0.1)
saveas(gcf,strcat('perturbedTargetingOptimization'),'png');
