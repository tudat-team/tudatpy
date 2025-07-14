%
% This script processes the results of the Earth-Mars Lambert targeting
% multi-objective  optimization, run by the
% multiObjectiveEarthMarsTransferExample.cpp Tudat/Pagmo2 example.
%

set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

clc
clear all
close all

% Load grid search data
saveFolder = '../SimulationOutput/';
load(strcat(saveFolder,'himmelBlauGridSearch.dat'))
load(strcat(saveFolder,'himmelBlauGridSearch_x_data.dat'))
load(strcat(saveFolder,'himmelBlauGridSearch_y_data.dat'))

% Create grid search plot
contour(himmelBlauGridSearch_x_data,himmelBlauGridSearch_y_data,himmelBlauGridSearch',50)
xlabel('x [-]')
ylabel('y [-]')
title('Himmelblau function value [-]')
colorbar

set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
set(gcf,'PaperPositionMode','auto');

pause(0.1)
saveas(figure(1),'himmelblauGridSearch','png');
%%
% Plot data for 6 generations
figure(2)
population = cell(1,1);
fitness = cell(1,1);
for k=1:6
    
    % Specify current generation
    if( k == 1 )
        indexToUse = 1;
    elseif( k == 2 )
        indexToUse = 5;
    elseif( k == 3 )
        indexToUse = 10;
    elseif( k == 4 )
        indexToUse = 20;
    elseif( k == 5 )
        indexToUse = 50;
    elseif( k == 6 )
        indexToUse = 100;
    end
    
    % Load data for current generation, for each of the three optimizers
    for j=indexToUse
        population{j} = load(strcat(saveFolder,'population_himmelblau_',num2str(j),'.dat'));
        fitness{j} = load(strcat(saveFolder,'fitness_himmelblau_',num2str(j),'.dat'));
    end
    
    % Plot population for current generation/optimizer on top of porkchop
    subplot(2,3,k)
    contour(himmelBlauGridSearch_x_data,himmelBlauGridSearch_y_data,himmelBlauGridSearch',50)
    xlabel('x [-]')
    ylabel('y [-]')
    hold on
    scatter(population{indexToUse}(:,1),population{indexToUse}(:,2),25,fitness{indexToUse}(:,1),'*')
    title(strcat('Iteration ',{' '},num2str(indexToUse),{' '},', Minimum=',num2str(min(fitness{j}(:,1)),4)));
    
end
%%
set(gcf, 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
set(gcf,'PaperPositionMode','auto');

pause(0.1)
saveas(gcf,strcat('himmelblauOptimization'),'png');
