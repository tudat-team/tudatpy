%
% This script processes the results of the CEC2013 test problems using
% the Pagmo algorithms, run by the cec2013OptimizerComparison.cpp 
% Tudat/Pagmo2 example.
%

set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

clc
close all
clear all

saveFolder = '../SimulationOutput/';

% Define colors.
optima_2_dim = cell(3,1);
plotOptions = cell(3,1);
plotOptions{1} = 'b';
plotOptions{2} = 'r';
plotOptions{3} = 'k';

% Load data for 3 different population sizes/number of generations; 2 dimensional problem optimization
for k=1:3
    optima_2_dim{k} = load(strcat(saveFolder,'cec2013Optima_0_',num2str(k-1),'.dat'));
end

%  Plot data for 2-dimensional optimization
figure(1)
for k=1:3
    
    % Determine minimum
    [Y,I]=min(optima_2_dim{k}');
    
    % Plot results for each test problem
    for j=1:28
        if( j == 1 )
            subplotIndex = j;
        else
            subplotIndex = j + floor((j-1)/7);
        end
        
        subplot(4,8,subplotIndex)
        if( j == 1)
            plotHandle{k}=plot(optima_2_dim{k}(j,:),strcat(plotOptions{k},'-*'));
        else
            plot(optima_2_dim{k}(j,:),strcat(plotOptions{k},'-*'))
        end
        hold on
        grid on
        
        % Add circle at optimum
        scatter(I(j),Y(j),50*k,strcat(plotOptions{k},'o'))
        xlim([0 12])
        title(strcat('Problem ',{' '},num2str(j)))
        if(rem(j,7)==1)
            ylabel('Minimum [-]')
        end
        if(j>21)
            xlabel('Optimizer [-]')
        end
        
    end
end

% Add legend
subplot(4,8,1)
legend([plotHandle{1} plotHandle{2} plotHandle{3}],'Pop: 16, Gen: 64', 'Pop: 32, Gen: 32', 'Pop: 64, Gen: 16'  )

% Resize and add super-title
set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 60 40]);
set(gcf,'PaperPositionMode','auto');
suptitle('Optima, CEC2013 2-D test problems')

subplot(1,8,8)
axis off
text(0.5,0.5,{'1: DE','2: SADE','3: DE1220','4: PSO','5: SEA','6: SGA','7: SA','8: BC','9: CMAES','10: IHS','11: XNES'},'FontSize',18);

pause(0.1)
figure(gcf)
saveas(gcf,strcat('cec2013_PagmoOptimization_2D'),'png');

%%
for k=1:3
    optima_5_dim{k} = load(strcat(saveFolder,'cec2013Optima_1_',num2str(k-1),'.dat'));
end

figure(2)
for k=1:3
    [Y,I]=min(optima_5_dim{k}');
    for j=1:28

        if( j == 1 )
            subplotIndex = j;
        else
            subplotIndex = j + floor((j-1)/7);
        end
        
        subplot(4,8,subplotIndex)
        
        if( j == 1)
            plotHandle{k}=plot(optima_5_dim{k}(j,:),strcat(plotOptions{k},'-*'));
        else
            plot(optima_5_dim{k}(j,:),strcat(plotOptions{k},'-*'))
        end
        
        hold on
        grid on
        scatter(I(j),Y(j),50*k,strcat(plotOptions{k},'o'))
        xlim([0 12])
        title(strcat('Problem ',{' '},num2str(j)))
        if(rem(j,7)==1)
            ylabel('Minimum [-]')
        end
        if(j>21)
            xlabel('Optimizer [-]')
        end
        
    end
end


subplot(4,8,1)
legend([plotHandle{1} plotHandle{2} plotHandle{3}],'Pop: 16, Gen: 64', 'Pop: 32, Gen: 32', 'Pop: 64, Gen: 16' )

suptitle('Optima, CEC2013 5-D test problems')

set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 60 40]);
set(gcf,'PaperPositionMode','auto');


subplot(1,8,8)
axis off
text(0.5,0.5,{'1: DE','2: SADE','3: DE1220','4: PSO','5: SEA','6: SGA','7: SA','8: BC','9: CMAES','10: IHS','11: XNES'},'FontSize',18);

pause(0.1)
figure(gcf)
saveas(gcf,strcat('cec2013_PagmoOptimization_5D'),'png');


%%
for k=1:3
    optima_10_dim{k} = load(strcat(saveFolder,'cec2013Optima_2_',num2str(k-1),'.dat'));
end

figure(3)
for k=1:3
    [Y,I]=min(optima_10_dim{k}');
    for j=1:28
        
        if( j == 1 )
            subplotIndex = j;
        else
            subplotIndex = j + floor((j-1)/7);
        end
        
        subplot(4,8,subplotIndex)
        
        if( j == 3)
            plotHandle{k}=plot(optima_10_dim{k}(j,:),strcat(plotOptions{k},'-*'));
        else
            plot(optima_10_dim{k}(j,:),strcat(plotOptions{k},'-*'))
        end
        hold on
        grid on
        scatter(I(j),Y(j),50*k,strcat(plotOptions{k},'o'))
        xlim([0 12])
        title(strcat('Problem ',{' '},num2str(j)))
        if(rem(j,7)==1)
            ylabel('Minimum [-]')
        end
        if(j>21)
            xlabel('Optimizer [-]')
        end
        
    end
end

subplot(4,8,1)
legend([plotHandle{1} plotHandle{2} plotHandle{3}],'Pop: 16, Gen: 64', 'Pop: 32, Gen: 32', 'Pop: 64, Gen: 16'  )

suptitle('Optima, CEC2013 10-D test problems')

set(gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 60 40]);
set(gcf,'PaperPositionMode','auto');

subplot(1,8,8)
axis off
text(0.5,0.5,{'1: DE','2: SADE','3: DE1220','4: PSO','5: SEA','6: SGA','7: SA','8: BC','9: CMAES','10: IHS','11: XNES'},'FontSize',18);

pause(0.1)
figure(gcf)
saveas(gcf,strcat('cec2013_PagmoOptimization_10D'),'png');
