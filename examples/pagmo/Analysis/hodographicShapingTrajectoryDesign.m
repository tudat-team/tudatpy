%
% This script processes the results of the hodographic shaping trajectory design optimization
% run by the hodographicShapingTrajectoryExample.cpp Tudat/Pagmo2 example.
%

set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

clc
clear all
close all

% Manually define the colorbar to obtain plots consistent with those from
% Gondelach, 2015 (Hodographic-Shaping Method for Low-Thrust Interplanetary Trajectory Design)
manuallyDefinedColorbar = zeros(45,3); 
manuallyDefinedColorbar(1:7,1:3)= [2/256 29/256 242/256].*ones(7,3); 
manuallyDefinedColorbar(8,1:3)= [0/256 136/256 255/256].*ones(1,3); 
manuallyDefinedColorbar(9:10,1:3)= [1/256 171/256 255/256].*ones(2,3); 
manuallyDefinedColorbar(11:15,1:3)= [1/256 213/256 253/256].*ones(5,3);
manuallyDefinedColorbar(16:20,1:3)= [55/256 255/256 200/256].*ones(5,3);
manuallyDefinedColorbar(21:40,1:3)= [158/256 255/256 98/256].*ones(20,3);
manuallyDefinedColorbar(41:45,1:3)= [203/256 2/256 0/256].*ones(5,3);

% Load results
hodographicShapingLowOrder = importdata('../SimulationOutput/hodographicShapingLowOrder.dat');
hodographicShapingHighOrder = importdata('../SimulationOutput/hodographicShapingHigherOrder.dat');
hodographicShapingLowOrderOneRev = importdata('../SimulationOutput/hodographicShapingLowOrderOneRevolution.dat');

% Retrieve low-order solution departure dates
for i = 1:401 
    departureDateLowOrder(i,1) = hodographicShapingLowOrder(i,3);
end
% Retrieve low-order solution time-of-flights
for i = 1:301 
    TOFlowOrder(i,1) = hodographicShapingLowOrder( ( i - 1 ) * 401 + 1,2);
end
% Retrieve low-order solution deltaVs
for i = 1:301
    for j = 1:401
        deltaVlowOrder(i,j) = hodographicShapingLowOrder( ( i - 1 ) * 401 + j, 4 );
    end
end

% Retrieve high-order solution departure dates
for i = 1:6 
    departureDateHighOrder(i,1) = hodographicShapingHighOrder(i,3);
end
% Retrieve high-order solution time-of-flights
for i = 1:21 
    TOFhighOrder(i,1) = hodographicShapingHighOrder( ( i - 1 ) * 6 + 1,2);
end
% Retrieve high-order solution deltaVs
for i = 1:21
    for j = 1:6
        deltaVhighOrder(i,j) = hodographicShapingHighOrder( ( i - 1 ) * 6 + j, 4 );
    end
end

% Retrieve departure dates for low-order solution, with number of revolutions set to 1, 
% for comparison with high order solution 
for i = 1:6 
    departureDateLowOrderOneRev(i,1) = hodographicShapingLowOrderOneRev(i,3);
end
% Retrieve time-of-flights for low-order solution, with number of revolutions set to 1
for i = 1:21 
    TOFlowOrderOneRev(i,1) = hodographicShapingLowOrderOneRev( ( i - 1 ) * 6 + 1,2);
end
% Retrieve deltaVs for low-order solution, with number of revolutions set to 1
for i = 1:21
    for j = 1:6
        deltaVlowOrderOneRev(i,j) = hodographicShapingLowOrderOneRev( ( i - 1 ) * 6 + j, 4 );
    end
end

%% Low-order solution

figure(1);
hold on;
contourf( departureDateLowOrder', TOFlowOrder', deltaVlowOrder / 1e3, 3000, 'LineColor', 'None' );
grid;
hold off;
ylabel('TOF [days]', 'fontSize', 11);
xlabel('Departure date [MJD2000]', 'fontSize', 11);
colormap(manuallyDefinedColorbar);
c = colorbar;
caxis([0 45]);
set(c,'YTick',[7 8 10 15 20 40 45])
title(c,'$\Delta V$ [km/s]', 'fontSize', 10, 'Interpreter', 'Latex');



%% High-order solution

figure(2);
subplot(1,2,1);
hold on;
contourf( departureDateHighOrder', TOFhighOrder', deltaVhighOrder / 1e3, 3000, 'LineColor', 'None' );
plot( hodographicShapingHighOrder(:,3), hodographicShapingHighOrder(:,2), '.k', 'markerSize', 10);
grid;
hold off;
ylabel('TOF [days]', 'fontSize', 11);
xlabel('Departure date [MJD2000]', 'fontSize', 11);

colormap(manuallyDefinedColorbar);
c = colorbar;
caxis([0 45]);
set(c,'YTick',[7 8 10 15 20 40 45])
title(c,'$\Delta V$ [km/s]', 'fontSize', 10, 'Interpreter', 'Latex');
title('High-order solution');

figure(2);
subplot(1,2,2);
hold on;
contourf( departureDateLowOrderOneRev', TOFlowOrderOneRev', deltaVlowOrderOneRev / 1e3, 3000, 'LineColor', 'None' );
plot( hodographicShapingLowOrderOneRev(:,3), hodographicShapingLowOrderOneRev(:,2), '.k', 'markerSize', 10);
grid;
hold off;
ylabel('TOF [days]', 'fontSize', 11);
xlabel('Departure date [MJD2000]', 'fontSize', 11);
title('Low-order solution');

colormap(manuallyDefinedColorbar);
c = colorbar;
caxis([0 45]);
set(c,'YTick',[7 8 10 15 20 40 45])
xlim([departureDateHighOrder(1) departureDateHighOrder(end)]);
ylim([TOFhighOrder(1) TOFhighOrder(end)])
title(c,'$\Delta V$ [km/s]', 'fontSize', 10, 'Interpreter', 'Latex');





