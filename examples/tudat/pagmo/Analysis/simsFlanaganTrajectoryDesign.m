%
% This script processes the results of the hodographic shaping trajectory design optimization
% run by the simsFlanaganTrajectoryExample.cpp Tudat/Pagmo2 example.
%

clear all;
clc;
close all;

% Load Sims-Flanagan results (trajectory, mass, thrust, and thrust acceleration profiles)
SimsFlanaganTrajectory = importdata('../SimulationOutput/SimsFlanaganTrajectory.dat');
SimsFlanaganMassProfile = importdata('../SimulationOutput/SimsFlanaganMassProfile.dat');
SimsFlanaganThrustProfile = importdata('../SimulationOutput/SimsFlanaganThrustProfile.dat');
SimsFlanaganThrustAccelerationProfile = importdata('../SimulationOutput/SimsFlanaganThrustAcceleration.dat');

% Load hodographic shaping results (trajectory, mass, thrust, and thrust acceleration profiles)
hodographicShapingTrajectory = importdata('../SimulationOutput/hodographicShapingTrajectory.dat');
hodographicShapingMassProfile = importdata('../SimulationOutput/hodographicShapingMassProfile.dat');
hodographicShapingThrustProfile = importdata('../SimulationOutput/hodographicShapingThrustProfile.dat');
hodographicShapingThrustAccelerationProfile = importdata('../SimulationOutput/hodographicShapingThrustAcceleration.dat');

% Retrieve initial guess used for Sims-Flanagan, as an approximation for
% the shape-based solution
initialGuessThrustProfile = importdata('../SimulationOutput/initialGuessThrustProfile.dat');
initialGuessThrustAccelerationProfile = importdata('../SimulationOutput/initialGuessThrustAccelerationProfile.dat');
initialGuessMassProfile = importdata('../SimulationOutput/initialGuessMassProfile.dat');

% Compute norm of the thrust and thrust acceleration vectors for hodographic shaping results
for i = 1:size(hodographicShapingTrajectory,1)
    hodographicShapingThrustProfile(i,5) = norm( hodographicShapingThrustProfile(i,2:4) );
    hodographicShapingThrustAccelerationProfile(i,5) = norm( hodographicShapingThrustAccelerationProfile(i,2:4) );
end

% Compute norm of the thrust and thrust acceleration vectors for Sims-Flanagan results
for i = 1:size(SimsFlanaganTrajectory,1)
    SimsFlanaganThrustProfile(i,5) = norm( SimsFlanaganThrustProfile(i,2:4) );
    SimsFlanaganThrustAccelerationProfile(i,5) = norm( SimsFlanaganThrustAccelerationProfile(i,2:4) );
end

% Compute norm of the thrust and thrust acceleration vectors for Sims-Flanagan initial guess 
for i = 1:size(initialGuessThrustProfile,1)
    initialGuessThrustProfile(i,5) = norm( initialGuessThrustProfile(i,2:4) );
    initialGuessThrustAccelerationProfile(i,5) = norm( initialGuessThrustAccelerationProfile(i,2:4) );
end


%% Plots

figure(1);

subplot(3,4,1);
% Plot the hodographic shaping trajectory, and the associated approximate initial guess for Sims-Flanagan
hold on;
plot3( hodographicShapingTrajectory(:,2), hodographicShapingTrajectory(:,3), hodographicShapingTrajectory(:,4), ...
    'b', 'lineWidth', 1.5 );
plot3( hodographicShapingTrajectory(:,2), hodographicShapingTrajectory(:,3), hodographicShapingTrajectory(:,4), ...
    'k', 'lineWidth', 0.8 );
hold off;
grid;
xlabel('x [m]', 'Interpreter', 'Latex');
ylabel('y [m]', 'Interpreter', 'Latex');
zlabel('z [m]', 'Interpreter', 'Latex');
title('Trajectory', 'fontSize', 12);
view(35,25);

subplot(3,4,2);
% Plot the hodographic shaping thrust profile, and the associated approximate initial guess for Sims-Flanagan
hold on;
plot( hodographicShapingThrustProfile(:,1), hodographicShapingThrustProfile(:,5), 'b', 'lineWidth', 1.5 );
plot( initialGuessThrustProfile(:,1), initialGuessThrustProfile(:,5), 'k', 'lineWidth', 0.8 );
hold off;
grid;
ylabel('Thrust [N]', 'Interpreter', 'Latex');
xlabel('Time [s]', 'Interpreter', 'Latex');
title('Thrust profile', 'fontSize', 12);

subplot(3,4,3);
% Plot the hodographic shaping thrust acceleration profile, and the associated approximate initial guess 
% for Sims-Flanagan
hold on;
plot( hodographicShapingThrustAccelerationProfile(:,1), hodographicShapingThrustAccelerationProfile(:,5), ...
    'b', 'lineWidth', 1.5 );
plot( initialGuessThrustAccelerationProfile(:,1), initialGuessThrustAccelerationProfile(:,5), 'k', ...
    'lineWidth', 0.8 );
hold off;
grid;
ylabel('Thrust acceleration [$\mathrm{m/s^2}$]', 'Interpreter', 'Latex');
xlabel('Time [s]', 'Interpreter', 'Latex');
title('Acceleration profile', 'fontSize', 12);

subplot(3,4,4);
% Plot the hodographic shaping mass profile, and the associated approximate initial guess for Sims-Flanagan
hold on;
plot( hodographicShapingMassProfile(:,1), hodographicShapingMassProfile(:,2), 'b', 'lineWidth', 1.5 );
plot( initialGuessMassProfile(:,1), initialGuessMassProfile(:,2), 'k', 'lineWidth', 0.8 );
hold off;
grid;
ylabel('Mass [kg]', 'Interpreter', 'Latex');
xlabel('Time [s]', 'Interpreter', 'Latex');
title('Mass profile', 'fontSize', 12);
legend({'Hodographic shaping', 'Initial guess for Sims-Flanagan'}, 'fontSize', 10, 'location', 'Best');

subplot(3,4,5);
% Plot the Sims-Flanagan trajectory
hold on;
plot3( SimsFlanaganTrajectory(:,2), SimsFlanaganTrajectory(:,3), SimsFlanaganTrajectory(:,4), ...
    'r', 'lineWidth', 1.5 );
hold off;
grid;
xlabel('x [m]', 'Interpreter', 'Latex');
ylabel('y [m]', 'Interpreter', 'Latex');
zlabel('z [m]', 'Interpreter', 'Latex');
view(35,25);

subplot(3,4,6);
% Plot the Sims-Flanagan thrust profile, and the one from hodographic shaping (that was used to derive the 
% initial guess)
hold on;
plot( SimsFlanaganThrustProfile(:,1), SimsFlanaganThrustProfile(:,5), 'r', 'lineWidth', 1.5 );
plot( hodographicShapingThrustProfile(:,1), hodographicShapingThrustProfile(:,5), 'b', 'lineWidth', 1.5);
hold off;
grid;
ylabel('Thrust [N]', 'Interpreter', 'Latex');
xlabel('Time [s]', 'Interpreter', 'Latex');

subplot(3,4,7);
hold on;
% Plot the Sims-Flanagan thrust acceleration profile, and the one from hodographic shaping 
% (that was used to derive the initial guess)
plot( SimsFlanaganThrustAccelerationProfile(:,1), SimsFlanaganThrustAccelerationProfile(:,5), ...
    'r', 'lineWidth', 1.5 );
plot( hodographicShapingThrustAccelerationProfile(:,1), hodographicShapingThrustAccelerationProfile(:,5), ...
    'b', 'lineWidth', 1.5);
hold off;
grid;
ylabel('Thrust acceleration [$\mathrm{m/s^2}$]', 'Interpreter', 'Latex');
xlabel('Time [s]', 'Interpreter', 'Latex');

subplot(3,4,8);
hold on;
% Plot the Sims-Flanagan mass profile, and the one from hodographic shaping (that was used to derive the 
% initial guess)
plot( SimsFlanaganMassProfile(:,1), SimsFlanaganMassProfile(:,2), 'r', 'lineWidth', 1.5 );
plot( hodographicShapingMassProfile(:,1), hodographicShapingMassProfile(:,2), 'b', 'lineWidth', 1.5 );
hold off;
grid;
ylabel('Mass [kg]', 'Interpreter', 'Latex');
xlabel('Time [s]', 'Interpreter', 'Latex');
legend({'Sims-Flanagan', 'Hodographic shaping'}, 'fontSize', 10, 'location', 'Best');

subplot(3,4,9);
% Plot the difference between the Sims-Flanagan and the hodographic trajectories
hold on;
plot3( hodographicShapingTrajectory(:,2) - SimsFlanaganTrajectory(:,2), ...
    hodographicShapingTrajectory(:,3) - SimsFlanaganTrajectory(:,3), ...
    hodographicShapingTrajectory(:,4) - SimsFlanaganTrajectory(:,4), 'k', 'lineWidth', 1.2 );
hold off;
grid;
xlabel('x [m]', 'Interpreter', 'Latex');
ylabel('y [m]', 'Interpreter', 'Latex');
zlabel('z [m]', 'Interpreter', 'Latex');
view(35,25);

subplot(3,4,10);
% Plot the difference between the Sims-Flanagan and the hodographic thrust profiles
hold on;
plot( hodographicShapingThrustProfile(:,1), ...
    hodographicShapingThrustProfile(:,5) - SimsFlanaganThrustProfile(:,5), 'k', 'lineWidth', 1.2 );
hold off;
grid;
ylabel('Thrust [N]', 'Interpreter', 'Latex');
xlabel('Time [s]', 'Interpreter', 'Latex');

subplot(3,4,11);
% Plot the difference between the Sims-Flanagan and the hodographic thrust acceleration profiles
hold on;
plot( hodographicShapingThrustAccelerationProfile(:,1),...
    hodographicShapingThrustAccelerationProfile(:,5) - SimsFlanaganThrustAccelerationProfile(:,5), ...
    'k', 'lineWidth', 1.2 );
hold off;
grid;
ylabel('Thrust acceleration [$\mathrm{m/s^2}$]', 'Interpreter', 'Latex');
xlabel('Time [s]', 'Interpreter', 'Latex');

subplot(3,4,12);
% Plot the difference between the Sims-Flanagan and the hodographic mass profiles
hold on;
plot( hodographicShapingMassProfile(:,1), hodographicShapingMassProfile(:,2) - SimsFlanaganMassProfile(:,2), ...
    'k', 'lineWidth', 1.2 );
hold off;
grid;
ylabel('Mass [kg]', 'Interpreter', 'Latex');
xlabel('Time [s]', 'Interpreter', 'Latex');
legend({'Diff shaping-SimsFlanagan'}, 'fontSize', 10, 'location', 'Best');

