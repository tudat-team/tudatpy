clc
close all

directory = '/home/dominic/.tudat/resource/';

for i=1:6

    residuals = load(strcat(directory,'ifms_residuals_utc_',num2str(i-1),'.dat'));
    times = load(strcat(directory,'ifms_times_utc_',num2str(i-1),'.dat'));

    if( i == 1 )
        first_time = times( 1 );
    end

    figure(4)

    scatter( times - first_time, residuals )
    hold on
    grid on
    
        figure(5)
    residuals_tdb = load(strcat(directory,'ifms_residuals_',num2str(i-1),'.dat'));
    times_tdb = load(strcat(directory,'ifms_times_',num2str(i-1),'.dat'));
    scatter( times_tdb - first_time, residuals_tdb, '*' )
        hold on
    grid on
    ylim([-0.1 0.1])

end    
