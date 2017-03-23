close all;
clear all;

all_data = load('C:\Users\uidg5371\udacity\CarND-Unscented-Kalman-Filter-Project\vs\unscented_kafi\unscented_kafi\output.txt');
laser_data = all_data(all_data(:,10)==1,:);
radar_data = all_data(all_data(:,10)==2,:);

x_est = all_data(:,1);
y_est = all_data(:,2);
radar_x = radar_data(:,6);
radar_y = radar_data(:,7);
laser_x = laser_data(:,6);
laser_y = laser_data(:,7);

count = 1;
x_gt = zeros(length(x_est),1);
y_gt = zeros(length(y_est),1);
vx_gt = zeros(length(x_est),1);
vy_gt = zeros(length(y_est),1);
radar_meas = zeros(length(radar_data),3);
laser_meas = zeros(length(laser_data),2);
NIS_Radar = radar_data(:,8);
NIS_Laser = laser_data(:,8);

fid = fopen('C:\Users\uidg5371\udacity\CarND-Unscented-Kalman-Filter-Project\data\sample-laser-radar-measurement-data-2.txt');
tline = fgets(fid); % read the next line of the data file
while ischar(tline)  % go through lines of data file
    if tline(1) == 'L' % laser measurement
        line_vector = textscan(tline,'%s %f %f %f %f %f %f %f %f %f');
        x_gt(count) = line_vector{5};
        y_gt(count) = line_vector{6};
        vx_gt(count) = line_vector{7};
        vy_gt(count) = line_vector{8};
        laser_meas(count,:) = [line_vector{2},line_vector{3}];
    elseif tline(1) == 'R' % radar measurement
        line_vector = textscan(tline,'%s %f %f %f %f %f %f %f %f %f %f');
        x_gt(count) = line_vector{6};
        y_gt(count) = line_vector{7};
        vx_gt(count) = line_vector{8};
        vy_gt(count) = line_vector{9};
        radar_meas(count,:) = [line_vector{2},line_vector{3},line_vector{4}];
    else
        disp('Error: not laser nor radar')
        return;
    end
    
    tline = fgets(fid); % read the next line of the data file
    count = count + 1;
end

% Calc RMSE
rmse_x = sqrt(mean((x_est-x_gt).^2))
rmse_y = sqrt(mean((y_est-y_gt).^2))
Rx_radar = sqrt(var(radar_x-x_gt(all_data(:,10)==2)))
Rx_laser = sqrt(var(laser_x-x_gt(all_data(:,10)==1)))
Ry_radar = sqrt(var(radar_y-y_gt(all_data(:,10)==2)))
Ry_laser = sqrt(var(laser_y-y_gt(all_data(:,10)==1)))

figure(1)
hold on
plot(radar_data(:,9) ,radar_data(:,8), '-r'); 
plot(radar_data(:,9), ones(size(radar_data(:,9))).*7.815, '-k');
title('NIS Radar')
hold off

figure(2)
hold on
plot(laser_data(:,9),laser_data(:,8), '-b'); 
plot(laser_data(:,9), ones(size(laser_data(:,9))).*5.991, '-k');
title('NIS Laser')
hold off

figure(3)
hold on;
plot(x_gt, y_gt, '-og');
plot(x_est, y_est, '-or');
plot(radar_x, radar_y, '-xb');
plot(laser_x, laser_y, '-xy');
axis equal
legend('GT', 'est', 'Radar', 'Laser')
hold off

