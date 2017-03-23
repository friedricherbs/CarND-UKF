close all;
clear all;

all_data = load('C:\Users\uidg5371\udacity\CarND-Unscented-Kalman-Filter-Project\vs\unscented_kafi\unscented_kafi\output.txt');
laser_data = all_data(all_data(:,10)==1,:);
radar_data = all_data(all_data(:,10)==2,:);

x_est = all_data(:,1);
y_est = all_data(:,2);
meas_x = all_data(:,6);
meas_y = all_data(:,7);

gt = load('C:\Users\uidg5371\udacity\CarND-Unscented-Kalman-Filter-Project\data\obj_pose-laser-radar-synthetic-gt.txt');
x_gt = gt(:,1);
y_gt = gt(:,2);
v_gt = gt(:,3);
yaw_gt = gt(:,4);
yaw_rate_gt = gt(:,5);
vx_gt = gt(:,6);
vy_gt = gt(:,7);
time_gt = gt(:,8);

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
 plot(meas_x, meas_y, '-xb');
 axis equal
 legend('GT', 'est', 'Meas')
 hold off

