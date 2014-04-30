%% Kalman Filter
clear all; clc;
addpath(['/Users/kevin/SkyDrive/KTH Work/Period 4 2014/GNSS/Labs/L4 ',...
'- Kalman filtering/']);
% Import measured values from sheet 1
data_meas = xlsread('Measurement.xlsx');
time_meas = data_meas(1:25,1);
east_meas = data_meas(1:25,2);
north_meas = data_meas(1:25,3);
speed_meas = data_meas(1:25,4);
% Import true values from sheet 2
data_true = xlsread('Measurement.xlsx','True values');
time_true = data_true(1:25,1);
east_true = data_true(1:25,2);
north_true = data_true(1:25,3);
vel_east_true = data_true(1:25,4);
vel_north_true = data_true(1:25,5);
%% A priori statistics
PSD = 0.01; % - PSD (power spectral density) of the random acceleration
sd_meas_coord = 3; % m - Standard error of measured coordinates
sd_meas_abs_vel = 0.5; % m/s - Standard error of measured abs. velocity
sd_ini_vel = 3; % m/s - Standard error of initial velocity
sd_ini_coord = 10; % m - Standard error of initial coordinates
vel_east = 3.53; % m/s
vel_north = 0.86; % m/s
%% 
