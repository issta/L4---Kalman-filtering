%% Kalman Filter
clear all; clc;
addpath(['/Users/kevin/SkyDrive/KTH Work/Period 4 2014/GNSS/Labs/L4 ',...
'- Kalman filtering/']);
data = xlsread('Measurement.xlsx');
time = data(1:25,1);
east = data(1:25,2);
north = data(1:25,3);
speed = data(1:25,4);
%% A priori statistics
psd = 0.1; % - PSD (power spectral density) of the random acceleration:
semc = 3; % - Standard error of measured coordinates:
semav = 0.5;% - Standard error of measured abs. velocity:
seiv = 0.01;% - Standard error of initial velocity:
seic = 0.1;% - Standard error of initial coordinates:
%% 
