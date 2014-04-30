%% Kalman Filter
clear all; clc;
%% A priori statistics
psd = 0.1; % - PSD (power spectral density) of the random acceleration:
semc = 3; % - Standard error of measured coordinates:
semav = 0.5;% - Standard error of measured abs. velocity:
seiv = 0.01;% - Standard error of initial velocity:
seic = 0.1;% - Standard error of initial coordinates:
% Measured coordinates: