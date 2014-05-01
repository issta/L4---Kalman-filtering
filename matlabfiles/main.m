%% Kalman Filter
clear all;clc;close all
addpath(['/Users/kevin/SkyDrive/KTH Work/',...
    'Period 4 2014/GNSS/Labs/L4 ',...
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
PSD = 0.01; % - PSD (power spectral density) of the
% random acceleration
sd_meas_coord = 3; % m - Standard error of measured coordinates
sd_meas_abs_vel = 0.5; % m/s - Standard error of measured
% abs. velocity
sd_ini_vel = 3; % m/s - Standard error of initial velocity
sd_ini_coord = 10; % m - Standard error of initial coordinates
ve = 3.53; % m/s
vn = 0.86; % m/s
dt = 2; % time difference -> 2 sec between measurements
%% State variables
xk = [east_meas(1) north_meas(1) ve vn]';

% Equation 4
F = zeros(4);
F(1,3) = 1;
F(2,4) = 1;
G = zeros(4,2);
G(3,1) = 1;
G(4,2) = 1;
% u = [acc_east; acc_north];
%% Equation 5
Tk = eye(length(F)) + dt * F;
%% Equation 9
Q = [ PSD 0 ; 0 PSD];
%% Equation 11
QG = G*Q*G';
%% Equation 12
Qk = QG * dt + (F*QG + QG*F')*dt^2/2 + F*QG*F'*dt^3/3;
%% Equation 15
Qx(1) = cov(xk);
%% FOR LOOP
%
%% Equation 25
vm = sqrt(ve^2 + vn^2); % should be equal to speed_meas(1)
%% Equation 26
Hk = [1,0, 0,        0;...
    0, 1, 0,        0;...
    0, 0, ve/vm,  vn/vm];
Lk = Hk * xk;
Rk = cov(Lk);
for i=1:25
    %% Equation 16
    % Time propagation
    xk(:,i+1) = Tk * xk(:,i);
    Qx = Tk * Qx * Tk' + Qk;
    %% Equation 17
    % Gain calculation
    Kk = Qx * Hk'/(Rk + Hk * Qx * Hk');
    %% Equation 18
    % Measurement update
    xk(:,i+1) = xk(:,i) + Kk*[Lk-Hk*xk(:,i)];
    %% Equation 19
    Qx = [eye(length(Kk*Hk))-Kk*Hk]*Qx;
    %% Equation 22
    Lk = [east_meas(i) north_meas(i) sqrt(ve^2 + vn^2)]';
    xplot(:,i+1) = xk(:,i);
end
x1 = xk(1,:); % final values
y1 = xk(2,:); % final values
x2 = east_meas;
y2 = north_meas;
x3 = east_true;
y3 = north_true;
finalplot = plot(x1,y1,'.','color','r');
hold on;
originalplot = plot(x2,y2,'-.','color','r');
trueplot = plot(x3,y3,'.');
legend([finalplot,originalplot,trueplot],...
    'Final','Original','True',...
    'location','best');
hold off;

% hold on plot(east_meas,north_meas,'r')
% plot(east_true,north_true,'b')