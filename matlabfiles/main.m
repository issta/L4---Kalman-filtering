%% Kalman Filter
clear all;clc;close all
addpath(['/Users/kevin/SkyDrive/KTH Work/',...
    'Period 4 2014/GNSS/Labs/L4 ',...
    '- Kalman filtering/']);
%% Import measured values from sheet 1
data.s1 = xlsread('Measurement.xlsx');
meas.time = data.s1(1:25,1);
meas.east = data.s1(1:25,2);
meas.north = data.s1(1:25,3);
meas.speed = data.s1(1:25,4);
% Import true values from sheet 2
data.s2 = xlsread('Measurement.xlsx','True values');
true.time = data.s2(1:25,1);
true.east = data.s2(1:25,2);
true.north = data.s2(1:25,3);
true.vel_east = data.s2(1:25,4);
true.vel_north = data.s2(1:25,5);
%% A priori statistics
PSD = 0.01; % - PSD (power spectral density) of the
% random acceleration
meas.sd_coord = 3; % m - Standard error of measured coordinates
meas.sd_abs_vel = 0.5; % m/s - Standard error of measured
% abs. velocity
sd_ini_vel = 3; % m/s - Standard error of initial velocity
sd_ini_coord = 10; % m - Standard error of initial coordinates
ve = 3.53; % m/s
vn = 0.86; % m/s
dt = 2; % time difference -> 2 sec between measurements
%% State variables
xk = [meas.east(1) meas.north(1) ve vn];
xk = padarray(xk,24,0,'post');
xk = xk';

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
Qx = cov(xk(:,1));
%% FOR LOOP
%
%% Equation 25
vm = sqrt(ve^2 + vn^2); % should be equal to speed_meas(1)
%% Equation 26
Hk = [1,0, 0,        0;...
    0, 1, 0,        0;...
    0, 0, ve/vm,  vn/vm];
% Lk = Hk * xk;
Lk = [meas.east(1) meas.north(1)...
    sqrt((xk(3,1))^2 + (xk(4,1))^2)]';
Rk = cov(Lk);
%%
for i=1:24
    xk(:,i+1) = Tk * xk(:,i); % Equation 16 Time propagation
    Qx = Tk * Qx * Tk' + Qk; % Equation 16
    Kk = Qx * Hk'*inv([Rk + Hk * Qx * Hk']); % Equation 17 Gain
    xk(:,i+1) = xk(:,i) + Kk*[ Lk - Hk * xk(:,i) ]; % Equation18 
%     Measurement update
    % Equation 19
    Qx = [eye(length(Kk*Hk))-Kk*Hk]*Qx;
    % Equation 22
    Lk = [meas.east(i) meas.north(i)...
        sqrt((xk(3,i))^2 + (xk(4,i))^2)]';
    Hk = inv(xk(:,i))*Lk;
    final.xplot(:,i+1) = xk(:,i);
end
%% Smoothing

%% Plot
final.x1 = xk(1,:)'; % final values
final.y1 = xk(2,:)'; % final values
meas.x2 = meas.east; % original
meas.y2 = meas.north; % original
true.x3 = true.east; % true
true.y3 = true.north; % true
final.plot = plot(final.x1,final.y1,'.','color','b');
hold on;
originalplot = plot(meas.x2,meas.y2,'.','color','r');
true.plot = plot(true.x3,true.y3,'-');
legend([final.plot,originalplot,true.plot],...
    'Final','Original','True',...
    'location','best');
hold off;

% hold on plot(east_meas,north_meas,'r')
% plot(east_true,north_true,'b')