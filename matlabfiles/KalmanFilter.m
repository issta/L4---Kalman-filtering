% Kalman filtering
% Created by Dominik Juchli, juchli@kth.se, 30-04-2014
function [x_s, e_p, n_p] = KalmanFilter
true = load('truevalues.mat');      % Load true values
meas = load('measurements.mat');    % Load measurements
len = size(meas.data,1);                 % number of measurements

%% Allocation and Initialisation

% Initial state
e_0 = meas.data(1,2);       % easting position [m]
n_0 = meas.data(1,3);       % northing position [m]
v_e_0 = 3.53;               % easting speed [m/s]
v_n_0 = 0.86;               % northing speed [m/s]
q_e = 0.01;                 % power spectral density
q_n = 0.01;

dt = 2;                     % time step between measurements [s]
x = zeros(4,len);           % State vector
F = [0 0 1 0;               % System dynamics
     0 0 0 1;
     0 0 0 0;
     0 0 0 0];
G = [0 0;                   % Distribution matrix
     0 0;
     1 0;
     0 1];
T = [1 0 dt 0;              % Transition matrix
     0 1 0 dt;
     0 0 1  0;
     0 0 0  1];
Q_k = [q_e*dt^3/3 0 q_e*dt^2/2 0;   % process noise covariance
       0 q_n*dt^3/3 0 q_n*dt^2/2;
       q_e*dt^2/2 0 q_e*dt 0;
       0 q_n*dt^2 0 q_n*dt];

% Standard devation
sigma_e = 3;                % easting position measurement [m]
sigma_n = 3;                % northing position measurement [m]
sigma_nv = 0.5;             % speed [m/s]
sigma_e0 = 10;              % initial easting position [m]
sigma_n0 = 10;              % initial northing position [m]
sigma_ve0 = 3;              % initial easting speed [m/s]
sigma_vn0 = 3;              % initial northing speed [m/s]
R = diag([sigma_e sigma_n sigma_nv]);   % Covariance matrix of measurements

% Prediction
x_p = zeros(4,len-1);
e_p = zeros(1,len-1);
n_p = zeros(1,len-1);
v_e_p = zeros(1,len-1);
v_n_p = zeros(1,len-1);
v_p = zeros(1,len-1);
Q_p = zeros(4,4,len-1);
Q_x = zeros(4,4,len);
H = zeros(3,4,len-1);
K = zeros(4,3,len-1);
L = zeros(3,len-1);
h = zeros(3,len-1);


 %% Kalman filtering
 
 % 1. Initialisation:
 x(:,1) = [e_0; n_0; v_e_0; v_n_0];      % initial state
 Q_x(:,:,1) = diag([sigma_e0^2 ...    % covariance matrix of initial state
   sigma_n0^2 sigma_ve0^2 sigma_vn0^2]);

for i=1:len-1
    
     % 2. Time propagation
     x_p(:,i) = T*x(:,i);                     % prediction
     e_p(i) = x_p(1,i);
     n_p(i) = x_p(2,i);
     v_e_p(i) = x_p(3,i);
     v_n_p(i) = x_p(4,i);
     v_p(i) = sqrt(v_e_p(i)^2+v_n_p(i)^2);
     Q_p(:,:,i) = T*Q_x(:,:,i)*T'+Q_k;          % predicted covariance
     
     % 3. Gain calculation
     H(:,:,i) = [1 0 0 0;                          % design matrix
                 0 1 0 0;
          0 0 v_e_p(i)/v_p(i) v_n_p(i)/v_p(i)];
      
     K(:,:,i) = Q_p(:,:,i)*H(:,:,1)'*inv(R ...   % Kalman gain
         + H(:,:,i)*Q_p(:,:,i)*H(:,:,i)');   
     
     % 4. Measurement update
     L(:,i) = meas.data(i+1,2:4)';          % vector of observation 
     h(:,i) = [e_p(i); n_p(i); v_p(i)];
     x(:,i+1) = x_p(:,i) + K(:,:,i)*(L(:,i)-h(:,i));
     
     % 5. Covariance update
     Q_x(:,:,i+1) = (eye(4) - K(:,:,i)*H(:,:,i))*Q_p(:,:,i);
end

%% Smoothing

% Allocation
x_s = zeros(4,len);
Q_x_s = zeros(4,4,len);

% Initialisation
x_s(:,len) = x(:,len);
Q_x_s(:,:,len) = Q_x(:,:,len);

for j=len:-1:2
    D = Q_x(:,:,j)*T'*inv(Q_p(:,:,j-1));
    a(:,j) = x_s(:,j) - x_p(:,j-1);
    x_s(:,j-1) = x(:,j-1) + D*(a(:,j));
    Q_x_s(:,:,j-1) = Q_x(:,:,j) + D*(Q_x_s(:,:,j) - Q_p(:,:,j-1))*D';
end

%% Visualization

% figure(1)
% plot(e_p,n_p,'b')
% hold on
% plot(true.data(:,2),true.data(:,3))
% hold on
% plot(x_s(1,:),x_s(2,:),'g')
% hold on
% plot(meas.data(:,2),meas.data(:,3),'.')