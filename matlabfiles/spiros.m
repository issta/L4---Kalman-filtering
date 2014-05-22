clear all;
clc;


TR= [0	0.00	0.00	4.33	2.50    % true values
    2	8.61	4.88	4.36	2.50
    4	17.34	10.03	4.32	2.55
    6	26.09	15.29	4.49	2.76
    8	35.08	20.55	4.59	2.58
    10	44.49	25.67	5.00	2.55
    12	54.72	30.82	5.14	2.61
    14	65.03	36.27	5.07	2.70
    16	75.15	41.81	5.08	2.80
    18	85.25	47.19	5.00	2.62
    20	95.08	52.69	4.86	2.78
    22	104.79	58.25	4.92	2.80
    24	114.90	63.75	4.99	2.66
    26	124.77	69.00	4.91	2.50
    28	134.50	74.13	4.90	2.56
    30	144.56	79.10	5.07	2.39
    32	154.86	83.70	5.15	2.26
    34	165.16	88.10	5.32	2.20
    36	175.38	92.61	5.13	2.23
    38	185.72	97.40	5.24	2.37
    40	196.07	102.12	5.13	2.40
    42	206.45	107.00	5.24	2.42
    44	216.89	111.60	5.12	2.26
    46	227.04	116.19	5.06	2.21
    48	237.12	120.52	5.11	2.16];

MRS=[0	-9.82	0.06	3.63            % measured values
    2	10.93	5.29	5.2
    4	16.29	9.71	4.98
    6	21.91	22.1	5.45
    8	34.44	24.68	5.7
    10	46.5	25.57	5.34
    12	52.85	32.46	5.7
    14	67.55	34.22	5.12
    16	78.44	41.59	6.44
    18	88.1	46.25	5.26
    20	99.6	48.8	5.7
    22	105.85	59.17	5.55
    24	116.73	67.02	5.44
    26	120.08	68.29	5.03
    28	135.96	70.02	5.4
    30	146.93	85.7	5.6
    32	159.62	84.83	5.12
    34	166.52	90.8	5.57
    36	178.12	93.98	5.09
    38	189.23	92.41	5.87
    40	192.56	102.25	5.33
    42	207.93	103.3	4.96
    44	216.83	110.51	5.59
    46	224.44	119.15	5.81
    48	236.97	120.5	5.32];


% standard deviations (given from exercise)

SDCE = 3;                % Standard deviation measured coordinates
SDCN= 3;                % Standard deviation measured coordinates
SDabsV = 0.5;             %>> >> absolute velocity
SDIvelE = 3;              % >>>> initial velocity EAST
SDIvelN = 3;              % >>>> initial velocity NORTH
SDICE = 10;              % >> >> initial coordinates EAST
SDICN = 10;              % >>>> initial coordinates NORTH

%Initial conditions
East0=-9.82;        %Initial east position
North0= 0.06;     % Initial north ponition
Veast0=3.53;         % Velocity horizontal
Vnorth0=0.86;       % Velocity perpendicular
qe=0.1;                 %PSD
qn=0.1;

%______________



G=[0 0          %distribution matrix
    0 0
    1 0
    0 1];

F=[0 0 1 0;               % system dynamic matrix
    0 0 0 1;
    0 0 0 0;
    0 0 0 0];


x=zeros(4,25);   % initial state,state vector

R=[SDCE,0,0
    0,SDCN,0
    0,0,SDabsV];

%  predicted

xPR=  zeros(4,24);
EASTPR = zeros(1,24);
NORTHPR = zeros(1,24);

VELeastPR = zeros(1,24);
VELnorthPR = zeros(1,24);
VELp = zeros(1,24);

Qp = zeros(4,4,24);
Qx = zeros(4,4,25);

H = zeros(3,4,24);    %(26)
K = zeros(4,3,24);    %(17)
L = zeros(3,24);    %(23)
h = zeros(3,24);    %(27)

% ______Kahlman filter loop___________

% 1 part
%define T and Qk
T= [1 0 2 0     % trasnition matrix where 2 is measurement step
    0 1 0 2
    0 0 1 0
    0 0 0 1];


Qk=[ (qe*8/3) 0 (qe*4/2) 0      %white noise covariance
    0 (qn*8/3) 0 (qn*4/2)
    (qe*4/2) 0 (qe*2) 0
    0 (qn*4/2) 0 (qn*2)];

%estimate x, Qx

Qx(:,:,1)=diag([SDICE*SDICE SDICN*SDICN SDIvelE*SDIvelE SDIvelN*SDIvelN]);

x(:,1)=[East0 North0 Veast0 Vnorth0];

E= [1     0     0     0
    0     1     0     0
    0     0     1     0
    0     0     0     1];

for k=1:24
    %______________
    xPR(:,k)=T*x(:,k);
    
    EASTPR(k) = xPR(1,k);
    NORTHPR(k) = xPR(2,k);
    
    VELeastPR(k) = xPR(3,k);
    VELnorthPR(k) = xPR(4,k);
    
    VELp(k) = sqrt( VELeastPR(k)^2+VELnorthPR(k)^2);
    
    Qp(:,:,k) = T*Qx(:,:,k)*T'+Qk;
    
    %_______________________
    
    H(:,:,k)=[1 0 0 0
        0 1 0 0
        0 0 VELeastPR(k)/VELp(k) VELnorthPR(k)/VELp(k)];
    
    K(:,:,k)= Qp(:,:,k)*H(:,:,1)'*inv(R+H(:,:,k)*Qp(:,:,k)*H(:,:,k)');
    
    %______________
    
    L(:,k) = MRS(k+1,2:4)'
    
    h(:,k)=[EASTPR(k),NORTHPR(k),VELp(k)];
    
    x(:,k+1)=xPR(:,k)+K(:,:,k)*(L(:,k)-h(:,k));
    
    %______________
    
    Qx(:,:,k+1)= (E - K(:,:,k)*H(:,:,k))*Qp(:,:,k);
    
end


%%__________
Xs=zeros(4,25);
QXs=zeros(4,4,25);

Xs(:,25) = x(:,25);
QXs(:,:,25) = Qx(:,:,25);

%

for o=25:-1:2;
    
    W = Qx(:,:,o)*T'*inv(Qp(:,:,o-1));
    c(:,o) = Xs(:,o) - xPR(:,o-1);
    
    Xs(:,o-1) = x(:,o-1) + W*(c(:,o));
    
    QXs(:,:,o-1) = Qx(:,:,o) + W*(QXs(:,:,o) - Qp(:,:,o-1))*W';
    
end

figure(2)
plot(...
    EASTPR,NORTHPR,'.',...  % 1
    TR(:,2),TR(:,3),'-',... % 2
    Xs(1,:),Xs(2,:),'.',... % 3
    MRS(:,2),MRS(:,3),'.'); % 4
legend('1','2','3','4')