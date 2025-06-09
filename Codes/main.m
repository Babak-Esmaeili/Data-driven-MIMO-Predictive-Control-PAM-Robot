clc;
clear;
close all;

tic;
%% Time
ts = 0;
dt = 0.001;
tf = 10;

t = ts:dt:tf;

nSteps = length(t);

%% Important sizes
nInputs = 2;
nOutputs = nInputs;

N = 5;  % Prediction horizon
Nu = 4; % Control horizon

%% Controller parameters
eta = 1;
mu = 1;

m = 0.01*diag([1,1]);

epsilon = 10^-5;

pp = 3;
qq = 5;

c_11 = 1;
c_12 = 1;
c_1 = diag([c_11,c_12]);

c_21 = 0.3;
c_22 = 0.06;
c_2 = diag([c_21,c_22]);

c_s1 = 0.001;
c_s2 = 0.001;
c_s = diag([c_s1,c_s2]);

omega_1 = .1;
omega_2 = .1;
omega = [omega_1 omega_2]';

% Predictive
GAMMA = ones(N,1);
PAI = tril(repmat(-c_1*c_s,N,N));
LAMBDA = tril(repmat(-c_1,N,N));
% XI = diag(diag(repmat(c_1,N,N)));
XI = kron(eye(N),c_1);

K = zeros(nInputs,Nu*nInputs);
K(:,1:nInputs) = eye(nInputs);

gamma = 20;

%% Initialization
q = [0 0]';
qdot = [0 0]';

U = zeros(nInputs,nSteps);
Y = zeros(nOutputs,nSteps);
Y_d = zeros(nOutputs,nSteps);
E = zeros(nOutputs,nSteps);
E_I = zeros(nOutputs,nSteps);
S = zeros(nOutputs,nSteps);

dist = zeros(nOutputs,nSteps);
delta_d = zeros(nOutputs,nSteps);
delta_d_hat = zeros(nOutputs,nSteps);
delta_d_tilde = zeros(nOutputs,nSteps);

delta_U_eq = zeros(nInputs,nSteps);
delta_U_sw = zeros(nInputs,nSteps);
delta_U_s = zeros(nInputs,nSteps);

delta_U_mpc = zeros(nInputs,nSteps);

delta_U = zeros(nInputs,nSteps);

delta_Y = zeros(nOutputs,nSteps);

SAI_hat = cell(1,nSteps);
SAI_hat{1} = 1*[ 1  0
                 0  1 ];
          
% Predictive
S_hat_p = zeros(N*nOutputs,1);
delta_Up = zeros(Nu*nOutputs,1);
R_hat_p = zeros(N*nOutputs,1);
L_hat_p = zeros(N*nOutputs,1);
PHI_hat_p = zeros(N*nOutputs,Nu*nInputs);

%% Simulation
for k = 2:nSteps
   
%     Y_d(:,k+1) = [1.5*cos((2*pi/4)*t(k))
%                       1.2*cos((2*pi/7)*t(k))];

%         y_d(:,k+1) = [(2*pi/4)*cos((2*pi/4)*t(k))
%                        (2*pi/7)*cos((2*pi/7)*t(k)) ];
%         
%         y_d(:,k+1) = [sat_func_yd(y_d(1,k+1),-0.4,0.4)
%                       sat_func_yd(y_d(2,k+1),-0.6,0.6)];

%         y_d(:,k+1) = [(2*pi/10)*cos((2*pi/10)*t(k))
%                        (2*pi/7)*cos((2*pi/7)*t(k)) ];

%         y_d(:,k+1) = [sat_func_yd(y_d(1,k+1),-0.4,0.4)
%                       sat_func_yd(y_d(2,k+1),-0.6,0.6)];

%         Y_d(:,k+1) = [0.15*(t(k)<=2*tf/5) + 0.3*(t(k)>2*tf/5 & t(k)<=4*tf/5) + 0.15*(t(k)>4*tf/5 & t(k)<=5*tf/5)
%                       0.20*(t(k)<=tf/5) + 0.4*(t(k)>tf/5 & t(k)<=2.5*tf/5) + 0.20*(t(k)>2.5*tf/5 & t(k)<=4.5*tf/5) + 0.05*(t(k)>4.5*tf/5 & t(k)<=5*tf/5)];

%         Y_d(:,k+1) = [(7*pi/10)*cos((2*pi/12)*t(k))
%                        (5*pi/8)*cos((3*pi/10)*t(k))];
%         
%         Y_d(:,k+1) = [sat_func_yd(Y_d(1,k+1),-0.5,0.5)
%                       sat_func_yd(Y_d(2,k+1),-0.4,0.4)];

    Y_d(:,k+1) = [(7*pi/10)*cos((7*pi/10)*t(k))
                   (5*pi/8)*cos((5*pi/8)*t(k))];

    Y_d(:,k+1) = [sat_func_yd(Y_d(1,k+1),-0.5,0.5)
                  sat_func_yd(Y_d(2,k+1),-0.4,0.4)];
     
    %% Disturbance estimation
%     delta_d_hat(:,k) = delta_Y(:,k) - (SAI_hat{k-1}+m)*delta_U(:,k-1);
    delta_d_hat(:,k) = delta_d(:,k-1);
    
    %% PJM estimation
    SAI_hat{k} = SAI_hat{k-1} + (eta*(delta_Y(:,k)-delta_d_hat(:,k)-(SAI_hat{k-1}+m)*delta_U(:,k-1))*delta_U(:,k-1)')/(mu+norm(delta_U(:,k-1),2)^2);
    
    for i = 1:nInputs
        for j = 1:nOutputs
            if(sign(SAI_hat{k}(i,j)) ~= sign(SAI_hat{1}(i,j)))    % better!
%             if(abs(SAI_hat{k}(i,j)) <= epsilon)
                SAI_hat{k}(i,j) = SAI_hat{1}(i,j);
            end
        end
    end
        
    %% Tracking error
    E(:,k) = Y(:,k) - Y_d(:,k);
    
    %% Sliding surface
    S(:,k) = c_1*E(:,k) + c_2*E_I(:,k-1);
    
    E_I(:,k) = E(:,k-1) + sig_func(E(:,k),pp/qq);
    
    %% Integral terminal sliding mode controller
    delta_U_eq(:,k) = ((SAI_hat{k}+m)^-1)*((c_1^-1)*S(:,k)-((c_1^-1)*c_2)*E_I(:,k)+Y_d(:,k+1)-Y(:,k)-delta_d_hat(:,k));
%     delta_U_sw(:,k) = ((SAI_hat{k}+m)^-1)*(-c_s*(sign(S(:,k))));
    delta_U_sw(:,k) = ((SAI_hat{k}+m)^-1)*(-c_s*(sat_func(S(:,k),omega)));
    
    delta_U_s(:,k) = delta_U_eq(:,k) + delta_U_sw(:,k);

    %% Predictive controller
    S_hat_p = repmat(S(:,k)',1,N)';
    R_hat_p = sat_func(S_hat_p,repmat(omega,N,1));
    L_hat_p = repmat(delta_d_tilde(:,k-1)',1,N)';
    PHI_hat_p = tril(repmat((SAI_hat{k}+m),N,Nu));
    
%     delta_Up = -(((XI*PHI_hat_p)'*(XI*PHI_hat_p)+gamma*eye(Nu*nInputs))^-1)*(XI*PHI_hat_p)'*(kron(GAMMA,S(:,k))+PAI*R_hat_p+LAMBDA*L_hat_p);
    delta_Up = -((XI*PHI_hat_p)'*(kron(GAMMA,S(:,k))+PAI*R_hat_p+LAMBDA*L_hat_p))/(gamma+norm(XI*PHI_hat_p,2)^2);
    
    delta_U_mpc(:,k) = K*delta_Up;
    
    %% Final controller
    delta_U(:,k) = 1*delta_U_s(:,k) + 1*delta_U_mpc(:,k);
    
    U(:,k) = U(:,k-1) + delta_U(:,k);
    
    %% Plant
    dist(:,k) = [0.10*sin(5*pi*t(k))
                 0.14*sin(2*pi*t(k))];
    
%     dist(:,k) = [(0.12*sin(t(k)))*(t(k)>=3 && t(k)<=6)
%                  (0.12*sin(t(k)))*(t(k)>=3 && t(k)<=6)];

%     dist(:,k) = [(0.1)*(t(k)>=3 && t(k)<=6)
%                  (0.1)*(t(k)>=3 && t(k)<=6)];
             
    delta_d(:,k) = dist(:,k) - dist(:,k-1);

    qddot = PAM_twoLink_robot_dynamics(1*t(k),q,qdot,U(:,k),dist(:,k));
    qdot = qdot + dt*qddot;
    q = q + dt*qdot;
    
    Y(:,k+1) = qdot;
    
    %%
    delta_Y(:,k+1) = Y(:,k+1) - Y(:,k);
    
%     delta_d_tilde(:,k) = (delta_Y(:,k)-delta_Y(:,k+1)) - (SAI_hat{k-1}+m)*delta_U(:,k-1) + (SAI_hat{k}+m)*delta_U(:,k);
    delta_d_tilde(:,k) = delta_d_hat(:,k) - delta_d(:,k);
    
end

SAI_hat_mat = cell2mat(SAI_hat);
sai_hat_11 = SAI_hat_mat(1,1:2:end);
sai_hat_12 = SAI_hat_mat(1,2:2:end);
sai_hat_21 = SAI_hat_mat(2,1:2:end);
sai_hat_22 = SAI_hat_mat(2,2:2:end);

Y_d = Y_d(:,1:end-1);
Y = Y(:,1:end-1);

toc;
%% Plot results
LineWidth = 1.5;

%% Dynamics of the PPD
figure(1);
plot(t,sai_hat_11,'r','Linewidth',LineWidth);
hold on;
plot(t,sai_hat_12,'g','Linewidth',LineWidth);
plot(t,sai_hat_21,'b','Linewidth',LineWidth);
plot(t,sai_hat_22,'c','Linewidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('PPD estimation');
lg_handle = legend('$\hat{\psi}_{11}$','$\hat{\psi}_{12}$',...
                   '$\hat{\psi}_{21}$','$\hat{\psi}_{22}$');
lg_handle.Interpreter = 'LATEX';
lg_handle.FontName = 'Times';
lg_handle.FontSize = 14;

ttl_handle = title('$\hat{\Psi}(k)$ for CFDL');
ttl_handle.Interpreter = 'LATEX';
ttl_handle.FontName = 'Times';
ttl_handle.FontSize = 14;

%% Disturbance estimation performance
figure(2);
% subplot(2,1,1);
% plot(t,delta_d(1,:),'r','LineWidth',LineWidth);
% hold on;
% plot(t,delta_d_hat(1,:),'b--','LineWidth',LineWidth);
% grid on;
% xlabel('Time (sec)');
% ylabel('Disturbance increment');
% legend('Real incremental disturbance','Estimated incremental disturbance');

subplot(2,1,1);
plot(t,delta_d_tilde(1,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('$\Delta{\tilde{d}}(k)$','Interpreter','LATEX');
legend('Incremental disturbance estimation error');

title('Disturbance estimation performance');

subplot(2,1,2);
plot(t,delta_d_tilde(2,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('$\Delta{\tilde{d}}(k)$','Interpreter','LATEX');
legend('Incremental disturbance estimation error');

%% Sliding surfaces
figure(3);
subplot(2,1,1);
plot(t,S(1,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('S_1(t)');

title('Sliding surfaces');

subplot(2,1,2);
plot(t,S(2,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('S_2(t)');

%% Tracking errors
figure(4);
subplot(2,1,1);
plot(t,E(1,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('E_1(t)');

title('Tracking errors');

subplot(2,1,2);
plot(t,E(2,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('E_2(t)');

%% Control inputs
figure(5);
subplot(2,1,1);
plot(t,U(1,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('U_1(t)');

title('Control inputs');

subplot(2,1,2);
plot(t,U(2,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('U_2(t)');

%% Plant outputs
figure(6);
subplot(2,1,1);
plot(t,Y_d(1,:),'r--','LineWidth',2);
hold on;
plot(t,Y(1,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('Y_1(t)');
legend('Set-point','Measured');

title('Plant outputs');

subplot(2,1,2);
plot(t,Y_d(2,:),'r--','LineWidth',2);
hold on;
plot(t,Y(2,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('Y_2(t)');
legend('Set-point','Measured');
