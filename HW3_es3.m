close all
clear
clc

%% Load workspace data
load('ws_homework_3_2025.mat');


%% System parameters
m = 1.5; % Nominal mass [kg]
J = diag([1.2416, 1.2416, 2*1.2416]); % Inertia matrix [kg·m²]
dt = 0.001; % Sampling time [s]
g = 9.81; % Gravitational acceleration [m/s²]

%% Extract measured data
thrust = thrust.signals.values; % Commanded thrust [N]
tau = tau.signals.values; % Commanded torques [Nm]
linear_vel = linear_vel.signals.values; % Linear velocity [m/s]
attitude = attitude.signals.values; % Euler angles [rad]
attitude_vel = attitude_vel.signals.values; % Euler angle rates [rad/s]
time = (0:length(thrust)-1)*dt; % Time vector
t_max=length(time);
% 
% R=[2,5,20,150];
% for RRR=1:length(R)
% 
%     %% ESTIMATOR
%     c0 = 10; 
%     k0 = c0;
%     r = R(RRR);  
%     % Transfer function
%     s = tf('s');
%     G = (k0/(s+c0))^r;
%     c = cell2mat(G.Denominator); % [c0 c_1 c_2 ... c_r-1 c_r]
% 
%     % K computation
%     K = zeros(r,1);
%     app = 1; 
%     for i = 1:r
%         K(i) = c(i+1)/app;
%         app = app*K(i);
%     end
%     K=flip(K);
% 
% 
%     %% Estimator variables 
%     external_wrench = zeros(6, length(time));
%     gamma = zeros(6, length(time), r); 
% 
%     for k = 1:length(time)-1
% 
% 
%         [C, Q, R_b, M] = compute_all_matrices(attitude(k,:), attitude_vel(k,:), J);
%         
%         q(:,k+1) = [m*eye(3) zeros(3,3); zeros(3,3) M]*[linear_vel(k+1,:) attitude_vel(k+1,:)]';
%         qdot_meno_ft=[m*g*[0 0 1]'-thrust(k)*R_b*[0 0 1]'; C'*attitude_vel(k,:)' + Q'*tau(k,:)'];
%         for i = 1:r
%             if i == 1
%                 gamma(:,k+1,1) = gamma(:,k,1) + K(1)*( (q(:,k+1) - q(:,k)) - (external_wrench(:,k) + qdot_meno_ft)*dt);
%             else
%                 gamma(:,k+1,i) = gamma(:,k,i) + K(i)*dt*(-external_wrench(:,k) + gamma(:,k,i-1));
%             end
%         end
%         external_wrench(:,k+1) = gamma(:,k+1,r);
%     end
%     external_wrench_R(RRR,:,:)=external_wrench(:,:);
% end
% %% Plots 
% 
% subplot(3,2,1)
% hold on
% plot(time, squeeze(external_wrench_R(1,1,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(2,1,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(3,1,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(4,1,:)), 'LineWidth', 1.5)
% plot(time, 1*ones(size(time)), '--', 'LineWidth', 1.5)
% hold off
% xlabel('t', 'Interpreter', 'latex', 'FontSize', 12)
% ylabel('fx', 'Interpreter', 'latex', 'FontSize', 12)
% legend('  c0=5','  c0=10',' c0=50'','  c0=100'','Expected dist.', 'Interpreter', 'latex')
% grid on
% title('X-axis Force Disturbance')
% 
% subplot(3,2,2)
% hold on
% plot(time, squeeze(external_wrench_R(1,4,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(2,4,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(3,4,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(4,4,:)), 'LineWidth', 1.5)
% hold off
% xlabel('t', 'Interpreter', 'latex', 'FontSize', 12)
% ylabel('$\tau_x$', 'Interpreter', 'latex', 'FontSize', 12)
% legend('  c0=5','  c0=10',' c0=50'','  c0=100'', 'Interpreter', 'latex')
% grid on
% title('x Torque Disturbance')
% 
% subplot(3,2,3)
% hold on
% plot(time, squeeze(external_wrench_R(1,2,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(2,2,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(3,2,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(4,2,:)), 'LineWidth', 1.5)
% plot(time, 1*ones(size(time)), '--', 'LineWidth', 1.5)
% hold off
% xlabel('t', 'Interpreter', 'latex', 'FontSize', 12)
% ylabel('fy', 'Interpreter', 'latex', 'FontSize', 12)
% legend('  c0=5','  c0=10',' c0=50'','  c0=100'','Expected dist.', 'Interpreter', 'latex')
% 
% grid on
% title('Y-axis Force Disturbance')
% 
% subplot(3,2,4)
% hold on
% plot(time, squeeze(external_wrench_R(1,5,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(2,5,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(3,5,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(4,5,:)), 'LineWidth', 1.5)
% hold off
% xlabel('t', 'Interpreter', 'latex', 'FontSize', 12)
% ylabel('$\tau_y$', 'Interpreter', 'latex', 'FontSize', 12)
% legend('  c0=5','  c0=10',' c0=50'','  c0=100'', 'Interpreter', 'latex')
% grid on
% title('y Torque Disturbance')
% 
% subplot(3,2,5)
% hold on
% plot(time, squeeze(external_wrench_R(1,3,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(2,3,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(3,3,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(4,3,:)), 'LineWidth', 1.5)
% hold off
% xlabel('t', 'Interpreter', 'latex', 'FontSize', 12)
% ylabel('fz', 'Interpreter', 'latex', 'FontSize', 12)
% legend('  c0=5','  c0=10',' c0=50'','  c0=100'', 'Interpreter', 'latex')
% grid on
% title('Z-axis Force Disturbance')
% 
% subplot(3,2,6)
% hold on
% plot(time, squeeze(external_wrench_R(1,6,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(2,6,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(3,6,:)), 'LineWidth', 1.5)
% plot(time, squeeze(external_wrench_R(4,6,:)), 'LineWidth', 1.5)
% plot(time, -0.4*ones(size(time)), '--', 'LineWidth', 1.5)
% hold off
% xlabel('t', 'Interpreter', 'latex', 'FontSize', 12)
% ylabel('$\tau_z$', 'Interpreter', 'latex', 'FontSize', 12)
% legend('  c0=5','  c0=10',' c0=50'','  c0=100'','Expected dist.', 'Interpreter', 'latex')
% grid on
% title('z Torque Disturbance')
% 
% sgtitle('Disturbance Estimation', 'FontSize', 20)


R=[5,10,50,100];
for RRR=1:length(R)

    %% ESTIMATOR
    c0 = R(RRR); 
    k0 = c0;
    r = 20;  
    % Transfer function
    s = tf('s');
    G = (k0/(s+c0))^r;
    c = cell2mat(G.Denominator); 
    
    % K computation
    K = zeros(r,1);
    app = 1; 
    for i = 1:r
        K(i) = c(i+1)/app;
        app = app*K(i);
    end
    K=flip(K);
    
    
    %% Estimator variables 
    external_wrench = zeros(6, length(time));
    gamma = zeros(6, length(time), r); 
    
    for k = 1:length(time)-1
    
        [C, Q, R_b, M] = compute_all_matrices(attitude(k,:), attitude_vel(k,:), J);

        %momentum
        q(:,k+1) = [m*eye(3) zeros(3,3); zeros(3,3) M]*[linear_vel(k+1,:) attitude_vel(k+1,:)]';
        qdot_meno_ft=[m*g*[0 0 1]'-thrust(k)*R_b*[0 0 1]'; C'*attitude_vel(k,:)' + Q'*tau(k,:)'];
        for i = 1:r
            if i == 1
                gamma(:,k+1,1) = gamma(:,k,1) + K(1)*( (q(:,k+1) - q(:,k)) - (external_wrench(:,k) + qdot_meno_ft)*dt);
            else
                gamma(:,k+1,i) = gamma(:,k,i) + K(i)*dt*(-external_wrench(:,k) + gamma(:,k,i-1));
            end
        end
        external_wrench(:,k+1) = gamma(:,k+1,r);
    end
    external_wrench_R(RRR,:,:)=external_wrench(:,:);
end
%% Plots 

subplot(3,2,1)
hold on
plot(time, squeeze(external_wrench_R(1,1,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(2,1,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(3,1,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(4,1,:)), 'LineWidth', 1.5)
plot(time, 1*ones(size(time)), '--', 'LineWidth', 1.5)
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('fx', 'Interpreter', 'latex', 'FontSize', 12)
legend('  c0=5','  c0=10',' c0=50','  c0=100','Expected dist.', 'Interpreter', 'latex')
grid on
title('X-axis Force Disturbance')

subplot(3,2,2)
hold on
plot(time, squeeze(external_wrench_R(1,4,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(2,4,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(3,4,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(4,4,:)), 'LineWidth', 1.5)
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\tau_x$', 'Interpreter', 'latex', 'FontSize', 12)
legend('  c0=5','  c0=10',' c0=50','  c0=100', 'Interpreter', 'latex')
grid on
title('x Torque Disturbance')

subplot(3,2,3)
hold on
plot(time, squeeze(external_wrench_R(1,2,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(2,2,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(3,2,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(4,2,:)), 'LineWidth', 1.5)
plot(time, 1*ones(size(time)), '--', 'LineWidth', 1.5)
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('fy', 'Interpreter', 'latex', 'FontSize', 12)
legend('  c0=5','  c0=10',' c0=50','  c0=100','Expected dist.', 'Interpreter', 'latex')

grid on
title('Y-axis Force Disturbance')

subplot(3,2,4)
hold on
plot(time, squeeze(external_wrench_R(1,5,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(2,5,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(3,5,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(4,5,:)), 'LineWidth', 1.5)
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\tau_y$', 'Interpreter', 'latex', 'FontSize', 12)
legend('  c0=5','  c0=10',' c0=50','  c0=100', 'Interpreter', 'latex')
grid on
title('y Torque Disturbance')

subplot(3,2,5)
hold on
plot(time, squeeze(external_wrench_R(1,3,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(2,3,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(3,3,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(4,3,:)), 'LineWidth', 1.5)
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('fz', 'Interpreter', 'latex', 'FontSize', 12)
legend('  c0=5','  c0=10',' c0=50','  c0=100', 'Interpreter', 'latex')
grid on
title('Z-axis Force Disturbance')

subplot(3,2,6)
hold on
plot(time, squeeze(external_wrench_R(1,6,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(2,6,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(3,6,:)), 'LineWidth', 1.5)
plot(time, squeeze(external_wrench_R(4,6,:)), 'LineWidth', 1.5)
plot(time, -0.4*ones(size(time)), '--', 'LineWidth', 1.5)
hold off
xlabel('t', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\tau_z$', 'Interpreter', 'latex', 'FontSize', 12)
legend('  c0=5','  c0=10',' c0=50','  c0=100','Expected dist.', 'Interpreter', 'latex')
grid on
title('z Torque Disturbance')

sgtitle('Disturbance Estimation', 'FontSize', 20)

% saveas(gcf, 'Hw3_es3_c0.png')
% saveas(gcf, 'Hw3_es3_c0.pdf')

% %% Real mass of the UAV
% real_mass = -external_wrench(3,end)/g + m % computation using the z-disturbance 
% 
% %%??
% u_tb = thrust(end)*compute_R_b(attitude(end,:))*[0 0 1]';
% mass_r_z_1=u_tb(3)/g; % not using disturbance but the steady-state equation
% 
% 
% saveas(gcf, 'Hw3_es3_r.png')
% saveas(gcf, 'Hw3_es3_r.pdf')
% 
% 
% 

%% Real mass of the UAV
real_mass = external_wrench_R(3,3,end)/g + m