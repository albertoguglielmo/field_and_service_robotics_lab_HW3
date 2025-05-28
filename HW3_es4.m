close all
clear all


out = sim('geometric_control_template');


% 1. Creazione del vettore tempo
time = out.tout;

% 2. Estrazione dei segnali (assumendo 13 segnali in uscita)
num_signals = 13+3;
for i = 1:num_signals
    data(:,i) = out.simout(:,i);
end

% 3. Assegnazione dei vettori specifici
% Errori lineari (posizione)
err_X = data(:,1);  % Errore asse X
err_Y = data(:,2);  % Errore asse Y
err_Z = data(:,3);  % Errore asse Z

% Errori rotazionali in SO(3) 
rot_err_X = data(:,4);  % Errore rotazionale X
rot_err_Y = data(:,5);  % Errore rotazionale Y
rot_err_Z = data(:,6);  % Errore rotazionale Z

% Coppie di controllo
tau_b_X = data(:,7);  % Coppia X body frame
tau_b_Y = data(:,8);  % Coppia Y body frame
tau_b_Z = data(:,9);  % Coppia Z body frame

% Spinta totale
uT = data(:,10);     % Spinta totale

% Accelerazioni lineari
acc_X = data(:,11);  % Accelerazione X
acc_Y = data(:,12);  % Accelerazione Y
acc_Z = data(:,13);  % Accelerazione Z


acc_X = data(:,11);  % Accelerazione X
acc_Y = data(:,12);  % Accelerazione Y
acc_Z = data(:,13);  % Accelerazione Z

yaw = data(:,14);  % Accelerazione X
pitch = data(:,15);  % Accelerazione Y
roll = data(:,16);  % Accelerazione Z


figure;

%% 1. Position Error (X, Y, Z)
subplot(2,2,1);
plot(time, err_X, 'r', 'LineWidth', 2.0); hold on;
plot(time, err_Y, 'g', 'LineWidth', 2.0);
plot(time, err_Z, 'b', 'LineWidth', 2.0);
xlabel('Time [s]');
ylabel('Error [m]');
title('Position Error');
legend('X', 'Y', 'Z');
grid on;

%% 2. Rotational Error (X, Y, Z)
subplot(2,2,2);
plot(time, rot_err_X, 'r', 'LineWidth', 2.0); hold on;
plot(time, rot_err_Y, 'g', 'LineWidth', 2.0);
plot(time, rot_err_Z, 'b', 'LineWidth', 2.0);
xlabel('Time [s]');
ylabel('Error [rad]');
title('Rotational Error');
legend('Rot X', 'Rot Y', 'Rot Z');
grid on;

%% 3. Control Torques
subplot(2,2,3);
plot(time, tau_b_X, 'r', 'LineWidth', 2.0); hold on;
plot(time, tau_b_Y, 'g', 'LineWidth', 2.0);
plot(time, tau_b_Z, 'b', 'LineWidth', 2.0);
xlabel('Time [s]');
ylabel('Torque [Nm]');
title('Control Torques (\tau)');
legend('\tau_X', '\tau_Y', '\tau_Z');
grid on;

%% 4. Total Thrust
subplot(2,2,4);
plot(time, uT, 'b', 'LineWidth', 2.0);
xlabel('Time [s]');
ylabel('Thrust [N]');
title('Total Thrust (u_T)');
legend('u_T');
grid on;


% saveas(gcf, 'Hw3_es4.png')
% saveas(gcf, 'Hw3_es4.pdf')

plot(time, roll*180/pi, 'r', 'LineWidth', 2.0); hold on;
plot(time, pitch*180/pi, 'g', 'LineWidth', 2.0);
plot(time, yaw*180/pi, 'b', 'LineWidth', 2.0);
xlabel('Time [s]');
ylabel('degrees  [deg]');
title('Orientation');
legend('roll', 'pitch', 'yaw');
grid on;

% saveas(gcf, 'Hw3_es4_or.png')
% saveas(gcf, 'Hw3_es4_or.pdf')