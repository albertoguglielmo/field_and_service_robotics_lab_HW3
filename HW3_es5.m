close all
clear all


out = sim('tilting_control_voliro_template_friends');


% 1. Creazione del vettore tempo
time = out.tout;

% 2. Estrazione dei segnali (assumendo 13 segnali in uscita)
num_signals = 20;
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


%POS 
pos_x = data(:,7);  
pos_y = data(:,8);  
pos_z = data(:,9); 

%ORIENTAMENTO 
yaw = abs(data(:,10));  
pitch = data(:,11);  
roll = data(:,12);  

% uw
uw1 = data(:,13); 
uw2 = data(:,14);  
uw3 = data(:,15);  
uw4 = data(:,16);

%alphaw 
alphaw1=data(:,17);
alphaw2=data(:,18);
alphaw3=data(:,19);
alphaw4=data(:,20);


figure;

% 1. Position Errors (X, Y, Z)
subplot(3,2,1);
plot(time, err_X, 'r', 'LineWidth', 1.8); hold on;
plot(time, err_Y, 'g', 'LineWidth', 1.8);
plot(time, err_Z, 'b', 'LineWidth', 1.8);
xlabel('Time [s]');
ylabel('Error [m]');
title('Position Errors');
legend('X', 'Y', 'Z');
grid on;

% 2. Rotational Errors (X, Y, Z)
subplot(3,2,2);
plot(time, rot_err_X, 'r', 'LineWidth', 1.8); hold on;
plot(time, rot_err_Y, 'g', 'LineWidth', 1.8);
plot(time, rot_err_Z, 'b', 'LineWidth', 1.8);
xlabel('Time [s]');
ylabel('Error [rad]');
title('Rotational Errors');
legend('Roll', 'Pitch', 'Yaw');
grid on;

% 3. UAV Position (X, Y, Z)
subplot(3,2,3);
plot(time, pos_x, 'r', 'LineWidth', 1.8); hold on;
plot(time, pos_y, 'g', 'LineWidth', 1.8);
plot(time, pos_z, 'b', 'LineWidth', 1.8);
xlabel('Time [s]');
ylabel('Position [m]');
title('UAV Position');
legend('X', 'Y', 'Z');
grid on;

% 4. Orientation Angles (Roll, Pitch, Yaw)
subplot(3,2,4);
plot(time, roll, 'r', 'LineWidth', 1.8); hold on;
plot(time, pitch, 'g', 'LineWidth', 1.8);
plot(time, yaw, 'b', 'LineWidth', 1.8);
xlabel('Time [s]');
ylabel('Angle [rad]');
title('Orientation (RPY)');
legend('Roll', 'Pitch', 'Yaw');
grid on;

% 5. Propeller Angular Speeds
subplot(3,2,5);
plot(time, uw1, 'r', 'LineWidth', 1.8); hold on;
plot(time, uw2, 'g', 'LineWidth', 1.8);
plot(time, uw3, 'b', 'LineWidth', 1.8);
plot(time, uw4, 'k', 'LineWidth', 1.8);
xlabel('Time [s]');
ylabel('Speed [rad/s]');
title('Propeller Speeds (uw)');
legend('uw1', 'uw2', 'uw3', 'uw4');
grid on;

% 6. Tilting Angles (Alpha)
subplot(3,2,6);
plot(time, alphaw1, 'r', 'LineWidth', 1.8); hold on;
plot(time, alphaw2, 'g', 'LineWidth', 1.8);
plot(time, alphaw3, 'b', 'LineWidth', 1.8);
plot(time, alphaw4, 'k', 'LineWidth', 1.8);
xlabel('Time [s]');
ylabel('Angle [rad]');
title('Tilting Angles (\alpha)');
legend('\alpha1', '\alpha2', '\alpha3', '\alpha4');
grid on;

% saveas(gcf, 'Hw3_es5.png')
% saveas(gcf, 'Hw3_es5.pdf')