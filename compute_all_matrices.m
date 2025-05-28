function [C, Q, R_b, M] = compute_all_matrices(eulero, eulero_vel, J)
    
    %% 1. Matrice di rotazione R_b 
    Rz = [cos(eulero(3)) -sin(eulero(3)) 0;  % Rotazione Z (yaw)
         sin(eulero(3))  cos(eulero(3)) 0;
         0       0      1];
     
    Ry = [cos(eulero(2))  0 sin(eulero(2));  % Rotazione Y (pitch)
         0       1 0;
         -sin(eulero(2)) 0 cos(eulero(2))];
     
    Rx = [1 0      0;        % Rotazione X (roll)
         0 cos(eulero(1)) -sin(eulero(1));
         0 sin(eulero(1))  cos(eulero(1))];
    
    R_b = Rz * Ry * Rx;  % Ordine moltiplicativo ZYX
    
    %% 2. Matrice Q 
    Q = [1,  0,        -sin(eulero(2));
         0,   cos(eulero(1)),    sin(eulero(1))*cos(eulero(2));
         0, -sin(eulero(1)),     cos(eulero(1))*cos(eulero(2))];

    Q_dot=[0 0                           -cos(eulero(2))*eulero_vel(2)
           0 -sin(eulero(1))*eulero_vel(1) -sin(eulero(2))*sin(eulero(1))*eulero_vel(2)+cos(eulero(2))*cos(eulero(1))*eulero_vel(1);
           0 -cos(eulero(1))*eulero_vel(1) -sin(eulero(2))*cos(eulero(1))*eulero_vel(2)-cos(eulero(2))*sin(eulero(1))*eulero_vel(1)];

    
    %% 3. Matrice C 
    w=Q*eulero_vel';
    if(numel(w)~= 1)
        S = [0    -w(3)  w(2);
             w(3)  0    -w(1);
            -w(2)  w(1)  0];
    else
        S= zeros(Q*eulero_vel');
    end
    
    C = Q'*S*J*Q+Q'*J*Q_dot;
    
    %% 4. Matrice M
     M = Q'*J*Q;
end