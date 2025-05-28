        % Definisci le variabili simboliche con nomi w1,w2,...,w8 e k1,k2,...,k8
    syms w1 w2 w3 w4 w5 w6 w7 w8
    syms k1 k2 k3 k4 k5 k6 k7 k8
    syms Cm Cf
    
    % Parametri fissi
    zPi = [0; 0; -1];  % Vettore zPi
    n = 8;            % Numero di elementi
    
    % Inizializzazione del vettore tau_b simbolico
    tau_b = sym(zeros(3, 1));
    
    % Crea array delle variabili w e k per comodità di accesso
    w = [w1, w2, w3, w4, w5, w6, w7, w8];
    k = [k1, k2, k3, k4, k5, k6, k7, k8];
    
    for i = 1:n
        % Calcolo dell'angolo (0, 45, 90, ..., 315 gradi)
        angle_rad = (i-1)*pi/4;
        
        % Costruzione del vettore pPi [cos(angolo), sin(angolo), 0]
        pPi = [cos(angle_rad); sin(angle_rad); 0];
        
        % Matrice skew-simmetrica di pPi
        S_pPi = [0       -pPi(3)  pPi(2);
                 pPi(3)   0       -pPi(1);
                -pPi(2)   pPi(1)  0];
        
        % Calcolo del termine per l'i-esimo elemento
        term = abs(w(i)) * w(i) * (-k(i) * Cm * zPi + Cf * S_pPi * zPi);
        if max(abs(coeffs(term(1)))) < 1e-10
            term(1) = 0;
        end
        if max(abs(coeffs(term(2)))) < 1e-10
            term(2) = 0;
        end
        if max(abs(coeffs(term(3)))) < 1e-10
            term(3) = 0;
        end
        % Accumulo nel risultato simbolico
        tau_b = tau_b + term;
    end
    
    % Semplifica e approssima i coefficienti con vpa
    %tau_b = vpa(simplify(tau_b), 6);  % 6 cifre significative
    
    % Semplifica l'espressione finale
    tau_b = simplify(tau_b);

    
    % Visualizza il risultato
    disp('Tau_b simbolico:');
    pretty(tau_b);