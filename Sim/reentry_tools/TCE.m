function [partP,X,C,RhoMix,S,H,MwMix] = TCE(P,T,NS,NH,N2S,N2H,NOS,NOH,OS,OH,O2S,O2H,trange)

iter = 50;                                      % Number of iterations
P = P/101325;                                   % Convert to ATM
R = 8.3144598;                                  % KJ/kmol/k
calToJoule = 4.184;                             % Converts cal/mol to J/mol
AtmToPa  = 101325;                              % Converts Atm to Pascal
MwN = 14.0067;                                  % Kg/Kmol
MwN2 = MwN*2;
MwO = 15.999;
MwO2 = MwO*2;
MwNO = MwN+MwO;

%% Kp Derivation

indexT = length(T);                             % Counting Index for loops

Kp2N = zeros(1,indexT);
Kp2O = zeros(1,indexT);
KpNO = zeros(1,indexT);
H_N = zeros(1,indexT); S_N = zeros(1,indexT);
H_O = zeros(1,indexT); S_O = zeros(1,indexT);
H_NO = zeros(1,indexT); S_NO = zeros(1,indexT);
H_N2 = zeros(1,indexT); S_N2 = zeros(1,indexT);
H_O2 = zeros(1,indexT); S_O2 = zeros(1,indexT);

for i = 1:indexT

    Ti = T(i);
    indexRow = find(trange == Ti);   % Row Index

    G_N = NH(indexRow)-Ti*NS(indexRow);        % Gibbs Free Energy of N
    G_N2 = N2H(indexRow)-Ti*N2S(indexRow);     % Gibbs Free Energy of N2
    G_O = OH(indexRow)-Ti*OS(indexRow);
    G_O2 = O2H(indexRow)-Ti*O2S(indexRow);
    G_NO = NOH(indexRow)-Ti*NOS(indexRow);

    H_N(i) = NH(indexRow); S_N(i) = NS(indexRow);
    H_O(i) = OH(indexRow); S_O(i) = OS(indexRow);
    H_NO(i) = NOH(indexRow); S_NO(i) = NOS(indexRow);
    H_N2(i) = N2H(indexRow); S_N2(i) = N2S(indexRow);
    H_O2(i) = O2H(indexRow); S_O2(i) = O2S(indexRow);

    deltaG_Kp2N = calToJoule*(2*G_N-G_N2);
    deltaG_KpN2 = calToJoule*(G_N2-2*G_N);
    deltaG_Kp2O = calToJoule*(2*G_O-G_O2);
    deltaG_KpO2 = calToJoule*(G_O2-2*G_O);
    deltaG_KpNO = calToJoule*(G_NO-G_N-G_O);

    Kp2N(i) = exp(-deltaG_Kp2N/(R*Ti));
    KpN2(i) = exp(-deltaG_KpN2/(R*Ti));
    Kp2O(i) = exp(-deltaG_Kp2O/(R*Ti));
    KpO2(i) = exp(-deltaG_KpO2/(R*Ti));
    KpNO(i) = exp(-deltaG_KpNO/(R*Ti));

end

%% Thermodynamic Properties Calculation

A = [2/3 1/6;  1/3 -1/6];       % Transformation matrix for high temp iterations
A2 = [4/5 1/10; 1/5 -1/10];     % Transformation matrix for low temp iterations
start = 1; er = 0;

PN = zeros(1,iter);P_N = zeros(1,indexT);X_N = zeros(1,indexT);C_N = zeros(1,indexT);
PO = zeros(1,iter);P_O = zeros(1,indexT);X_O = zeros(1,indexT);C_O = zeros(1,indexT);
PNO = zeros(1,iter);P_NO = zeros(1,indexT);X_NO = zeros(1,indexT);C_NO = zeros(1,indexT);
PN2 = zeros(1,iter);P_N2 = zeros(1,indexT);X_N2 = zeros(1,indexT);C_N2 = zeros(1,indexT);
PO2 = zeros(1,iter);P_O2 = zeros(1,indexT);X_O2 = zeros(1,indexT);C_O2 = zeros(1,indexT);
MwMix = zeros(1,indexT);RMix = zeros(1,indexT);RhoMix = zeros(1,indexT);
H = zeros(1,indexT);S = zeros(1,indexT);

while start == 1 || er == 1
    start = 0;
    for ii = 1:indexT
        Ti = T(ii);
        Pi = P(ii);
        for jj = 1:iter
        j_1 = jj+1;
        if er == 0              % Lower temperature partial pressure iterations
            PN(j_1) = sqrt(PN2(jj)/KpN2(ii));
            PO(j_1) = sqrt(PO2(jj)/KpO2(ii));
            PNO(j_1) = KpNO(ii)*PN(jj)*PO(jj);
            B = A2*[Pi-PNO(j_1)-PO(j_1)-PN(j_1);...
            3*PNO(j_1)+4*PO(j_1)-PN(j_1)];
            PN2(j_1) = B(1);
            PO2(j_1)= B(2);
        else                    % High temperature partial pressure iterations
            PN(j_1) = sqrt(PN2(jj)*Kp2N(ii));
            PO2(j_1) = PO(jj)^2/Kp2O(ii);
            PNO(j_1) = KpNO(ii)*PN(jj)*PO(jj);
            B = A*[Pi-PNO(j_1)-PO2(j_1)-PN(j_1);...
            3*PNO(j_1)+8*PO2(j_1)-PN(j_1)];
            PN2(j_1) = B(1);
            PO(j_1)= B(2);
        end
    end
    end

    for ii = 1:indexT
    P_N(ii) = PN(iter+1);
    P_O(ii) = PO(iter+1);
    P_NO(ii) = PNO(iter+1);
    P_N2(ii) = PN2(iter+1);
    P_O2(ii) = PO2(iter+1);

    X_N(ii) = P_N(ii)/Pi;
    X_O(ii) = P_O(ii)/Pi;
    X_NO(ii) = P_NO(ii)/Pi;
    X_N2(ii) = P_N2(ii)/Pi;
    X_O2(ii) = P_O2(ii)/Pi;

    MwMix(ii) = X_N(ii)*MwN+X_O(ii)*MwO+X_NO(ii)*MwNO+X_N2(ii)*MwN2+...
        X_O2(ii)*MwO2;
    RMix(ii) = R/MwMix(ii)*1e3;
    RhoMix(ii) = Pi*AtmToPa/(RMix(ii)*Ti);

    H(ii) = (X_N(ii)*H_N(ii)+X_O(ii)*H_O(ii)+X_NO(ii)*H_NO(ii)+X_N2(ii)*...
        H_N2(ii)+X_O2(ii)*H_O2(ii))*calToJoule/MwMix(ii)*1e3;
    S(ii) = (X_N(ii)*S_N(ii)+X_O(ii)*S_O(ii)+X_NO(ii)*S_NO(ii)+X_N2(ii)*...
        S_N2(ii)+X_O2(ii)*S_O2(ii))*calToJoule/MwMix(ii)*1e3;

    C_N(ii) = X_N(ii)*MwN/MwMix(ii);
    C_O(ii) = X_O(ii)*MwO/MwMix(ii);
    C_NO(ii) = X_NO(ii)*MwNO/MwMix(ii);
    C_N2(ii) = X_N2(ii)*MwN2/MwMix(ii);
    C_O2(ii) = X_O2(ii)*MwO2/MwMix(ii);

    partP = [P_N(ii) P_NO(ii) P_N2(ii) P_O(ii) P_O2(ii)];
    X = [X_N(ii) X_NO(ii) X_N2(ii) X_O(ii) X_O2(ii)];
    C = [C_N(ii) C_NO(ii) C_N2(ii) C_O(ii) C_O2(ii)];
    end

    notNum = sum(isnan(partP));
    neg = any(partP<0);
    imaginary = isreal(partP);

    if imaginary == 0
        check = sum(imag(partP));
        if check<1e-10
            imaginary = 1;
            partP = real(partP);
        end
    end
    if notNum == 1 || neg == 1 || imaginary == 0
        er = 1;
    else
        er = 0;
    end
end

end