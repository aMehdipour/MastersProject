%% AE 280: Hypersonics
% Final Project
clc;clear;close all;

%% Inputs

% answer = inputdlg({...
%         'Please enter one altitude (meters):',...
%         'Please enter one freestream Mach Number:',...
%         'Please enter the Radius of Curvature (meters) of the blunt body:',...
%         'Please enter the maximum acceptable temperature  (K) of the material at the stagnation point:'...
%         },'Input Array',[1 1e2;1 1e2;1 1e2;1 1e2]);
% alt = str2double(answer{1});
% minf = str2double(answer{2});
% radius = str2double(answer{3});
% tw = str2double(answer{4});

alt = 80000;
minf = 9;
radius = 2;
tw = 3000;

[rhoinf,ainf,tinf,pinf,~,~] = atmos(alt);
vinf = minf*ainf;

%% Constants

trange = 100:1:6000;
N = csvread('N.csv');
N2 = csvread('N2.csv');
NO = csvread('NO.csv');
NO(21,1) = 1900;
O = csvread('O.csv');
O2 = csvread('O2.csv');

NS = interp1(N(2:end,1),N(2:end,4),trange,'spline')';
NH = interp1(N(2:end,1),N(2:end,6),trange,'spline')';
N2S = interp1(N2(2:end,1),N2(2:end,4),trange,'spline')';
N2H = interp1(N2(2:end,1),N2(2:end,6),trange,'spline')';
NOS = interp1(NO(2:end,1),NO(2:end,4),trange,'spline')';
NOH = interp1(NO(2:end,1),NO(2:end,6),trange,'spline')';
OS = interp1(O(2:end,1),O(2:end,4),trange,'spline')';
OH = interp1(O(2:end,1),O(2:end,6),trange,'spline')';
O2S = interp1(O2(2:end,1),O2(2:end,4),trange,'spline')';
O2H = interp1(O2(2:end,1),O2(2:end,6),trange,'spline')';

gainf = 1.4;
tStart = tinf*(1+2*gainf/(gainf+1)*(minf^2-1))*((2+(gainf-1)*minf^2)/...
    ((gainf+1)*minf^2));
if tStart>6000
    tStart = 6000;
else
    tStart = round(tStart);
end

[~,~,~,~,~,hinf,MwMinf] = TCE(round(pinf),round(tinf),NS,NH,N2S,...
            N2H,NOS,NOH,OS,OH,O2S,O2H,trange);

%% Calculation of properties downstream of normal shock

zeta = [0.1;0.101];
count1 = 0; count2 =0; count3 = 0;
RhoMix1 = 0; RhoMix2 = RhoMix1;

if tStart>4500
    tol = 1e-5;
else
    tol = 1e-4;
end

tic
while abs(zeta(2)-zeta(1))>tol
    p2z1 = pinf+rhoinf*vinf^2*(1-zeta(1));
    rho2z1 = rhoinf/zeta(1);
    h2z1 = hinf+vinf^2/2*(1-zeta(1)^2);
    p2z2 = pinf+rhoinf*vinf^2*(1-zeta(2));
    rho2z2 = rhoinf/zeta(2);
    h2z2 = hinf+vinf^2/2*(1-zeta(2)^2);
    count1 = 0; count2 =0;

    while abs(rho2z1-RhoMix1)>0.001
        tIter1 = tStart-count1;
        [partP1,X1,C1,RhoMix1,S1,H2_z1,MwMix1] = TCE(p2z1,tIter1,NS,NH,N2S,...
            N2H,NOS,NOH,OS,OH,O2S,O2H,trange);
        count1 = count1+1;
    end

    while abs(rho2z2-RhoMix2)>0.001
        tIter2 = tStart-count2;
        [partP2,X2,C2,RhoMix2,S2,H2_z2,MwMix2] = TCE(p2z2,tIter2,NS,NH,N2S,...
            N2H,NOS,NOH,OS,OH,O2S,O2H,trange);
        count2 = count2+1;
    end

    fz1 = h2z1-H2_z1;
    fz2 = h2z2-H2_z2;
    zeta3 = zeta(2)-fz2*(zeta(2)-zeta(1))/(fz2-fz1)
    zeta(1) = zeta(2);
    zeta(2) = zeta3;
end
toc

r2 = RhoMix2;
v2 = vinf*zeta(2);
t2 = tIter2;
p2 = p2z2;
h2 = H2_z2;
s2 = S2;
pComp2 = partP2;
xComp2 = X2;
cComp2 = C2;
po2 = p2+r2*v2^2/2;
ht = hinf+vinf^2/2;

Ttinf = tinf*(1+(gainf-1)/2*minf^2)
H2_z3 = 0;

while abs(H2_z3-ht)>1e3
        tIter3 = tStart-count3;
        [partP3,X3,C3,RhoMix3,S3,H2_z3,MwMix1] = TCE(po2,tIter3,NS,NH,N2S,...
            N2H,NOS,NOH,OS,OH,O2S,O2H,trange);
        count3 = count3+1;
end

tIter3



