clc
clear all
close all
clear memory

V= 0.1;
a=0.1;
e=0.1;
Iz=1;

Fz=9000;
c=-100000;
cFa=20;
cMa=-2;

k=0;
kappa=-270;
sigma=0.3;

alpha_g=10;
delta=5;

%Initialization
psi=0.1;
dpsi=0.1;
y1=0.1;

alpha = y1/sigma;

Mz = Fz*cMa*alpha_g/180*sin(180*alpha/alpha_g);

Fy = cFa*delta*Fz*(alpha/abs(alpha));

M1 = c * psi;

M2 = k * dpsi;

M3 = Mz - e * Fy;

M4 = kappa/V * dpsi;


%c1 = c/Iz;
%c2 = k/Iz + kappa/(V*Iz);
%c3 = Fz*(cMa-e*cFa)/Iz/sigma;
%c4 = V*psi;
%c5 = e - a;
%c6 = -V/sigma;

t=60;
dt = 1;
vec_time = (0:dt:t-1);
psi_vals = zeros(1,t/dt);
dpsi_vals = zeros(1,t/dt);
y1_vals = zeros(1,t/dt);
psi_vals(1) = psi;
dpsi_vals(1) = psi;
y1_vals(1) = psi;
V_vals=0:80/(t-1):80;

for i=2:t

    alpha = y1/sigma;

    if(alpha<=delta)
        Fy=cFa*alpha*Fz;
    else
        Fy=cFa*delta*Fz*(alpha/abs(alpha));
    end
    if(abs(alpha)<=alpha_g)
        Mz = Fz*cMa*alpha_g/180*sin(180*alpha/alpha_g);
    else
        Mz = 0;
    end

    c1 = c/Iz;
    c2 = k/Iz + kappa/(V_vals(i)*Iz);
    c3 = Fz*(cMa-e*cFa)/Iz/sigma;
    c4 = V_vals(i)*psi;
    c5 = e - a;
    c6 = -V_vals(i)/sigma;
    
    M1 = c * psi;
    
    M2 = k * dpsi;
    
    M3 = Mz - e * Fy;
    
    M4 = kappa/V_vals(i) * dpsi;
    
    y1 = c4 + c5*dpsi + c6*y1;
    y1_vals(i) = y1;
    
    psi = dpsi;
    psi_vals(i) = psi;

    dpsi = c1*psi + c2*dpsi + c3*y1_vals(i-1);
    dpsi_vals(i) = dpsi;
end

plot(vec_time,psi_vals,'--r','LineWidth',2)