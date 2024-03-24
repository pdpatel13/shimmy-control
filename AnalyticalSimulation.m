clc;clear;close all;
syms t tau x1 x2 x3;
c = -100000;  % torsional spring constant Nm/rad
I = 1;        % z moment of inertia kg m^2
k = -50;      % torsional damping constant Nm/rad/s
kappa = -270; % tread width moment constant Nm^2 /rad
V = 80;       % velocity m/s
C_ma = -2;    % m/rad
e = 0.1;      % caster length
C_fa = 20;    % side force derivative 1/rad
F_z = 9000;   % vertical force N
a = 0.1;      % half the contact length m
sigma = 3*a;  % relaxation length m

A = [0,1,0;c/I, k/I + kappa/(V*I), (C_ma - e*C_fa)*F_z/(sigma*I); V, e-a, -V/sigma];

B = [0;1/I;0];
C = eye(3);
D = [0;0;0];
sys = ss(A,B,C,D);

x0 = [0, 20, 0];
initial(sys,x0)

[V, D] = eig(A);

phi = V*exp(D*t)*inv(V);%Not nice equation try to find by hand
vpa(phi*[x1;x2;x3],3)%Not nice equation try to find by hand

syms z a b c d e n
a = A(2,1);
b = A(2,2);
c = A(2,3);
d = A(3,1);
e = A(3,3);
F2 = (1-e*z^-1)/((1 - a*z^-2 - b*z^-1)*(1-e*z^-1) - c*d*z^-3);
F3 = (d*z^-2)/((1 - a*z^-2 - b*z^-1)*(1-e*z^-1) - c*d*z^-3);
f2 = iztrans(F2);
f3 = iztrans(F3);
f1 = subs(f2,n,n-1);
vpa(f2)
E1 = ztrans(f1/factorial(n));
E2 = ztrans(f2/factorial(n));
E3 = ztrans(f3/factorial(n));

E1 = vpa(subs(E1,z,t^-1),4);
E2 = vpa(subs(E2,z,t^-1),4);
E3 = vpa(subs(E3,z,t^-1),4);

% E1 = - 0.000242*exp(-318.7*t) + exp(t*(- 0.6758 - 337.3i))*(0.000121 + 0.001368i) + exp(t*(- 0.6758 + 337.3i))*(0.000121 - 0.001368i)
E1 = - 0.000242*exp(-318.7*t) + 2*exp(t*(- 0.6758))*(0.000121*cos(337.3*t) + 0.001368*sin(337.3*t));
E2 = - 0.07713*exp(-318.7*t)*exp(-318.7*t) + 2*exp(t*(- 0.6758))*(0.4614*cos(337.3*t) - 0.04175*sin(337.3*t));
E3 = - 0.0003722*exp(-318.7*t) + 2*exp(t*(- 0.6758))*(0.0001861*cos(337.3*t) - 0.0001754*sin(337.3*t));
