% Autotune a MIMO PID for the nominal case by LMIs
% [S. Boyd, M. Hast, and J. Astrom, "MIMO PID tuning via iterated LMI
% restiction", Inter. Jur. of Robust & Nonlinear Control, 2016.]

clc
clear variables
addlib yalmip sdpt3

s = tf('s');
P = [ 12.8*exp(-s)/(16.7*s+1)   -18.9*exp(-3*s)/(21*s+1)
      6.6*exp(-7*s)/(10.9*s+1)  -19.4*exp(-3*s)/(14.2*s+1) ];
N=300;
w = logspace(-3,3,N);
P0=dcgain(P);

tau=0.3;

% -------------------------------------
% specs: 
% sensitivity (4db))
Smax = 1.4*ones(1,N);               
% low frequency sensitivity (4db)
Tmax = 1.4*ones(1,N);               
% cost of feedback 
Qmax = 3/min(svd(dcgain(P)))*ones(1,N);   

% Option.Structure.Kp=[1 0 ; 0 1];
% Option.Structure.Ki=[1 0 ; 0 1];
% Option.Structure.Kd=[1 0 ; 0 1];
% Option.Initial.Kp=[0.001 0 ; 0 -0.001];
% Option.Initial.Ki=[0.001 0 ; 0 -0.001];
% Option.Initial.Kd=[0 0 ; 0 0];
Option=[];

G = auto_mimo_pid( P,w,Smax,Tmax,Qmax,tau,Option )

%% analysis
L=frd(P*G,w);
S=(eye(2)+P*G)^-1;
T=P*G*S;
Q=G*S;
figure(1), sigma(S,w); hold on; semilogx(w,20*log10(Smax),'r');
figure(2), sigma(T,w); hold on; semilogx(w,20*log10(Tmax),'r');
figure(3), sigma(Q,w); hold on; semilogx(w,20*log10(Qmax),'r');
