clear,clc,close all
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
% Order of polymomials used for approximation 

NODETOL = 1e-10;
FinalTime = 10;


xmin = 0.0;
xmax = 2*pi;
A = 1;
T = 2*pi;
g = 9.81;
omega = 2*pi/FinalTime;
k = 2*pi/(abs(xmax-xmin));
h0 = omega^2/(k^2*g);
etasol = @(x,t,A,omega,k) A*cos(omega*t-k*x); 
usol = @(x,t,A,omega,k,h0) omega/(k*h0)*A*cos(omega*t-k*x);
D = [sqrt(g*h0),0;0,-sqrt(g*h0)];
param.g      = g;
param.A      = A;
param.h0     = h0;
param.omega  = omega;
param.k      = k;
param.etasol = etasol;
param.usol   = usol;
param.D      = D;
param.xmin   = xmin;
param.xmax   = xmax;

%% Periodic P convergence

periodic = 1;
nidx = 0;
kidx = 0;
factor = 0.5;
% Polynomial konvergence
K = 1;
Narray = [1,2,4,6,8,12,20,30];

errorPC = zeros(1,length(Narray));
alpha = 1;
for N = Narray
    nidx = nidx + 1;
    init
    struct_pack
    eta0 = etasol(x,0,A,omega,k);
    u0   = usol  (x,0,A,omega,k,h0);
    [u,eta] = LinSWE1D(u0,eta0,x,FinalTime,1,param,maps);
    utrue = usol(x,FinalTime,A,omega,k,h0);
    etatrue = etasol(x,FinalTime,A,omega,k);
    q = [u,eta];
    qtrue = [utrue,etatrue];

    errorPC(nidx) = max(max(abs(qtrue-q)));
%     errorP(nidx) = max(max(abs(u-utrue)));
%     
end

errorPU = zeros(1,length(Narray));
alpha = 0;
nidx = 0;
for N = Narray
    nidx = nidx + 1;
    init
    struct_pack
    eta0 = etasol(x,0,A,omega,k);
    u0   = usol  (x,0,A,omega,k,h0);
    [u,eta] = LinSWE1D(u0,eta0,x,FinalTime,1,param,maps);
    utrue = usol(x,FinalTime,A,omega,k,h0);
    etatrue = etasol(x,FinalTime,A,omega,k);
    q = [u,eta];
    qtrue = [utrue,etatrue];

    errorPU(nidx) = max(max(abs(qtrue-q)));
%     errorP(nidx) = max(max(abs(u-utrue)));
%     
end
%%
loglog(Narray,errorPU,'k-*')
hold on
loglog(Narray,errorPC,'k-o')
hold off
legend('Upwind','Central')
grid('minor')
xlabel('Polynomial Order')
ylabel('E = ||q-q_{true}||_\infty')
title('Periodic Domain')

%% Non-periodic P convergence

periodic = 0;
nidx = 0;
factor = 0.5;
% Polynomial konvergence
K = 1;
Narray = [1,2,4,6,8,12,20,30];

errorPC = zeros(1,length(Narray));
alpha = 1;
for N = Narray
    nidx = nidx + 1;
    init
    struct_pack
    eta0 = etasol(x,0,A,omega,k);
    u0   = usol  (x,0,A,omega,k,h0);
    [u,eta] = LinSWE1D(u0,eta0,x,FinalTime,1,param,maps);
    utrue = usol(x,FinalTime,A,omega,k,h0);
    etatrue = etasol(x,FinalTime,A,omega,k);
    q = [u,eta];
    qtrue = [utrue,etatrue];

    errorPC(nidx) = max(max(abs(qtrue-q)));
%     errorP(nidx) = max(max(abs(u-utrue)));
%     
end

errorPU = zeros(1,length(Narray));
alpha = 0;
nidx = 0;
for N = Narray
    nidx = nidx + 1;
    init
    struct_pack
    eta0 = etasol(x,0,A,omega,k);
    u0   = usol  (x,0,A,omega,k,h0);
    [u,eta] = LinSWE1D(u0,eta0,x,FinalTime,1,param,maps);
    utrue = usol(x,FinalTime,A,omega,k,h0);
    etatrue = etasol(x,FinalTime,A,omega,k);
    q = [u,eta];
    qtrue = [utrue,etatrue];

    errorPU(nidx) = max(max(abs(qtrue-q)));
%     errorP(nidx) = max(max(abs(u-utrue)));
%     
end
%%
loglog(Narray,errorPU,'k-*')
hold on
loglog(Narray,errorPC,'k-o')
hold off
legend('Upwind','Central')
grid('minor')
xlabel('Polynomial Order')
ylabel('E = ||q-q_{true}||_\infty')
title('Non-periodic Domain')




