clear,clc,close all
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\Linear';
% Order of polymomials used for approximation 

NODETOL = 1e-10;
FinalTime = 2;

N = 6;
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


%% Periodic K convergence

periodic = 1;
factor = 0.5;
% Polynomial konvergence
nk = [1,2,3,4,5,6];

Karray = 2.^nk;

nidx = 0;
errorKCper = zeros(1,length(Karray));
alpha = 1;
for K = Karray
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

    errorKCper(nidx) = max(max(abs(qtrue-q)));
%     errorP(nidx) = max(max(abs(u-utrue)));
%     
end

errorKUper = zeros(1,length(Karray));
alpha = 0;
nidx = 0;
for K = Karray
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

    errorKUper(nidx) = max(max(abs(qtrue-q)));
%     errorP(nidx) = max(max(abs(u-utrue)));
%     
end
%%

loglog(Karray,errorKUper,'k-*')
hold on
loglog(Karray,errorKCper,'k-o')
loglog(Karray,k_line,'k--')
hold off
legend('Upwind','Central')
grid('minor')
xlabel('K')
ylabel('E = ||q-q_{true}||_\infty')
title('Periodic Domain')


%% Non-periodic K convergence

periodic = 0;
factor = 0.5;
% Polynomial konvergence

nidx = 0;
errorKC = zeros(1,length(Karray));
alpha = 1;
for K = Karray
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

    errorKC(nidx) = max(max(abs(qtrue-q)));
%     errorP(nidx) = max(max(abs(u-utrue)));
%     
end

errorKU = zeros(1,length(Karray));
alpha = 0;
nidx = 0;
for K = Karray
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

    errorKU(nidx) = max(max(abs(qtrue-q)));
%     errorP(nidx) = max(max(abs(u-utrue)));
%     
end
%%
figure
loglog(Karray,errorKU,'k-*')
hold on
loglog(Karray,errorKC,'k-o')
loglog(Karray,k_line,'k--')
hold off
legend('Upwind','Central')
grid('minor')
xlabel('K')
ylabel('E = ||q-q_{true}||_\infty')
title('Non-periodic Domain')
%%
figure
loglog(Karray,errorKU,'k-*')
hold on
loglog(Karray,errorKC,'k-o')
loglog(Karray,errorKUper,'r-*')
loglog(Karray,errorKCper,'r-o')
loglog(Karray,k_line,'b--')
hold off
grid('minor')
legend('Upwind','Central','Upwind Periodic','Central Periodic','1/K^N')

xlabel('K')
ylabel('E = ||q-q_{true}||_\infty')
