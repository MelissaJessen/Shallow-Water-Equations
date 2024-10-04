clear,clc,close all
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
% Order of polymomials used for approximation 

NODETOL = 1e-10;
alpha = 0;
P = [1,2,3,4];


Narray = [2,4,6,8,12,16,24];
Karray = [2,4,6,8,12,16,24];
% a = 2*pi;
a = 2;
f = @(x) sin(x);
usol = @(x,t,a) f(x-a*t);
xmin = 0.0;
xmax = 4.0;
xplot = linspace(xmin,xmax);
% Generate simple mesh
%%
nidx = 0;
kidx = 0;
factor = 0.5;
% Polynomial konvergence
K = 1;
errorP = zeros(1,length(Narray));
for N = Narray
    nidx = nidx + 1;
    init
    struct_pack
    u = Advec1D(u0,x,FinalTime,factor,param,maps);
    utrue = usol(x,FinalTime,a);
    errorP(nidx) = norm(utrue-u,'inf');
end
 
loglog(Narray,errorP,'-*')
grid('minor')

%% K convergence
errorK = zeros(3,length(Karray));
kidx = 0;
factor = 0.5;
for K = Karray
    kidx = kidx +1;
    nidx = 0;
    for N = [2,4,8]
        nidx = nidx + 1;
        init
        struct_pack
        u = Advec1D(u0,x,FinalTime,factor,param,maps);
        utrue = usol(x,FinalTime,a);
        errorK(nidx,kidx) = norm(utrue-u,'inf');
    end
end
 %%
 symbs = {'k-*','k-o','k-p'};
for i = 1:3
loglog(Karray,errorK(i,:),symbs{i})
hold on
end

hold off

grid minor