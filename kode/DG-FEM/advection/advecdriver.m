% Driver script for solving the 1D advection equations
clear,clc,close all
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
% Order of polymomials used for approximation 

NODETOL = 1e-10;
alpha = 1;
N = 1;
K = 8;
% a = 2*pi;
a = 2;
f = @(x) sin(x);
periodic = 1;
FinalTime = 4;
usol = @(x,t,a) f(x-a*t);
xmin = 0.0;
xmax = 2*pi;

% Run initialization script
initialize_solver;

if ~periodic
    param.periodic = periodic;
else
    param.periodic = 1;
    vmapP(1) = vmapM(end);
    vmapP(end) = vmapM(1);
end

struct_pack

u0 = usol(x,0,a);

param.a  = a;
param.usol = usol;



%% 


%
u = Advec1D(u0,x,FinalTime,0.5,param,maps);
% hesteu = Advec1D(u,FinalTime);
% [u] = Advec1D(u,2);
% 
plot(x,u)
grid on
% hold on
% plot(x(:),usol(x(:),FinalTime,a)+1)
% xlim([xmin,xmax])
