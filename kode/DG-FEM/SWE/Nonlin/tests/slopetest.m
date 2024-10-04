%Driver script for solving the 1D wave equation using a monotone scheme
clear all,clc,close all
%Set problem parameters
L=2;FinalTime=4.0;N=2048;h=L/N;CFL=0.90;

clear,clc,close all
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\analytical_solver');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\tests')
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin')
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\nonlin\tests\';
defplotname = 'BCtest1';

% Order of polymomials used for approximation 
% load('dambreakdata.mat')
% htrue  = h;
% h0true = h0;
% hutrue = hu;
% utrue  = u;
% xtrue  = x;
% x0true = x0;

Globals1D;
NODETOL = 1e-10;
periodic = false;
BC = true;
flux = 'LF';
param.SL = 0;
dt_factor = 0.1;

N = 3;
K = 200;
FinalTime = 5;


xmin = -1;
xmax = 1;
A = 1;
g = 9.81;
b = 0;
param.b = 0;
gate = 20;
initialize_solver
hu_inflow = 1;
h_inflow = 0;
hl = 3;
hr = 3;
hul = 0;
hur = 0;
ul = hul/hl;


chalen = xmax-xmin;

% Pack structs
struct_pack

% Pack analytical parameters

[h]=wavetest(x,0.5,-0.7,0.005,10,log(2)/(36*0.005^2));
% h = sin(2*pi*x);

% plot(x,h)
type = 5;
param.SLtype = type;
[zlim,xmid,aa,bb] = SlopeLim(x,h,param);
% plot(zlim)
zlim1 = SlopeLimit1(x,h,param);
zlimN = SlopeLimitN(x,h,param);
hold on
plot(x(:),zlim(:),'b--')

plot(x(:),zlimN(:),'r--')
plot(x(:),zlim1(:),'g--')

plot(x(:),h(:),'k')

hold off
%Definedomainandinitialconditions
