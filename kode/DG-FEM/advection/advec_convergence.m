clear,clc,close all
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
% Order of polymomials used for approximation 
set(0,'defaultfigurecolor',[1 1 1])
set(groot,'DefaultAxesTickLabelInterpreter','Tex');
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


nk = [2,3,4,5,6,7];
Narray = [1,2,3];
Karray = 2.^nk;
% 
% k_line = 1./(Karray.^Narray);
tic
nidx       = 0;
kidx       = 0;
error_u_max     = zeros(length(Karray),length(Narray));
error_u_l2  = zeros(length(Karray),length(Narray));

% k_lines    = [];
alpha = 1;
for K = Karray
    kidx = kidx + 1;
    nidx = 0;
        for N = Narray
%             N
        nidx = nidx + 1;
        initialize_solver
        struct_pack; param.x = x;param.b = zeros(size(x));param.bx = zeros(size(x));
        u0  = usol(x,0,a);
        u = Advec1D(u0,x,FinalTime,0.5,param,maps);
        utrue = usol(x,FinalTime,a);
        eu    = u-utrue;
        
        erru    = zeros(1,K);
        for i = 1:K
        erru(i)    = (u(:,i)-utrue(:,i))'*Mk*(u(:,i)-utrue(:,i));       
        end
        error_u_l2(kidx,nidx)   = sqrt(sum(erru));
        error_u_max(kidx,nidx)  = max(max(abs(eu)));
   
        end
end

k_lines = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
    k_lines(:,i) = (1./Karray').^(Narray(i)+1);
end
%% plot

set(0,'defaultfigurecolor',[1 1 1])
set(groot,'DefaultAxesTickLabelInterpreter','Tex');

t = tiledlayout(1,2,'TileSpacing','compact');

% u max

k_lines_shift = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
k_lines_shift(:,i) = k_lines(:,i) * (error_u_max(1,i)/k_lines(1,i));
end

tilearray = nexttile;
plot1 = loglog(Karray,k_lines_shift,'--','Linewidth',1.5);

hold on
for i = 1:length(Narray)
    loglog(Karray,error_u_max(:,i),'^','Color',plot1(i).Color,'Linewidth',1.5);
end

% xticklabels([])
xlabel('$K$','interpreter','latex')
% xlim([2,])
ylim([1e-9,1e1])
grid minor
% ylabel('$||\textrm{error}||_{\infty}$','interpreter','latex','fontsize',14)
title('$||\textrm{error}||_{\infty}$','interpreter','latex','fontsize',14)

% u l2
k_lines_shift = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
k_lines_shift(:,i) = k_lines(:,i) * (error_u_l2(1,i)/k_lines(1,i));
end

tilearray = nexttile;
plot1 = loglog(Karray,k_lines_shift,'--','Linewidth',1.5);

hold on
for i = 1:length(Narray)
    loglog(Karray,error_u_l2(:,i),'^','Color',plot1(i).Color,'Linewidth',1.5);
end
title('$||\textrm{error}||_{L^2}$','interpreter','latex','fontsize',14)

ylim([1e-9,1e1])
grid minor
% ylabel('$||\textrm{error}||_{L^2}$','interpreter','latex','fontsize',14)
xlabel('$K$','interpreter','latex')

% title('$h(\textrm{m})$','interpreter','latex')
% xticklabels([])
% yticklabels([])




set(gcf, 'Position',  [300, 300, 1100, 500])
% set(gca,'LooseInset',get(gca,'TightInset'));
