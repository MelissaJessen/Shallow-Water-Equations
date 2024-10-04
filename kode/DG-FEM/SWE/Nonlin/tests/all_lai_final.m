clear,clc,close all
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\analytical_solver');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\tests')
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin')
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\nonlin\LaiKhan\';
defplotname = 'LaiKhanTest1';

set(groot,'DefaultAxesTickLabelInterpreter','Tex');
set(groot, 'DefaultTextInterpreter', 'Tex')
set(groot, 'DefaultLegendInterpreter', 'Tex')
% Order of polymomials used for approximation 
% load('dambreakdata.mat')
% htrue  = h;
% h0true = h0;
% hutrue = hu;
% utrue  = u;
% xtrue  = x;
% x0true = x0;
NODETOL = 1e-10;

periodic = false;
BC = 1;
flux = @Roe_swe;
% flux = @LaxFswe;

param.SL = 1;
dt_factor = 1;
param.SLtype = 2;
param.Limiter = @SlopeLim;

N = 2;
K = 100;
FinalTime = 25;

xmin = 0;
xmax = 25;
A = 1;
g = 9.81;

gate = 20;

initialize_solver
hl = 0.33;
hr = 0.33;
ul = 0;
ur = 0;
% u_inflow = 0.18/hl;
h_inflow = @(x,t,h,hu,param,maps) h(maps.vmapI);
hu_inflow =  @(x,t,h,hu,param,maps) 0;

h_outflow = @(x,t,h,hu,param,maps) h(maps.vmapO);
hu_outflow =  @(x,t,h,hu,param,maps) 0;
u_inflow = 0;
u_outflow = 0;
% u_inflow = 0;
h  = zeros(size(x));
u = zeros(size(x));

h(x<=gate) = hl;
h(x>gate) = hr;


u(x<=gate) = ul;
u(x>gate) = ur;

% parabolic bump;
alpha = 0.05;


% [minValue,a8] = min(abs(x(:)-8))
b = 0.2 - alpha*(x-10).^2;
a8 = find(x==8);
a12 = find(x==12);
b(x<8) = 0;
b(12<x) = 0;


h = h-b;% Det virker med at trække b fra


bx = -2*alpha*(x-10);
bx(x<8) = 0;
bx(12<x) = 0;
bx(a8(1)) = 0;
bx(a12(2)) = 0;

% 

hu = u.*h;

% Pack structs
struct_pack
param.g = g;
param.xmin   = xmin;
param.xmax   = xmax;
param.BC     = BC;
param.hl     = hl;
param.hr     = hr;
param.ul     = ul;
param.ur     = ur;
param.flux   = flux;
param.bx     = bx;
param.Minv   = Minv;
param.u_inflow = u_inflow;
param.h_inflow = h_inflow;
param.hu_inflow = hu_inflow;
param.h_outflow = h_outflow;
param.hu_outflow = hu_outflow;
param.u_outflow = 0;
param.Limiter = @SlopeLim;
param.SLtype = 2;
param.b = b;
param.x = x;


% Pack analytical parameters

nrtol  = 1e-10;
niter = 1e3;
cl = sqrt(g*hl);
cr = sqrt(g*hr);
% alpha = 0.0*pi;
mcells = length(x(:));
chalen = xmax;
timout = FinalTime;

%% Lai & Khan test 1

h0 = h;
hu0 = hu;
flux = @LaxFswe;
param.flux = flux;
% [h_lf1,hu_lf1,ht_lf1,hut_lf1,tt_lf1] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);
[h_lf1,hu_lf1] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @HLLswe;
param.flux = flux;
% [h_hll1,hu_hll1,ht_hll1,hut_hll1,tt_hll1] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);
[h_hll1,hu_hll1] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @HLLCswe;
param.flux = flux;
% [h_hllc1,hu_hllc1,ht_hllc1,hut_hllc1,tt_hllc1] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);
[h_hllc1,hu_hllc1] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @Roe_swe;
param.flux = flux;
% [h_roe1,hu_roe1,ht_roe1,hut_roe1,tt_roe1] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);
[h_roe1,hu_roe1] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

%% Init test 2

hl = 0.5;
hr = 0.5;
hu_inflow_cst = 0.18;

h_inflow = @(x,t,h,hu,param,maps) h(maps.vmapI);
hu_inflow =  @(x,t,h,hu,param,maps) 0.18;

h_outflow = @(x,t,h,hu,param,maps) h(maps.vmapO);
hu_outflow =  @(x,t,h,hu,param,maps) hu(maps.vmapO);
ul = hu_inflow_cst/hl;
ur = hu_inflow_cst/hl;
u_inflow = 0.18/hl;


u_outflow = 0;
h  = zeros(size(x));
hu = zeros(size(x));

initialize_solver
param.u_inflow = u_inflow;
param.u_outflow = 0;
param.u_inflow = u_inflow;
param.h_inflow = h_inflow;
param.hu_inflow = hu_inflow;
param.h_outflow = h_outflow;
param.hu_outflow = hu_outflow;


hu(x<=gate) = hu_inflow_cst;
hu(x>gate) = hu_inflow_cst;

h(x<=gate) = hl;
h(x>gate) = hr;


% u(x<=gate) = ul;
% u(x>gate) = ur;


h = h-b;

%% test 2 compute
h0 = h;
hu0 = hu;
flux = @LaxFswe;
param.flux = flux;
[h_lf2,hu_lf2] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @HLLswe;
param.flux = flux;
[h_hll2,hu_hll2] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @HLLCswe;
param.flux = flux;
[h_hllc2,hu_hllc2] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @Roe_swe;
param.flux = flux;
[h_roe2,hu_roe2] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);


%% Init test 3

hl = 0.33;
hr = 0.33;
hu_inflow_cst = 0.18;

h_inflow = @(x,t,h,hu,param,maps) h(maps.vmapI);
hu_inflow =  @(x,t,h,hu,param,maps) 0.18;

h_outflow = @(x,t,h,hu,param,maps) h(maps.vmapO);
hu_outflow =  @(x,t,h,hu,param,maps) hu(maps.vmapO);
u_inflow = 0.18/hl;


h  = zeros(size(x));
hu = zeros(size(x));

initialize_solver
param.u_inflow = u_inflow;
param.u_outflow = 0;
param.u_inflow = u_inflow;
param.h_inflow = h_inflow;
param.hu_inflow = hu_inflow;
param.h_outflow = h_outflow;
param.hu_outflow = hu_outflow;


hu(x<=gate) = hu_inflow_cst;
hu(x>gate) = hu_inflow_cst;

h(x<=gate) = hl;
h(x>gate) = hr;


h = h-b;

%% Test 3 compute

h0 = h;
hu0 = hu;
flux = @LaxFswe;
param.flux = flux;
[h_lf3,hu_lf3] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @HLLswe;
param.flux = flux;
[h_hll3,hu_hll3] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @HLLCswe;
param.flux = flux;
[h_hllc3,hu_hllc3] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @Roe_swe;
param.flux = flux;
[h_roe3,hu_roe3] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

%% Test 4 init


hl = 2;
hr = 2;
hu_inflow_cst = 25.0567;

h_inflow = @(x,t,h,hu,param,maps) h(maps.vmapI);
hu_inflow =  @(x,t,h,hu,param,maps) 25.0567;

h_outflow = @(x,t,h,hu,param,maps) h(maps.vmapO);
hu_outflow =  @(x,t,h,hu,param,maps) hu(maps.vmapO);
ul = hu_inflow_cst/hl;
ur = hu_inflow_cst/hl;
u_inflow = 25.0567/hl;


u_outflow = 0;
h  = zeros(size(x));
hu = zeros(size(x));

initialize_solver
param.u_inflow = u_inflow;
param.u_outflow = 0;
param.u_inflow = u_inflow;
param.h_inflow = h_inflow;
param.hu_inflow = hu_inflow;
param.h_outflow = h_outflow;
param.hu_outflow = hu_outflow;


hu(x<=gate) = hu_inflow_cst;
hu(x>gate) = hu_inflow_cst;

h(x<=gate) = hl;
h(x>gate) = hr;


% u(x<=gate) = ul;
% u(x>gate) = ur;


h = h-b;


%%
hl = 2;
hr = 2;

hu_inflow_cst = 25.0567;

h_inflow = @(x,t,h,hu,param,maps) h(maps.vmapI);
hu_inflow =  @(x,t,h,hu,param,maps) hu_inflow_cst;

h_outflow = @(x,t,h,hu,param,maps) h(maps.vmapO);
hu_outflow =  @(x,t,h,hu,param,maps) hu(maps.vmapO);
ul = hu_inflow_cst/hl;
ur = hu_inflow_cst/hl;
u_inflow = 0.18/hl;


u_outflow = 0;
h  = zeros(size(x));
hu = zeros(size(x));

initialize_solver
param.u_inflow = u_inflow;
param.u_outflow = 0;
param.u_inflow = u_inflow;
param.h_inflow = h_inflow;
param.hu_inflow = hu_inflow;
param.h_outflow = h_outflow;
param.hu_outflow = hu_outflow;
param.b = b;

hu(x<=gate) = hu_inflow_cst;
hu(x>gate) = hu_inflow_cst;

h(x<=gate) = hl;
h(x>gate) = hr;


% u(x<=gate) = ul;
% u(x>gate) = ur;


% flow_factor = h./(h-b);
h = h-b;

u = hu./h;

% Pack structs
struct_pack
param.g = g;
param.xmin   = xmin;
param.xmax   = xmax;
param.BC     = BC;
param.hl     = hl;
param.hr     = hr;
param.ul     = ul;
param.ur     = ur;
param.flux   = flux;
param.bx     = bx;
param.Minv   = Minv;
param.x = x;


% Pack analytical parameters

nrtol  = 1e-10;
niter = 1e3;
cl = sqrt(g*hl);
cr = sqrt(g*hr);
% alpha = 0.0*pi;
mcells = length(x(:));
chalen = xmax;
timout = FinalTime;
%%


%% test 4 compute
h0 = h;
hu0 = hu;
flux = @LaxFswe;
param.flux = flux;
[h_lf4,hu_lf4] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @HLLswe;
param.flux = flux;
[h_hll4,hu_hll4] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @HLLCswe;
param.flux = flux;
[h_hllc4,hu_hllc4] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @Roe_swe;
param.flux = flux;
[h_roe4,hu_roe4] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

%% h with tile plots

figure
t = tiledlayout(2,2,'TileSpacing','compact');
h1 = nexttile;
plot(x(:),h_lf1(:)+b(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),h_hll1(:)+b(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_hllc1(:)+b(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_roe1(:)+b(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),b(:),'k--','Linewidth',1.5)
xlim([0,25])
% axes_factor = 0.01;
% xlims =  get(gca,'Xlim');
% xlim([xlims(1)-axes_factor*xlims(2),xlims(2) + axes_factor*xlims(2)])

% pos = [7,0.31,13-7,0.35-0.31];
% rectangle('Position',pos);
% a1 = axes('Position',[0,0,pos(3)/25,pos(4)/0.4])
hold off
grid minor
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)
% xticklabels([])
ylim([0,0.6])


h2 = nexttile;

plot(x(:),h_lf2(:)+b(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),h_hll2(:)+b(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_hllc2(:)+b(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_roe2(:)+b(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),b(:),'k--','Linewidth',1.5)
ylim([0,0.6])
hold off
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)



grid minor
% yticklabels([])
% xticklabels([])
axes_factor = 0.01;
xlim([0,25])
xticks(0:5:25)

xlims =  get(gca,'Xlim');
xlim([xlims(1)-axes_factor*xlims(2),xlims(2) + axes_factor*xlims(2)])


h3 = nexttile;
plot(x(:),h_lf3(:)+b(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),h_hll3(:)+b(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_hllc3(:)+b(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_roe3(:)+b(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),b(:),'k--','Linewidth',1.5)

hold off
grid minor
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)
% xticklabels([])
ylim([0,0.6])

axes_factor = 0.01;
xlim([0,25])

xlims =  get(gca,'Xlim');
xlim([xlims(1)-axes_factor*xlims(2),xlims(2) + axes_factor*xlims(2)])
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)
xticks(0:5:25)


h4 = nexttile;
plot(x(:),h_lf4(:)+b(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),h_hll4(:)+b(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_hllc4(:)+b(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_roe4(:)+b(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),b(:),'k--','Linewidth',1.5)

hold off
% yticklabels([])
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)
axes_factor = 0.01;
xlim([0,25])
xticks(0:5:25)

xlims =  get(gca,'Xlim');
xlim([xlims(1)-axes_factor*xlims(2),xlims(2) + axes_factor*xlims(2)])
ylim([0,2.5])
grid minor
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)

hold on
LFleg = plot(nan, nan, 'ko','DisplayName','LF','MarkerSize',6,'Linewidth',1.5);
HLLleg = plot(nan, nan, 'ks','DisplayName','HLL','MarkerSize',6,'Linewidth',1.5);
HLLCleg = plot(nan, nan, 'k^','DisplayName','HLLC','MarkerSize',6,'Linewidth',1.5);
Roeleg = plot(nan, nan, 'k+','DisplayName','Roe','MarkerSize',6,'Linewidth',1.5);
bleg = plot(nan, nan, 'k--','DisplayName','Bottom','Linewidth',1.5,'Linewidth',1.5);
lgd = legend([LFleg,HLLleg,HLLCleg,Roeleg,bleg],'interpreter','latex');
lgd.FontSize = 9;
lgd.Layout.Tile = 'south';
lgd.Orientation = 'horizontal';
% set(lgd.Position,[0.6520 0.1405 0.1387 0.1730])

hold off
set(gcf, 'Position',  [300, 300, 600, 300])
set(gca,'LooseInset',get(gca,'TightInset'));

%% hu plots
figure
t = tiledlayout(2,2,'TileSpacing','compact');

h1 = nexttile;
plot(x(:),hu_lf1(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),hu_hll1(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_hllc1(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_roe1(:), 'k+','MarkerSize',4,'Linewidth',1.5)
% plot(x(:),b(:),'k--','Linewidth',1.5)   
xlim([0,25])
xlims =  get(gca,'Xlim');
xlim([xlims(1)-axes_factor*xlims(2),xlims(2) + axes_factor*xlims(2)])
% pos = [7,0.31,13-7,0.35-0.31];
% rectangle('Position',pos);
% a1 = axes('Position',[0,0,pos(3)/25,pos(4)/0.4])
hold off
grid minor
ylabel('$hu[\mathrm m^2/s]$','interpreter','latex','FontSize',11)
% % xticklabels([])
% ylim([0,0.4])
xticks(0:5:25)


h2 = nexttile;
plot(x(:),hu_lf2(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),hu_hll2(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_hllc2(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_roe2(:), 'k+','MarkerSize',4,'Linewidth',1.5)
% plot(x(:),b(:),'k--','Linewidth',1.5)
hold off
% ylim([0,0.6])
grid minor
% yticklabels([])
% xticklabels([])
ylabel('$hu[\mathrm m^2/s]$','interpreter','latex','FontSize',11)

xticks(0:5:25)
xlim([0,25])
xlims =  get(gca,'Xlim');
xlim([xlims(1)-axes_factor*xlims(2),xlims(2) + axes_factor*xlims(2)])
h3 = nexttile;
plot(x(:),hu_lf3(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),hu_hll3(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_hllc3(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_roe3(:), 'k+','MarkerSize',4,'Linewidth',1.5)
% plot(x(:),b(:),'k--','Linewidth',1.5)

hold off
grid minor
ylabel('$hu[\mathrm m^2/s]$','interpreter','latex','FontSize',11)
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)
xticks(0:5:25)

% xticklabels([])
% ylim([0,0.5])
xlim([0,25])

xlims =  get(gca,'Xlim');
xlim([xlims(1)-axes_factor*xlims(2),xlims(2) + axes_factor*xlims(2)])
h4 = nexttile;
plot(x(:),hu_lf4(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),hu_hll4(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_hllc4(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_roe4(:), 'k+','MarkerSize',4,'Linewidth',1.5)
% plot(x(:),b(:),'k--','Linewidth',1.5)
xlim([0,25])
xlims =  get(gca,'Xlim');
xlim([xlims(1)-axes_factor*xlims(2),xlims(2) + axes_factor*xlims(2)])
hold off
% yticklabels([])
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)

grid minor
ylabel('$hu[\mathrm m^2/s]$','interpreter','latex','FontSize',11)

hold on
LFleg = plot(nan, nan, 'ko','DisplayName','LF','MarkerSize',6,'Linewidth',1.5);
HLLleg = plot(nan, nan, 'ks','DisplayName','HLL','MarkerSize',6,'Linewidth',1.5);
HLLCleg = plot(nan, nan, 'k^','DisplayName','HLLC','MarkerSize',6,'Linewidth',1.5);
Roeleg = plot(nan, nan, 'k+','DisplayName','Roe','MarkerSize',6,'Linewidth',1.5);
lgd = legend([LFleg,HLLleg,HLLCleg,Roeleg],'interpreter','latex');
lgd.Layout.Tile = 'south';
lgd.Orientation = 'horizontal';

lgd.FontSize = 9;
set(gcf, 'Position',  [300, 300, 600, 300])
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)
xticks(0:5:25)


%%
figure

t = tiledlayout(2,2,'TileSpacing','compact');

h1 = nexttile;
plot(x(:),h_lf1(:)+b(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),h_hll1(:)+b(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_hllc1(:)+b(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_roe1(:)+b(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),b(:),'k--','Linewidth',1.5)   
xlim([7,13])
ylim([0.3,0.35])

% pos = [7,0.31,13-7,0.35-0.31];
% rectangle('Position',pos);
% a1 = axes('Position',[0,0,pos(3)/25,pos(4)/0.4])
hold off
grid minor
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)
% xticklabels([])
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)


h2 = nexttile;

plot(x(:),h_lf2(:)+b(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),h_hll2(:)+b(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_hllc2(:)+b(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_roe2(:)+b(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),b(:),'k--','Linewidth',1.5)
hold off
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)

xlim([7,13])
ylim([0.45,0.55])

grid minor
% yticklabels([])
% xticklabels([])


h3 = nexttile;
plot(x(:),h_lf3(:)+b(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),h_hll3(:)+b(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_hllc3(:)+b(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_roe3(:)+b(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),b(:),'k--','Linewidth',1.5)

hold off
grid minor
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)
% xticklabels([])
xlim([11,17])
ylim([0.16,0.34])
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)


h4 = nexttile;
plot(x(:),h_lf4(:)+b(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),h_hll4(:)+b(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_hllc4(:)+b(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_roe4(:)+b(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),b(:),'k--','Linewidth',1.5)

hold off
% yticklabels([])
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)
xlim([7,13])
ylim([1.9,2.3])
grid minor
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)

hold on
LFleg = plot(nan, nan, 'ko','DisplayName','LF','MarkerSize',6,'Linewidth',1.5);
HLLleg = plot(nan, nan, 'ks','DisplayName','HLL','MarkerSize',6,'Linewidth',1.5);
HLLCleg = plot(nan, nan, 'k^','DisplayName','HLLC','MarkerSize',6,'Linewidth',1.5);
Roeleg = plot(nan, nan, 'k+','DisplayName','Roe','MarkerSize',6,'Linewidth',1.5);
bleg = plot(nan, nan, 'k--','DisplayName','Bottom','Linewidth',1.5);
lgd = legend([LFleg,HLLleg,HLLCleg,Roeleg]);
lgd.Layout.Tile = 'south';
lgd.FontSize = 11;
lgd.Orientation = 'horizontal';

% set(lgd.Position,[0.6520 0.1405 0.1387 0.1730])

hold off
set(gcf, 'Position',  [300, 300, 600, 300])
set(gca,'LooseInset',get(gca,'TightInset'));



%% Initial conditions

figure

t = tiledlayout(2,2,'TileSpacing','compact');

h1 = nexttile;
plot(x(:),ones(size(x(:)))*0.33, 'k--','MarkerSize',4,'Linewidth',1.5)
hold on
plot(x(:),b(:),'k-','Linewidth',1.5)   
xlim([0,25])
ylim([0,0.6])
title('Parabolic bump case 1')

% pos = [7,0.31,13-7,0.35-0.31];
% rectangle('Position',pos);
% a1 = axes('Position',[0,0,pos(3)/25,pos(4)/0.4])
hold off
grid minor
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)
% xticklabels([])
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)


h2 = nexttile;

plot(x(:),ones(size(x(:)))*0.5, 'k--','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),b(:),'k-','Linewidth',1.5)
hold off
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)

xlim([0,25])
ylim([0,0.6])
title('Parabolic bump case 2')

grid minor
% yticklabels([])
% xticklabels([])


h3 = nexttile;
plot(x(:),ones(size(x(:)))*0.33, 'k--','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),b(:),'k-','Linewidth',1.5)

hold off
grid minor
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)
% xticklabels([])
xlim([0,25])
ylim([0,0.6])
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)
title('Parabolic bump case 2')


h4 = nexttile;
plot(x(:),ones(size(x(:)))*2, 'k--','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),b(:),'k-','Linewidth',1.5)

hold off
% yticklabels([])
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)
xlim([0,25])
ylim([0,2.5])
grid minor
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)

hold on
hleg = plot(nan, nan, 'k-','DisplayName','$h_0$','MarkerSize',6,'Linewidth',1.5);
bleg = plot(nan, nan, 'k--','DisplayName','Bottom','Linewidth',1.5);
lgd = legend([hleg,bleg],'interpreter','latex');
lgd.Layout.Tile = 'south';
lgd.FontSize = 11;
lgd.Orientation = 'horizontal';

% set(lgd.Position,[0.6520 0.1405 0.1387 0.1730])
title('Parabolic bump case 4')

hold off
set(gcf, 'Position',  [300, 300, 600, 300])
set(gca,'LooseInset',get(gca,'TightInset'));


