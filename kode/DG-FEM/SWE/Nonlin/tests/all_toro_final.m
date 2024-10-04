
clear,clc
% close all
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\analytical_solver');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin')
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\nonlin\toro-cases\';
defplotname = 'torotest1';
set(0,'defaultfigurecolor',[1 1 1])
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
BC = false;


param.SL = 1;
param.SLtype = 2;
dt_factor = 1;


N = 2;
K = 100;
FinalTime = 2.5;


xmin = 0;
xmax = 50;
A = 1;
g = 9.81;

gate = 20;
initialize_solver
hl = 1;
hr = 0.1;
hlt = hl;
hrt = hr;
ul = 2.5;
ur = 0;  
ult = ul;
urt = ur;
u_inflow = 0;
h  = zeros(size(x));
u = zeros(size(x));

h(x<=gate) = hl;
h(x>gate) = hr;


u(x<=gate) = ul;
u(x>gate) = ur;
hu = u.*h;


% % No slope
alpha = 0;
bx = 0;
b = 0;

% Linear slope
% alpha = 0;
% bx = ones(size(x))*alpha;

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
% param.flux   = flux;
param.bx     = bx;
param.Minv   = Minv;
param.u_inflow = u_inflow;
param.Limiter = @SlopeLimit1;
param.b = b;



param.EToV = EToV;
% Pack analytical parameters

nrtol  = 1e-10;
niter = 1e3;
cl = sqrt(g*hl);
cr = sqrt(g*hr);
mcells = length(x(:));
chalen = xmax;
timout = FinalTime;

pack_analytical_parameters
dcrit = (ur - ul) - 2.0*(cl + cr);
h0 = h;
hu0 = hu;


if(hl<=0 || hr<=0 || dcrit>=0)
%     disp(['Dry',' BC: ',num2str(BC),' flux: ',flux])
    [htrue1,utrue1] = drybed(h,u,analytical_param);
%     plot_title = ['Dry',' x_0=',num2str(gate),' T=',num2str(FinalTime),' h_L=',num2str(hl),' h_R=',num2str(hr),' u_L=',num2str(ul),' u_R=',num2str(ur), ' \alpha = ',num2str(alpha),', t = ',num2str(FinalTime)];

else
%     disp(['Wet',' BC: ',num2str(BC),' flux: ',flux])
    [htrue1,utrue1,converged_flag] = wetbed(h,u,niter,analytical_param);
%     plot_title = ['Wet',' x_0=',num2str(gate),' T=',num2str(FinalTime),' h_L=',num2str(hl),' h_R=',num2str(hr),' u_L=',num2str(ul),' u_R=',num2str(ur), ' \alpha = ',num2str(alpha),', t = ',num2str(FinalTime)];
end

%% Toro 1
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



%% Toro 2 Init
xmin = 0;
xmax = 50;
A = 1;
g = 9.81;
% xplot(xmin,xmax,1000);

gate = 25;
initialize_solver
hl = 1;
hr = 1;
hlt = hl;
hrt = hr;
ul = -1;
ur = 1;  
ult = ul;
urt = ur;
u_inflow = 0;
h  = zeros(size(x));
u = zeros(size(x));

h(x<=gate) = hl;
h(x>gate) = hr;


u(x<=gate) = ul;
u(x>gate) = ur;
hu = u.*h;


% % No slope
alpha = 0;
bx = 0;


% Linear slope
% alpha = 0;
% bx = ones(size(x))*alpha;

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


% Pack analytical parameters

nrtol  = 1e-10;
niter = 1e3;
cl = sqrt(g*hl);
cr = sqrt(g*hr);
mcells = length(x(:));
chalen = xmax;
timout = FinalTime;

pack_analytical_parameters
dcrit = (ur - ul) - 2.0*(cl + cr);



if(hl<=0 || hr<=0 || dcrit>=0)
    [htrue2,utrue2] = drybed(h,u,analytical_param);
else
    [htrue2,utrue2,converged_flag] = wetbed(h,u,niter,analytical_param);
end
%% Toro 2 Compute
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


%% Toro 3 Init


gate = 20;
initialize_solver
hl = 1;
hr = 0.001;
hlt = 1;
hrt = 0;
ul = 0;
ult = 0;
ur = 0; 
urt = 0;
u_inflow = 0;
h  = zeros(size(x));
h_ana = zeros(size(x));
u = zeros(size(x));
u_ana = zeros(size(x));

h(x<=gate) = hl;
h(x>gate) = hr;
h_ana(x<=gate) = hlt;
h_ana(x>gate) = hrt;


u(x<=gate) = ul;
u(x>gate) = ur;
hu = u.*h;


% % No slope
alpha = 0;
bx = 0;


% Linear slope
% alpha = 0;
% bx = ones(size(x))*alpha;

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


% Pack analytical parameters

nrtol  = 1e-10;
niter = 1e3;
cl = sqrt(g*hlt);
cr = sqrt(g*hrt);
mcells = length(x(:));
chalen = xmax;
timout = FinalTime;

pack_analytical_parameters
dcrit = (urt - ult) - 2.0*(cl + cr);



if(hlt<=0 || hrt<=0 || dcrit>=0)
    [htrue3,utrue3] = drybed(h_ana,u,analytical_param);
else
    [htrue3,utrue3,converged_flag] = wetbed(h_ana,u,niter,analytical_param);
end
%% Toro 3 Compute
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

%% Toro 4 Init


gate = 30;
initialize_solver
hl = 0.005;
hr = 1;
hlt = 0;
hrt = 1;
ul = 0;
ult = 0;
ur = 0; 
urt = 0;
u_inflow = 0;
h  = zeros(size(x));
h_ana = zeros(size(x));
u = zeros(size(x));
u_ana = zeros(size(x));

h(x<=gate) = hl;
h(x>gate) = hr;
h_ana(x<=gate) = hlt;
h_ana(x>gate) = hrt;


u(x<=gate) = ul;
u(x>gate) = ur;
hu = u.*h;


% % No slope
alpha = 0;
bx = 0;


% Linear slope
% alpha = 0;
% bx = ones(size(x))*alpha;

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


% Pack analytical parameters

nrtol  = 1e-10;
niter = 1e3;
cl = sqrt(g*hlt);
cr = sqrt(g*hrt);
mcells = length(x(:));
chalen = xmax;
timout = FinalTime;

pack_analytical_parameters
dcrit = (urt - ult) - 2.0*(cl + cr);



if(hlt<=0 || hrt<=0 || dcrit>=0)
    [htrue4,utrue4] = drybed(h_ana,u,analytical_param);
else
    [htrue4,utrue4,converged_flag] = wetbed(h_ana,u,niter,analytical_param);
end

%% Toro 4 compute
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

%% Toro 5 Init


gate = 25;
initialize_solver
hl = 0.1;
hr = 0.1;
hlt = hl;
hrt = hr;
ul = -3;
ur = 3; 
ult = ul;
urt = ur;
u_inflow = 0;
h  = zeros(size(x));
u = zeros(size(x));

h(x<=gate) = hl;
h(x>gate) = hr;


u(x<=gate) = ul;
u(x>gate) = ur;
hu = u.*h;


% % No slope
alpha = 0;
bx = 0;


% Linear slope
% alpha = 0;
% bx = ones(size(x))*alpha;

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


% Pack analytical parameters

nrtol  = 1e-10;
niter = 1e3;
cl = sqrt(g*hl);
cr = sqrt(g*hr);
mcells = length(x(:));
chalen = xmax;
timout = FinalTime;

pack_analytical_parameters
dcrit = (ur - ul) - 2.0*(cl + cr);



if(hl<=0 || hr<=0 || dcrit>=0)

    [htrue5,utrue5] = drybed(h,u,analytical_param);

else
    [htrue5,utrue5,converged_flag] = wetbed(h,u,niter,analytical_param);
end
%% Toro 5 Compute
h0 = h;
hu0 = hu;
flux = @LaxFswe;
param.flux = flux;
[h_lf5,hu_lf5] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @HLLswe;
param.flux = flux;
[h_hll5,hu_hll5] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @HLLCswe;
param.flux = flux;
[h_hllc5,hu_hllc5] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

flux = @Roe_swe;
param.flux = flux;
[h_roe5,hu_roe5] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

%% Plots
t = tiledlayout(3,2,'TileSpacing','tight');

h1 = nexttile;
plot(x(:),h_lf1(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),h_hll1(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_hllc1(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_roe1(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),htrue1(:),'r--','Linewidth',1.5)
hold off
grid minor
ylabel('$h[\mathrm m]$','interpreter','latex','FontSize',11)
% xticklabels([])

h2 = nexttile;
plot(x(:),h_lf2(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),h_hll2(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_hllc2(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_roe2(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),htrue2(:),'r--','Linewidth',1.5)
hold off
grid minor
ylabel('$h[\mathrm m]$','interpreter','latex','FontSize',11)
% xticklabels([])



h3 = nexttile;
plot(x(:),h_lf3(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),h_hll3(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_hllc3(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_roe3(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),htrue3(:),'r--','Linewidth',1.5)
hold off
grid minor
ylabel('$h[\mathrm m]$','interpreter','latex','FontSize',11)
% xticklabels([])


h4 = nexttile;
plot(x(:),h_lf4(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),h_hll4(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_hllc4(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_roe4(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),htrue4(:),'r--','Linewidth',1.5)
hold off
ylabel('$h[\mathrm m]$','interpreter','latex','FontSize',11)
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)

grid minor

h5 = nexttile;
plot(x(:),h_lf5(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),h_hll5(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_hllc5(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),h_roe5(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot1 = plot(x(:),htrue5(:),'r--','Linewidth',1.5,'DisplayName','True');

ylabel('$h[\mathrm m]$','interpreter','latex','FontSize',11)

hold off
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)

hold on
LFleg = plot(nan, nan, 'ko','DisplayName','LF','MarkerSize',6,'Linewidth',1.5);
HLLleg = plot(nan, nan, 'ks','DisplayName','HLL','MarkerSize',6,'Linewidth',1.5);
HLLCleg = plot(nan, nan, 'k^','DisplayName','HLLC','MarkerSize',6,'Linewidth',1.5);
Roeleg = plot(nan, nan, 'k+','DisplayName','Roe','MarkerSize',6,'Linewidth',1.5);
Trueleg = plot(nan, nan, 'r--','DisplayName','True','Linewidth',1.5);
lgd = legend([plot1,LFleg,HLLleg,HLLCleg,Roeleg],'interpreter','latex');
% lgd.Layout.Tile = 'southeast';
lgd.FontSize = 9;
% set(lgd.Position,[0.6520 0.1405 0.1387 0.1730])

hold off
grid minor
set(gcf, 'Position',  [300, 300, 600, 300])
set(gca,'LooseInset',get(gca,'TightInset'));
% lgd.Layout.Tile = 'west';
%%

%% Plots hu
figure
t = tiledlayout(3,2,'TileSpacing','tight');

h1 = nexttile;
plot(x(:),hu_lf1(:)./h_lf1(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),hu_hll1(:)./h_hll1(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_hllc1(:)./h_hllc1(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_roe1(:)./h_roe1(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),utrue1(:),'r--','Linewidth',1.5)
hold off
grid minor
ylabel('$u[\mathrm m/\mathrm s]$','interpreter','latex','FontSize',11)
xticklabels([])

h2 = nexttile;
plot(x(:),hu_lf2(:)./h_lf2(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),hu_hll2(:)./h_hll2(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_hllc2(:)./h_hllc2(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_roe2(:)./h_roe2(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),utrue2(:),'r--','Linewidth',1.5)
hold off
grid minor
ylabel('$u[\mathrm m/\mathrm s]$','interpreter','latex','FontSize',11)
xticklabels([])


h3 = nexttile;
plot(x(:),hu_lf3(:)./h_lf3(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),hu_hll3(:)./h_hll3(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_hllc3(:)./h_hllc3(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_roe3(:)./h_roe3(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),utrue3(:),'r--','Linewidth',1.5)
hold off
grid minor
xticklabels([])

ylabel('$u[\mathrm m/\mathrm s]$','interpreter','latex','FontSize',11)



h4 = nexttile;

plot(x(:),hu_lf4(:)./h_lf4(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),hu_hll4(:)./h_hll4(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_hllc4(:)./h_hllc4(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_roe4(:)./h_roe4(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot(x(:),utrue4(:),'r--','Linewidth',1.5)
hold off
grid minor
ylabel('$u[\mathrm m/\mathrm s]$','interpreter','latex','FontSize',11)
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)

% ylabel('$h(\mathrm m)$','interpreter','latex')

h5 = nexttile;

plot(x(:),hu_lf5(:)./h_lf5(:), 'ko','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),hu_hll5(:)./h_hll5(:), 'ks','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_hllc5(:)./h_hllc5(:), 'k^','MarkerSize',4,'Linewidth',1.5)
plot(x(:),hu_roe5(:)./h_roe5(:), 'k+','MarkerSize',4,'Linewidth',1.5)
plot1 = plot(x(:),utrue5(:),'r--','Linewidth',1.5,'DisplayName','True');
hold off
grid minor
% yticklabels([])
ylabel('$u[\mathrm m/\mathrm s]$','interpreter','latex','FontSize',11)
xlabel('$x[\mathrm m]$','interpreter','latex','FontSize',11)


hold on
LFleg = plot(nan, nan, 'ko','DisplayName','LF','MarkerSize',6,'Linewidth',1.5);
HLLleg = plot(nan, nan, 'ks','DisplayName','HLL','MarkerSize',6,'Linewidth',1.5);
HLLCleg = plot(nan, nan, 'k^','DisplayName','HLLC','MarkerSize',6,'Linewidth',1.5);
Roeleg = plot(nan, nan, 'k+','DisplayName','Roe','MarkerSize',6,'Linewidth',1.5);
Trueleg = plot(nan, nan, 'r--','DisplayName','True','Linewidth',1.5);
lgd = legend([plot1,LFleg,HLLleg,HLLCleg,Roeleg],'interpreter','latex');
lgd.FontSize = 9;
% lgd.Layout.Tile = 'south';
% lgd.Orientation = 'horizontal';
% lgd.Layout.Tile = 'southeast';
% set(lgd.Position,[0.6520 0.1405 0.1387 0.1730])
ylabel('$u[\mathrm m/\mathrm s]$','interpreter','latex','FontSize',11)

hold off
% ylabel('$h(\mathrm m)$','interpreter','latex')

% grid minor

set(gcf, 'Position',  [300, 300, 600, 300])
set(gca,'LooseInset',get(gca,'TightInset'));  