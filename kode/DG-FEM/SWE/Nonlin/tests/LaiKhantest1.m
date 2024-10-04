clear,clc,close all
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\analytical_solver');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\tests')
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin')
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\nonlin\LaiKhan\';

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
flux = @HLLswe;
% flux = @LaxFswe;

param.SL = 1;
dt_factor = 0.5;
param.SLtype = 5;
param.Limiter = @SlopeLim;

N = 2;
K = 100;
FinalTime = 1;

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

% b(a8(1)) = 0;
% b(a12(2)) = 0;

h = h-b;% Det virker med at trække b fra


bx = -2*alpha*(x-10);
bx(x<8) = 0;
bx(12<x) = 0;
bx(a8(1)) = 0;
bx(a12(2)) = 0;


figure
subplot(211)
hold on
plot(x,b,'Linewidth',2)
% plot(x(:),ones(size(x(:)))*hr,'k--')
ylabel('Bottom height ($\mathrm m$)','interpreter','latex')
grid minor

hold off
subplot(212)
plot(x,bx,'Linewidth',2)
ylabel('Bottom slope','interpreter','latex')
xlabel('$x$','interpreter','latex')
grid minor
% sgtitle(['T = ',num2str(FinalTime),'$s$'],'interpreter','latex')
set(gcf,'color','w')
set(gcf, 'Position',  [400, 400, 600, 300])
set(gca,'LooseInset',get(gca,'TightInset'));

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
% param.b      = b;
param.Minv   = Minv;
param.u_inflow = u_inflow;
param.h_inflow = h_inflow;
param.hu_inflow = hu_inflow;
param.h_outflow = h_outflow;
param.hu_outflow = hu_outflow;
param.u_outflow = 0;
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


%%
[h,hu,ht,hut,tt] = NSWE1D(h,hu,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);


% h = h+b;


%%
% Plot_hhu(x,h,hu,b)
% 
% saveas(gcf,PlotName(plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat<','.eps'),'epsc');
% 
% saveas(gcf,PlotName(plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat','.png'),'png');

%%
% SWE2gif(x,ht,hut,b,tt,plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat',100,[],[],false);

%% plots
figure
defplotname = 'LaiKhanTest1_h_and_b_';

subplot(121)
hold on
plot(x,h+b,'ko','MarkerSize',3,'Linewidth',1)
hold off
grid minor
title('$h+b$ [m]','interpreter','latex')
xlim([xmin,xmax])
% yticks(0:0.5:2.5)
xlabel('$x$ [m]','interpreter','latex')

subplot(1,2,2)
hold on
plot(x,hu,'ko','MarkerSize',3,'Linewidth',1)
% plot(x,b,'k--')
hold off
title('$hu$ [m$^2$s$^{-1}$]','interpreter','latex')
grid minor
xlim([xmin,xmax])
xlabel('$x$ [m]','interpreter','latex')

sgtitle(['T = ',num2str(FinalTime),'$s$'],'interpreter','latex')
set(gcf,'color','w')

set(gcf, 'Position',  [300, 300, 600, 300])

% saveas(gcf,[plotpath,defplotname,'T_',num2str(FinalTime)],'epsc')
% 
% saveas(gcf,[plotpath,defplotname,'T',num2str(FinalTime)],'png')

%%
figure
defplotname = 'LaiKhanTest1_only_h_';

subplot(121)
hold on
plot(x,h,'ko','MarkerSize',3,'Linewidth',1)
hold off
grid minor
title('$h$ [m]','interpreter','latex')
xlim([xmin,xmax])
% yticks(0:0.5:2.5)
xlabel('$x$ [m]','interpreter','latex')

subplot(1,2,2)
hold on
plot(x,hu,'ko','MarkerSize',3,'Linewidth',1)
% plot(x,b,'k--')
hold off
title('$hu$ [m$^2$s$^{-1}$]','interpreter','latex')
grid minor
xlim([xmin,xmax])
xlabel('$x$ [m]','interpreter','latex')

sgtitle(['T = ',num2str(FinalTime),'$s$'],'interpreter','latex')
set(gcf,'color','w')

set(gcf, 'Position',  [300, 300, 600, 300])

% saveas(gcf,[plotpath,defplotname,'T_',num2str(FinalTime)],'epsc')
% 
% saveas(gcf,[plotpath,defplotname,'T',num2str(FinalTime)],'png')


%%

% %% GIF
% clear mov;
% figure(2)
% printcount = 0;
% filename  = [plotpath,'LaiKhantest2.gif'];
% frameidx = 0;
% for i = 1:length(tt)
%     
%     printcount = printcount + 1;
%     if (printcount == 10||i==length(tt))
%         frameidx = frameidx + 1;
% %         title(['T = ',num2str(tt(i))])
% 
%         subplot(2,1,1)
%         plot(x,squeeze(ht(i,:,:)))
%         hold on
%         plot(x,b,'k--')
%         hold off
%         title('Free surface height')
% 
%         xlim([xmin,xmax])
% %         ylim([-0.5,5])
% 
%         subplot(2,1,2)
%         plot(x,squeeze(hut(i,:,:)))
%         title('hu')
% 
% %         ylim([-3,8])
%         xlim([xmin,xmax])
%         set(gcf,'color','w')
% 
%         drawnow
%         mov(frameidx) = getframe(gcf);
%         printcount = 0;
%         pause(0.01)
%     end
% end
% movie2gif(mov, filename, 'LoopCount', 3, 'DelayTime', 0)

datapath = 'C:\Users\Bruger\Documents\GitHub\Thesis\data\dgfem_data\';

dat_filename = ['LK1',num2str(FinalTime)];
% tab = table(x, tt,ht,hut,b,bx, 'VariableNames', { 'x', 't','h','hu','b','bx'} );
filename = fullfile(datapath,dat_filename);

save(filename,'x','h','hu','b','bx')
