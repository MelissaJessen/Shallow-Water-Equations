clear,clc,close all
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\analytical_solver');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\tests')
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin')
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\nonlin\LaiKhan\';
defplotname = 'LaiKhanTest2';
Globals1D;
NODETOL = 1e-10;
periodic = false;
BC = true;
flux = @HLLswe;
param.SL = true;
dt_factor = 1;
param.SLtype = 3;
param.Limiter = @SlopeLim;


N = 2;
K = 100;
FinalTime = 2;



xmin = 0;
xmax = 25;
A = 1;
g = 9.81;

gate = 20;




initialize_solver

% parabolic bump;
alpha = 0.05;
b = 0.2 - alpha*(x-10).^2;
a8 = find(x==8);
a12 = find(x==12);
b(x<8) = 0;
b(12<x) = 0;




% title('bottom geometry')
% % 
% alpha = 0.05;
bx = -2*alpha*(x-10);
bx(x<8) = 0;
bx(12<x) = 0;
bx(a8(1)) = 0;
bx(a12(2)) = 0;

% bx = sgolayfilt(bx,3,3);
hl = 0.5;
hr = 0.5;

hu_inflow_cst = 0.18;

h_inflow = @(x,t,h,hu,param,maps) h(maps.vmapI);
hu_inflow =  @(x,t,h,hu,param,maps) 0.18;

h_outflow = @(x,t,h,hu,param,maps) h(maps.vmapO);
% h_outflow = @(x,t,h,hu,param,maps) 0.5;

hu_outflow =  @(x,t,h,hu,param,maps) hu(maps.vmapO);
ul = hu_inflow_cst/hl;
ur = hu_inflow_cst/hl;
u_inflow = 0.18/hl;


u_outflow = 0;
h  = zeros(size(x));
hu = zeros(size(x));
u = zeros(size(x));

initialize_solver
% param.u_inflow = u_inflow;
param.u_outflow = 0;
% param.u_inflow = u_inflow;
param.h_inflow = h_inflow;
param.hu_inflow = hu_inflow;
param.h_outflow = h_outflow;
param.hu_outflow = hu_outflow;





hu(:) = hu_inflow_cst;
h(:) = hl;
% u(:) = hu_inflow_cst/hl;
h = h-b;
% hu = h.*u;

u = hu./h;
% u = u.*flow_factor;


figure
subplot(311)
hold on
plot(x,h+b,'Linewidth',2)
plot(x,b,'k--','Linewidth',2)
hold off
grid minor
ylabel('$h+b$','interpreter','latex')
yticks(0:0.1:0.6)
subplot(312)
hold on
plot(x,u,'Linewidth',2)
% plot(x,b,'k--')
hold off
ylabel('u','interpreter','latex')
grid minor

subplot(313)
hold on
plot(x,hu,'Linewidth',2)
% plot(x,b,'k--')
hold off
ylabel('hu','interpreter','latex')
grid minor
sgtitle('Initial conditions','interpreter','latex')
set(gcf,'color','w')
% 
% saveas(gcf,[plotpath,defplotname,'IC'],'epsc')
% 
% saveas(gcf,[plotpath,defplotname,'IC'],'png')

set(gcf, 'Position',  [300, 300, 1100, 500])

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
param.x      = x;
param.b      = b;


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



% Plot_hhu(x,h,hu,b)
% 
% saveas(gcf,PlotName(plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat','.eps'),'epsc');
% 
% saveas(gcf,PlotName(plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat','.png'),'png');
% 
% 
% 
% %% GIF
% clear mov;
% figure(2)
% printcount = 0;
% filename  = [plotpath,'nswe.gif'];
% frameidx = 0;
% for i = 1:length(tt)
%     
%     printcount = printcount + 1;
%     if (printcount == 10||i==length(tt))
%         frameidx = frameidx + 1;
% %         title(['T = ',num2str(tt(i))])
% 
%         subplot(2,1,1)
%         plot(x,squeeze(ht(i,:,:))+b)
%         title('Free surface height')
% 
%         xlim([xmin,xmax])
% %         ylim([-0.5,5])
% 
%         subplot(2,1,2)
%         plot(x,squeeze(hut(i,:,:)./ht(i,:,:)))
%         title('hu')
% 
% %         ylim([-3,8])
%         xlim([xmin,xmax])
%         drawnow
%         mov(frameidx) = getframe(gcf);
%         printcount = 0;
%         pause(0.01)
%     end
% end
% %%
% SWE2gif(x,ht,hut,b,tt,plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat',100,true,true);
% 
% %%

%% plots
figure
defplotname = 'LaiKhanTest2_h_and_b_';

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
defplotname = 'LaiKhanTest2_only_h_';

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

dat_filename = ['LK2_T_',num2str(FinalTime)];
% tab = table(x, tt,ht,hut,b,bx, 'VariableNames', { 'x', 't','h','hu','b','bx'} );
filename = fullfile(datapath,dat_filename);

save(filename,'x','h','hu','b','bx')
