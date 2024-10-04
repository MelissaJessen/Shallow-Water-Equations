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
flux = @HLLCswe;
param.SL = true;
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

% parabolic bump;
alpha = 0.05;
b = 0.2 - alpha*(x-10).^2;
a8 = find(x==8);
a12 = find(x==12);
b(x<=8) = 0;
b(12<=x) = 0;
% b
% b(a8(1)) = 0;
% b(a12(2)) = 0;
% b
% plot(x,b)



% title('bottom geometry')
% % 
% alpha = 0.05;
bx = -2*alpha*(x-10);
bx(x<8) = 0;
bx(12<x) = 0;
% bx(a8(1)) = 0;
% bx(a12(2)) = 0;
hl = 2;
hr = 2;

hu_inflow_cst = 25.0567;

h_inflow = @(x,t,h,hu,param,maps) h(maps.vmapI);
hu_inflow =  @(x,t,h,hu,param,maps) hu_inflow_cst;

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
% u = u.*flow_factor;

% 
% figure
% subplot(311)
% hold on
% plot(x,h+b)
% plot(x,b,'k--')
% hold off
% grid minor
% ylabel('$h+b$','interpreter','latex')
% yticks(0:0.1:0.6)
% subplot(312)
% hold on
% plot(x,u)
% plot(x,b,'k--')
% hold off
% ylabel('u','interpreter','latex')
% grid minor
% 
% subplot(313)
% hold on
% plot(x,hu)
% plot(x,b,'k--')
% hold off
% ylabel('hu','interpreter','latex')
% grid minor
% sgtitle('Initial conditions','interpreter','latex')
% set(gcf,'color','w')
% 
% saveas(gcf,[plotpath,defplotname,'IC'],'epsc')
% 
% saveas(gcf,[plotpath,defplotname,'IC'],'png')


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
[h,hu,ht,hut,tt] = NSWE1D(h,hu,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);


%% plots
figure

subplot(311)
hold on
plot(x,h+b)
plot(x,b,'k--')
hold off
grid minor
ylabel('$h+b$','interpreter','latex')
% ylim([0,2.5])
yticks(0:0.5:2.5)
subplot(312)
hold on
plot(x,hu./h)
plot(x,b,'k--')
hold off
ylabel('u','interpreter','latex')
grid minor

subplot(313)
hold on
plot(x,hu)
plot(x,b,'k--')
hold off
ylabel('hu','interpreter','latex')
grid minor
sgtitle(['T = ',num2str(FinalTime),'$s$'],'interpreter','latex')
set(gcf,'color','w')

saveas(gcf,[plotpath,defplotname,'T',num2str(FinalTime)],'epsc')

saveas(gcf,[plotpath,defplotname,'T',num2str(FinalTime)],'png')



%% GIF
% close all
clear mov;
figure
printcount = 0;
filename  = [plotpath,defplotname,'.gif'];
frameidx = 0;
for i = 1:length(tt)
    
    printcount = printcount + 1;
    if (printcount == floor(length(tt)/25)||i==length(tt))
        frameidx = frameidx + 1;
%         title(['T = ',num2str(tt(i))])


        subplot(211)
        plot(x,squeeze(ht(i,:,:))+b)
        hold on
        plot(x,b,'k--')
        hold off
%         title('Free surface height')
        ylabel('$h$','interpreter','latex')
        xlim([xmin,xmax])
%         ylim([0,2.5])
%         yticks(0:0.5:2.5)
        grid minor
        subplot(212)
        plot(x,squeeze(hut(i,:,:)))
%         ylim([0.16,0.19])
%         hold on
% %         plot(x,b,'k--')
%         hold off
        ylabel('$hu$','interpreter','latex')
        grid minor

        ylim([20,30])
        xlim([xmin,xmax])
%         ylim([0,30])
        sgtitle(['T = ',num2str(tt(i)),'$s$'],'interpreter','latex');
        set(gcf,'color','w')

        drawnow
        mov(frameidx) = getframe(gcf);
        printcount = 0;
%         pause(0.01)
        end
end
% movie2gif(mov, filename, 'LoopCount', 3, 'DelayTime', 0)



