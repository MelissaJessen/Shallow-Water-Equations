clear,clc,close all
addpath('C:\Users\Matteo\Documents\Speciale\kode\DG-FEM\utilities');
addpath('C:\Users\Matteo\Documents\Speciale\kode\DG-FEM\SWE\Nonlin\analytical_solver');
addpath('C:\Users\Matteo\Documents\Speciale\kode\DG-FEM\SWE');
addpath('C:\Users\Matteo\Documents\Speciale\kode\DG-FEM\SWE\Nonlin')
plotpath = 'C:\Users\Matteo\Documents\Speciale';
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
param.Limiter = @SlopeLim;
flux = @LaxFswe;

param.SL = 1;
param.SLtype = 2;
dt_factor = 1;


N = 1;
K = 100;
FinalTime = 5.0;

xmin = 0;
xmax = 50;
A = 1;
g = 9.81;

gate = 25.0;
initialize_solver
hl = 0.1;
hr = 0.1;
hlt = hl;
hrt = hr;
ul = -3.0;
ur = 3.0;  
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
param.flux   = flux;
param.bx     = bx;
param.Minv   = Minv;
param.u_inflow = u_inflow;
param.Limiter = @SlopeLimit1;
param.b = b;
param.flux = @LaxFswe;



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
    [htrue,utrue] = drybed(h,u,analytical_param);
%     plot_title = ['Dry',' x_0=',num2str(gate),' T=',num2str(FinalTime),' h_L=',num2str(hl),' h_R=',num2str(hr),' u_L=',num2str(ul),' u_R=',num2str(ur), ' \alpha = ',num2str(alpha),', t = ',num2str(FinalTime)];

else
%     disp(['Wet',' BC: ',num2str(BC),' flux: ',flux])
    [htrue,utrue,converged_flag] = wetbed(h,u,niter,analytical_param);
%     plot_title = ['Wet',' x_0=',num2str(gate),' T=',num2str(FinalTime),' h_L=',num2str(hl),' h_R=',num2str(hr),' u_L=',num2str(ul),' u_R=',num2str(ur), ' \alpha = ',num2str(alpha),', t = ',num2str(FinalTime)];
end
%%
[h,hu,ht,hut,tt] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);

%% export htrue and utrue

% Define the file name
%filename = 'C:\Users\Matteo\Shallow-Water-Equations\data\torotest5.mat';

% Save the variables htrue and utrue into the .mat file
%save(filename, 'htrue', 'utrue', 'x');


%% plots
figure
subplot(2,1,1)
a = plot(x(:),h(:),'ko','MarkerSize',4);
ylabel('h','interpreter','latex')
hold on
plot(x(:),htrue(:),'r--','Linewidth',1.5)
legend('DG-FEM','True','location','southwest')
% plot(x,b,'k--')
grid minor

hold off
subplot(2,1,2)
plot(x(:),hu(:)./h(:),'ko','MarkerSize',4)
ylabel('u','interpreter','latex')
hold on
% plot(x,utrue+alpha*g*FinalTime,'k--')

plot(x(:),utrue(:),'r--','Linewidth',1.5)
hold off
grid minor
xlabel('x','interpreter','latex')
legend('DG-FEM','True','location','southwest')

% sgtitle(['T = ',num2str(FinalTime),'$s$'],'interpreter','latex')
set(gcf,'color','w')

% 
% saveas(gcf,[plotpath,defplotname,'_K_',num2str(K),'_flux_',flux,'.eps'],'epsc');
% 
% saveas(gcf,[plotpath,defplotname,'_K_',num2str(K),'_flux_',flux,'.png'],'png');

%% plot 2
figure
subplot(3,1,1)
a = plot(x,h,'ko','MarkerSize',4);
ylabel('h','interpreter','latex')
hold on
plot(x,htrue,'r--','Linewidth',1.5)

% plot(x,b,'k--')
grid minor

hold off
subplot(3,1,2)
plot(x,hu./h,'ko','MarkerSize',4)
ylabel('u','interpreter','latex')
hold on
% plot(x,utrue+alpha*g*FinalTime,'k--')

plot(x,utrue,'r--','Linewidth',1.5)
hold off
grid minor


subplot(3,1,3)
plot(x,hu,'ko','MarkerSize',4)
ylabel('hu','interpreter','latex')
hold on
% plot(x,utrue+alpha*g*FinalTime,'k--')

plot(x,htrue.*utrue,'r--','Linewidth',1.5)
hold off
grid minor

xlabel('x','interpreter','latex')
% sgtitle(['T = ',num2str(FinalTime),'$s$'],'interpreter','latex')
set(gcf,'color','w')

%% GIF
clear mov;
figure(2)
printcount = 0;
filename  = [plotpath,'nswe.gif'];
frameidx = 0;
for i = 1:length(tt)
    
    printcount = printcount + 1;
    if (printcount == 10||i==length(tt))
        frameidx = frameidx + 1;
%         title(['T = ',num2str(tt(i))])

        subplot(2,1,1)
        plot(x,squeeze(ht(i,:,:)))
        title('Free surface height')

        xlim([xmin,xmax])
%         ylim([-0.5,5])

        subplot(2,1,2)
        plot(x,squeeze(hut(i,:,:)./ht(i,:,:)))
        title('hu')

%         ylim([-3,8])
        xlim([xmin,xmax])
        drawnow
        mov(frameidx) = getframe(gcf);
        printcount = 0;
        pause(0.01)
    end
end
% movie2gif(mov, filename, 'LoopCount', 3, 'DelayTime', 0)


