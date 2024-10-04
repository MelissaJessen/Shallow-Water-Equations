clear,clc,close all
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');


addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\analytical_solver');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\tests')
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin')
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\nonlin\pipe\';
defplotname = 'pipe_case_presentation';

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
% if BC
%     param.bc_type = 'dirichlet';
%     
% end
flux = @LaxFswe;
param.SL = 1;
param.Limiter = @SlopeLim;
param.SLtype = 1;
dt_factor = 0.9;

N = 2;
K = 200;
FinalTime = 15;


xmin = 0;
xmax = 50;
A = 1;
g = 9.81;

gate = 2;

initialize_solver

x_intersection = [10,20,30,40];
alphas = -[0.02*pi,0.07*pi,0.01*pi,0.025*pi,0.01*pi];

[bx,b] = Build_slopes(x,x_intersection,alphas);
b = b-b(end);
figure
subplot(211)
plot(x,bx,'Linewidth',2)
ylabel('$b_x$','interpreter','latex')
subplot(212)
plot(x,b,'Linewidth',2)
ylabel('$b$','interpreter','latex')

set(gcf,'color','w')



h  = zeros(size(x));
% u = zeros(size(x));
hu = zeros(size(x));

hl = 0.33;
hr = 0.33;
ul = 0;
ur = 0;
% u_inflow = 0.18/hl;
h_inflow = @(x,t,h,hu,param,maps) h(maps.vmapI);
hu_inflow =  @(x,t,h,hu,param,maps) sin(10*t)+1.5;

h_outflow = @(x,t,h,hu,param,maps) h(maps.vmapO);
hu_outflow =  @(x,t,h,hu,param,maps) hu(maps.vmapO);
u_inflow = 0;
u_outflow = 0;
hu(x<=gate) = 0;
hu(x>gate) = 0;

h(x<=gate) = hl;
h(x>gate) = hr;

hrt = 0;
hlt = 0;
ult = 0;
urt = 0;
% u(x<=gate) = ul;
% u(x>gate) = ur;



u = hu./h;
% u = u.*flow_factor;
% inbetween = [zeros(size(x(:))),flipud(b(:))];
% bottomcolor = [1,1,1]*0.7;
% x2 = [x(:),flipud(x(:))];
% 


figure
subplot(311)
hold on
plot(x,h+b)
plot(x,b,'k--')
hold off
grid minor
ylabel('$h+b$','interpreter','latex')
subplot(312)
hold on
plot(x,u)
plot(x,b,'k--')

% fill(x2, inbetween, bottomcolor);

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
sgtitle('Initial conditions','interpreter','latex')
set(gcf,'color','w')
set(gcf,'color','w')

saveas(gcf,[plotpath,defplotname,'IC'],'epsc')

saveas(gcf,[plotpath,defplotname,'IC'],'png')



% figure
% subplot(211)
% hold on
% plot(x,b,'Linewidth',2)
% plot(x,h+b,'k--')
% hold off
% ylim([0,0.6])
% subplot(212)
% plot(x,bx,'Linewidth',2)
% % 



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
param.invV   = invV;
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

pack_analytical_parameters
dcrit = (ur - ul) - 2.0*(cl + cr);


% 
% if(hl<=0 || hr<=0 || dcrit>=0)
%     disp(['Dry',' BC: ',num2str(BC),' flux: ',flux])
%     [htrue,utrue] = drybed(h,u,analytical_param);
%     plot_title = ['Dry',' x_0=',num2str(gate),' T=',num2str(FinalTime),' h_L=',num2str(hl),' h_R=',num2str(hr),' u_L=',num2str(ul),' u_R=',num2str(ur), ' \alpha = ',num2str(alpha),', t = ',num2str(FinalTime)];
% 
% else
%     disp(['Wet',' BC: ',num2str(BC),' flux: ',flux])
%     [htrue,utrue,converged_flag] = wetbed(h,u,niter,analytical_param);
%     plot_title = ['Wet',' x_0=',num2str(gate),' T=',num2str(FinalTime),' h_L=',num2str(hl),' h_R=',num2str(hr),' u_L=',num2str(ul),' u_R=',num2str(ur), ' \alpha = ',num2str(alpha),', t = ',num2str(FinalTime)];
% end
%%
[h,hu,ht,hut,tt] = NSWE1D(h,hu,x,FinalTime,dt_factor,@NSWE_RHS1D,param,maps);



%% Plots
% 
Plot_hhu(x,h,hu,b)
% 
% saveas(gcf,PlotName(plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat','.eps'),'epsc');
% 
% saveas(gcf,PlotName(plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat','.png'),'png');

%%
% SWE2gif(x,ht,hut,b,tt,plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat',100,false,false)

%% GIF
close all
clear mov;
figure
printcount = 0;
filename  = [plotpath,defplotname,'.gif'];
frameidx = 0;
for i = 1:length(tt)
    
    printcount = printcount + 1;
%     if (printcount == 10*FinalTime||i==length(tt))
    if (printcount == floor(length(tt)/50)||i==length(tt))

        frameidx = frameidx + 1;
%         title(['T = ',num2str(tt(i))])


        subplot(211)
        plot(x,squeeze(ht(:,:,i))+b)
        hold on
        plot(x,b,'k--')
        hold off
%         title('Free surface height')
        ylabel('$h$','interpreter','latex')
        xlim([xmin,xmax])
%         ylim([0,2.5])
%         yticks(0:0.5:2.5)
        grid minor
%         ylim([-0.5,5])
        
%         subplot(312)
%         plot(x,squeeze(hut(i,:,:)./ht(i,:,:)))
%         hold on
%         plot(x,b,'k--')
%         hold off
%         ylabel('$u$','interpreter','latex')
%         grid minor
%         ylim([0,15])

        subplot(212)
        plot(x,squeeze(hut(:,:,i)))
        hold on
        plot(x,b,'k--')
        hold off
        ylabel('$hu$','interpreter','latex')
        grid minor

%         ylim([-3,8])
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
movie2gif(mov, filename, 'LoopCount', 1, 'DelayTime', 0)

% 
% %% Gif with true solution
% 
% 
% clear mov;
% figure(2)
% printcount = 0;
% filename  = [plotpath,'nswe_andtrue_bumb.gif'];
% frameidx = 0;
% for i = 1:length(tt)
%     
%     printcount = printcount + 1;
%     if (printcount == 10||i==length(tt))
%         frameidx = frameidx + 1;
% %         title(['T = ',num2str(tt(i))])
% 
%         analytical_param.timout = tt(i);
%         [htrue,utrue,converged_flag] = wetbed(h,u,niter,analytical_param);
% 
%         subplot(2,1,1)
% %         hold on
%         plot(x,squeeze(ht(i,:,:)))
% %         plot(x,htrue,'k--');
% %         title('Free surface height')
%         hold off
% 
%         xlim([xmin,xmax])
%         ylim([0,inf])
% 
%         subplot(2,1,2)
% %         hold on
%         plot(x,squeeze(hut(i,:,:)./ht(i,:,:)))
% %         plot(x,utrue,'k--')
% %         hold off
%         
%         title('u')
% 
% %         ylim([-3,8])
% %         xlim([xmin,xmax])
%         drawnow
%         mov(frameidx) = getframe(gcf);
%         printcount = 0;
% %         pause(0.01)
%     end
% end
% % movie2gif(mov, filename, 'LoopCount', 3, 'DelayTime', 0)
