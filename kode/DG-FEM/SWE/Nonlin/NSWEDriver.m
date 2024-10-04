clear,clc,close all
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\analytical_solver');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE');
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\nonlin\random_test\';
defplotname = 'BC_test';

Globals1D;

NODETOL = 1e-10;
periodic = false;
BC = true;
flux = 'LF';
param.SL = false;
N = 2;
K = 100;
FinalTime = 10;


xmin = 0;
xmax = 50;
A = 1;
g = 9.81;

gate = 20;
initialize_solver
hl = 3.5;
hr = 1.25;

hu_inflow = 2;
% ul = hu_inflow/hl;
ul = 0;

ur = 0;
h_inflow = hl;
u_inflow = hu_inflow/h_inflow;
h  = zeros(size(x));
u = zeros(size(x));

h(x<=gate) = hl;
h(x>gate) = hr;


u(x<=gate) = ul;
u(x>gate) = ur;


% No slope
bx = 0;
alpha = 0;
b = 0;

% alpha =0;
% x_intersection = [25,35];
% alphas = [0.02*pi,0.3*pi,0.05*pi];
% 
% [bx,b] = Build_slopes(x,x_intersection,alphas);

% % Linear sloped bottom
% alpha = 0.005*pi;% bottom slope
% % alpha = 0
% center = 0;
% bx = ones(size(x))*alpha;
% b = (x-center)*alpha;
% b = b - b(end);
% 

% bx(25>=x) = ones(size(x))*alpha2;
% 

% figure
% subplot(211)
% plot(x,bx)
% subplot(212)
% plot(x,b)

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
param.invV   = invV;
param.u_inflow = u_inflow;
param.hu_inflow = hu_inflow;
param.h_inflow = h_inflow;


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



if(hl<=0 || hr<=0 || dcrit>=0)
    disp(['Dry',' BC: ',num2str(BC),' flux: ',flux])
    [htrue,utrue] = drybed(h,u,analytical_param);
    plot_title = ['Dry',' x_0=',num2str(gate),' T=',num2str(FinalTime),' h_L=',num2str(hl),' h_R=',num2str(hr),' u_L=',num2str(ul),' u_R=',num2str(ur), ' \alpha = ',num2str(alpha),', t = ',num2str(FinalTime)];

else
    disp(['Wet',' BC: ',num2str(BC),' flux: ',flux])
    [htrue,utrue,converged_flag] = wetbed(h,u,niter,analytical_param);
    plot_title = ['Wet',' x_0=',num2str(gate),' T=',num2str(FinalTime),' h_L=',num2str(hl),' h_R=',num2str(hr),' u_L=',num2str(ul),' u_R=',num2str(ur), ' \alpha = ',num2str(alpha),', t = ',num2str(FinalTime)];
end
 %%
[h,hu,ht,hut,tt] = NSWE1D(h,hu,x,FinalTime,0.1,param,maps);

hu(1)

%%
figure
subplot(2,1,1)
plot(x,h)
ylabel('h','interpreter','latex')
hold on
plot(x,htrue,'k--')
% plot(x,b,'k--')
hold off
subplot(2,1,2)
plot(x,hu)
ylabel('hu','interpreter','latex')
hold on
% plot(x,utrue+alpha*g*FinalTime,'k--')
plot(x,htrue.*utrue,'k--')
hold off

saveas(gcf,PlotName(plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat','.epsc'));

saveas(gcf,PlotName(plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat','.png'));


%% GIF
clear mov;
figure(2)
printcount = 0;
filename  = [plotpath,defplotname,'T',num2str(FinalTime),'N',num2str(N),flux,'SL',num2str(param.SL),'.gif'];
frameidx = 0;
for i = 1:length(tt)
    
    printcount = printcount + 1;
    if (printcount == floor(length(tt)/25)||i==length(tt))

        frameidx = frameidx + 1;
%         title(['T = ',num2str(tt(i))])


        subplot(211)
        plot(x,squeeze(ht(i,:,:))+b)
        hold on
%         plot(x,b,'k--')
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
        plot(x,squeeze(hut(i,:,:)))
        hold on
%         plot(x,b,'k--')
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
        pause(0.01)
    end
end
movie2gif(mov, PlotName(plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat','.gif'), 'LoopCount', 3, 'DelayTime', 0)


%% Gif with true solution


clear mov;
figure(2)
printcount = 0;
filename  = [plotpath,'nswe_andtrue_bumb.gif'];
frameidx = 0;
for i = 1:length(tt)
    
    printcount = printcount + 1;
    if (printcount == 10||i==length(tt))
        frameidx = frameidx + 1;
%         title(['T = ',num2str(tt(i))])

        analytical_param.timout = tt(i);
        [htrue,utrue,converged_flag] = wetbed(h,u,niter,analytical_param);

        subplot(2,1,1)
%         hold on
        plot(x,squeeze(ht(i,:,:)))
%         plot(x,htrue,'k--');
%         title('Free surface height')
        hold off

        xlim([xmin,xmax])
%         ylim([-0.5,5])

        subplot(2,1,2)
%         hold on
        plot(x,squeeze(hut(i,:,:)./ht(i,:,:)))
%         plot(x,utrue,'k--')
%         hold off
        
        title('u')

%         ylim([-3,8])
%         xlim([xmin,xmax])
        drawnow
        mov(frameidx) = getframe(gcf);
        printcount = 0;
        pause(0.01)
    end
end
% movie2gif(mov, filename, 'LoopCount', 3, 'DelayTime', 0)
