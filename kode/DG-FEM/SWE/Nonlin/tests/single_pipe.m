clear,clc,close all
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');


addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\analytical_solver');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\tests')
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin')
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\nonlin\single_pipe\';
% defplotname = 'single_pipe_fixed_height';
defplotname = 'single_pipe_sinflow';
datapath = 'C:\Users\Bruger\Documents\GitHub\Thesis\data\dgfem_data\';
% Should be run without slope limiting on h + b
Globals1D;

gif_toggle = 0;
NODETOL = 1e-10;
periodic = false;
BC = true;
% if BC
%     param.bc_type = 'dirichlet';
%     
% end
numframes = 200;
flux          = @HLLswe;
param.SL      = 1;
param.Limiter = @SlopeLim;
param.SLtype  = 5;
dt_factor     = 0.025;
N             = 2;
K             = 100;

A             = 1;
g             = 9.81;
slope = -0.005;

% FinalTime     = 10;
% xmin          = 0;
% xmax          = 25;
% hu_inflow =  @(x,t,h,hu,param,maps) 0.5;


FinalTime     = 1200;
xmin          = 0;
xmax          = 400;
hu_inflow =  @(x,t,h,hu,param,maps) 0.5*sin(4*pi/30*t)+1;
initialize_solver

b  = slope*x - xmax*slope;
bx = slope*ones(size(x));
h  = zeros(size(x)); hu_0 = 0.0; hu = ones(size(x))*hu_0; h0 = 0.1; ul = 0; ur = 0;
h_inflow = @(x,t,h,hu,param,maps) h(maps.vmapI);

h_outflow = @(x,t,h,hu,param,maps) h(maps.vmapO);hu_outflow =  @(x,t,h,hu,param,maps) hu(maps.vmapO);u_inflow = 0;u_outflow = 0;

h(:) = h0;
hrt = 0;hlt = 0;ult = 0;urt = 0;
u = hu./h;


%% Pack structs
struct_pack
param.g = g;
param.xmin   = xmin;
param.xmax   = xmax;
param.BC     = BC;
param.flux   = flux;
param.bx     = bx;
param.Minv   = Minv;
param.invV   = invV;
struct_pack
param.g = g;
param.xmin   = xmin;
param.xmax   = xmax;
param.BC     = BC;
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
% alpha = 0.0*pi;
mcells = length(x(:));
chalen = xmax;
timout = FinalTime;

figure
subplot(211)
hold on
plot(x(:),b(:),'k-','Linewidth',1.5)
plot(x(:),ones(size(x(:)))*h0+b(:),'k--')
ylabel('[$h(x,0)+b] (\mathrm m)$','interpreter','latex')
grid minor

hold off
subplot(212)
plot(x(:),hu(:),'k--','Linewidth',1)
ylabel('${hu(x,0)}$','interpreter','latex')
xlabel('$x$','interpreter','latex')
grid minor
% sgtitle(['T = ',num2str(FinalTime),'$s$'],'interpreter','latex')
set(gcf,'color','w')
set(gcf, 'Position',  [400, 400, 600, 250])
set(gca,'LooseInset',get(gca,'TightInset'));

% 


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

% %% Plots
% 
% Plot_hhu(x,h,hu,b)
% 
% saveas(gcf,PlotName(plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat','.eps'),'epsc');
% 
% saveas(gcf,PlotName(plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat','.png'),'png');
%% GIF
% close all
% clear mov;
% figure
% printcount = 0;
% filename  = [plotpath,defplotname,'T',num2str(FinalTime),'N',num2str(N),'.gif'];
% frameidx = 0;
% for i = 1:length(tt)
%     
%     printcount = printcount + 1;
% %     if (printcount == 10*FinalTime||i==length(tt))
%     if (printcount == floor(length(tt)/200)||i==length(tt))
% 
%         frameidx = frameidx + 1;
% %         title(['T = ',num2str(tt(i))])
% 
% 
%         subplot(211)
% %         plot(x,squeeze(ht(i,:,:)))
%         plot(x,squeeze(ht(i,:,:))+b)
%         hold on
%         plot(x,b,'k--')
%         hold off
%         title('Free surface height')
%         ylabel('$h$','interpreter','latex')
%         xlim([xmin,xmax])
% %         ylim([0,2.5])
% %         yticks(0:0.5:2.5)
%         grid minor
%         subplot(212)
%         plot(x,squeeze(hut(i,:,:)))
% %         hold on
% %         plot(x,b,'k--')
% %         hold off
%         ylabel('$hu$','interpreter','latex')
%         grid minor
% 
%         ylim([0,10])
%         xlim([xmin,xmax])
% %         ylim([0,30])
%         sgtitle(['T = ',num2str(tt(i)),'$s$'],'interpreter','latex');
%         set(gcf,'color','w')
% 
%         drawnow
%         mov(frameidx) = getframe(gcf);
%         printcount = 0;
% %         pause(0.1)
%     end
% end
% movie2gif(mov, filename, 'LoopCount', 3, 'DelayTime', 0)

%%
% 
if gif_toggle
close all
clear mov;
figure
printcount = 0;
filename  = [plotpath,defplotname,'T',num2str(FinalTime),'N',num2str(N),'.gif'];
frameidx = 0;
for i = 1:length(tt)
    
    printcount = printcount + 1;
%     if (printcount == 10*FinalTime||i==length(tt))
    if (printcount == floor(length(tt)/numframes)||i==length(tt))

        frameidx = frameidx + 1;
%         title(['T = ',num2str(tt(i))])


        subplot(211)
%         plot(x,squeeze(ht(i,:,:)))
        plot(x,squeeze(ht(:,:,i))+b)
        hold on
        plot(x,b,'k--')
        hold off
        title('Free surface height')
        ylabel('$h$','interpreter','latex')
        xlim([xmin,xmax])
        ylim([0,2.5])
        grid minor
        subplot(212)
        plot(x,squeeze(hut(:,:,i)))
%         hold on
%         plot(x,b,'k--')
%         hold off
        ylabel('$hu$','interpreter','latex')
        grid minor

        ylim([0,10])
        xlim([xmin,xmax])
%         ylim([0,30])
        sgtitle(['T = ',num2str(tt(i)),'$s$'],'interpreter','latex');
        set(gcf,'color','w')

        drawnow
        mov(frameidx) = getframe(gcf);
        printcount = 0;
%         pause(0.1)
    end
end
movie2gif(mov, filename, 'LoopCount', 3, 'DelayTime', 1)

end
%% spatial plots at different times (hu)
t_values = [100,400,700,950];

figure
t = tiledlayout(2,2,'TileSpacing','compact');
h1 = nexttile;
t1 = t_values(1);
[t1val,t1_idx] = min(abs(tt-t1));
plot(x,ht(:,:,t1_idx)+b, 'k-','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),b(:),'k--','Linewidth',1.5)
hold off
grid minor
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)
title(['t = ',num2str(tt(t1_idx))])
% xticklabels([])


h2 = nexttile;

t2 = t_values(2);
[t2val,t2_idx] = min(abs(tt-t2));
plot(x,ht(:,:,t2_idx)+b, 'k-','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),b(:),'k--','Linewidth',1.5)
hold off
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)
title(['t = ',num2str(tt(t2_idx))])
grid minor


% grid minor


h3 = nexttile;
t3 = t_values(3);
[t3val,t3_idx] = min(abs(tt-t3));
plot(x,ht(:,:,t3_idx)+b, 'k-','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),b(:),'k--','Linewidth',1.5)
hold off
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)
title(['t = ',num2str(tt(t3_idx))])
grid minor
h4 = nexttile;
t4 = t_values(4);
[t4val,t4_idx] = min(abs(tt-t4));
plot(x,ht(:,:,t4_idx)+b, 'k-','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),b(:),'k--','Linewidth',1.5)
hold off
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)
title(['t = ',num2str(tt(t4_idx))])
set(gcf, 'Position',  [300, 300, 600, 200])
grid minor
set(gca,'LooseInset',get(gca,'TightInset'));

%%t_values = [100,400,700,950];
%% HU
figure
t = tiledlayout(2,2,'TileSpacing','compact');
h1 = nexttile;
t1 = t_values(1);
[t1val,t1_idx] = min(abs(tt-t1));
plot(x,hut(:,:,t1_idx)+b, 'k-','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),b(:),'k--','Linewidth',1.5)
hold off
grid minor
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)
title(['t = ',num2str(tt(t1_idx))])
% xticklabels([])


h2 = nexttile;

t2 = t_values(2);
[t2val,t2_idx] = min(abs(tt-t2));
plot(x,hut(:,:,t2_idx)+b, 'k-','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),b(:),'k--','Linewidth',1.5)
hold off
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)
title(['t = ',num2str(tt(t2_idx))])
grid minor


% grid minor


h3 = nexttile;
t3 = t_values(3);
[t3val,t3_idx] = min(abs(tt-t3));
plot(x,hut(:,:,t3_idx)+b, 'k-','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),b(:),'k--','Linewidth',1.5)
hold off
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)
title(['t = ',num2str(tt(t3_idx))])
grid minor
h4 = nexttile;
t4 = t_values(4);
[t4val,t4_idx] = min(abs(tt-t4));
plot(x,hut(:,:,t4_idx)+b, 'k-','MarkerSize',4,'Linewidth',1.5)
hold on 
plot(x(:),b(:),'k--','Linewidth',1.5)
hold off
ylabel('$h+b[\mathrm m]$','interpreter','latex','FontSize',11)
title(['t = ',num2str(tt(t4_idx))])
set(gcf, 'Position',  [300, 300, 600, 200])
grid minor
set(gca,'LooseInset',get(gca,'TightInset'));
%% Fixed position varying time slots
figure
subplot(211)
xwhere = 322;
[x13,xidx] = min(abs(xx(:)-xwhere));

hold on
plot(tt,reshape(ht(2,81,:),[]),'k--')
ylabel('[$h(x,0)+b] (\mathrm m)$','interpreter','latex')
grid minor
hold off

subplot(212)
plot(x(:),hu(:),'k--','Linewidth',1)
ylabel('${hu(x,0)}$','interpreter','latex')
xlabel('$x$','interpreter','latex')
grid minor
% sgtitle(['T = ',num2str(FinalTime),'$s$'],'interpreter','latex')
set(gcf,'color','w')
set(gcf, 'Position',  [400, 400, 600, 200])
set(gca,'LooseInset',get(gca,'TightInset'));



%%
% figure
% hold on
% plot(tt,hu_inflow(x,tt,ht,hut,param,maps),'k','Linewidth',1.5)
% plot(tt,reshape(hut(1,1,:),[1,length(tt)]),'r--','Linewidth',1.5)
% hold off
% 
% xlim([0,FinalTime])
% 
% %%
% figure
% n0plot = 1;
% [XX,TT]= meshgrid(x(:),tt(n0plot:end));
% 
% hutplot = reshape(hut(n0plot:end,:,:),size(TT));
% htplot = reshape(ht(n0plot:end,:,:),size(TT));
% contourf(XX,TT,hutplot)
% colorbar
% 
% 
% %% Get element average
% h_modal = pagemtimes(invV,ht);
% h_modal(2:Np,:,:) = 0;
% h_avg = pagemtimes(V,h_modal);
% 
% h_avg = squeeze(h_avg(1,:,:));
% % h_avg = reshape(h_avg(1,:,:),[K,length(tt)]);
% 
% 
% hu_modal = pagemtimes(invV,hut);
% hu_modal(2:Np,:,:) = 0;
% hu_avg = pagemtimes(V,hu_modal);
% hu_avg = squeeze(hu_avg(1,:,:));

%%
% dat_filename = ['singlepipeCFL0025_T_',num2str(FinalTime),'_xmax_',num2str(xmax),'_',date];
% input_series = hu_inflow(x,tt,h,hu,param,maps);
% 
% % tab = table(x, tt,ht,hut,b,bx, 'VariableNames', { 'x', 't','h','hu','b','bx'} );
% filename = fullfile(datapath,dat_filename);
% % writetable(tab, [data_path,dat_filename])
% save(filename,'x','tt','ht','hut','h_avg','hu_avg','b','bx','input_series')
% % save [data_path,dat_filename] x tt ht hut b bx
