clear,clc,close all
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\analytical_solver');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\tests')
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin')
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\nonlin\MMS\';
set(0,'defaultfigurecolor',[1 1 1])
set(groot,'DefaultAxesTickLabelInterpreter','Tex');
% Order of polymomials used for approximation 
NODETOL = 1e-10;

periodic = false;
BC = 1;
flux = @Roe_swe;
param.SL = 1;
dt_factor = 0.1;
param.SLtype = 3;
% param.Limiter = @SlopeLim;
param.Limiter = @SlopeLimitN;





N = 10;
K = 200;
FinalTime = 0.1;


xmin = 0;
xmax = 1;
A = 1;
g = 9.81;

gate = 20;


initialize_solver
hl = 0.33;
hr = 0.33;
ul = 0;
ur = 0;
% u_inflow = 0.18/hl;
h_inflow = hl;
hu_inflow = 0;
u_inflow = 0;
u_outflow = 0;
% u_inflow = 0;
hfun  = @(x,t) exp(x/(xmax-xmin))*(cos(t) + 3);
hufun = @(x,t) exp(x/(xmax-xmin))*(sin(t) + 3);
h  = hfun(x,0);
hu = hufun(x,0);
u = hu./h;

% parabolic bump;
alpha = 0.05;
b = 0.2 - alpha*(x-10).^2;
% a8 = find(x==8);
% a12 = find(x==12);

[~,a8] = min(abs(x(:)-8));

[~,a12] = min(abs(x(:)-12));

b = 0*b;
h = h-b;% Det virker med at trække b fra


bx = -2*alpha*(x-10);


bx = 0*bx;


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
param.u_inflow  = u_inflow;
param.h_inflow  = @(x,t,h,hu,param,varargin) hfun(x,t);
param.hu_inflow = @(x,t,h,hu,param,varargin) hufun(x,t);

param.h_outflow  = @(x,t,h,hu,param,varargin) hfun(x,t);
param.hu_outflow = @(x,t,h,hu,param,varargin) hufun(x,t);
param.u_outflow = 0;
param.Limiter = @SlopeLim;
param.SLtype = 2;
param.b = b;
param.L = xmax-xmin;
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


%% max and L2 norms
% nk = [1,2,3,4,5,6,7,8,9];


% nk = [1,2,3,4,5,6,7,8,9];
% Narray = [1,2,3,4];

nk = [1,2,3,4,5,6,7];
Narray = [1,2,3,4];
Karray = 2.^nk;
% 
% k_line = 1./(Karray.^Narray);
tic
nidx       = 0;
kidx       = 0;
error_h_max_lf     = zeros(length(Karray),length(Narray));
error_hu_max_lf    = zeros(length(Karray),length(Narray));
error_h_l2_lf  = zeros(length(Karray),length(Narray));
error_hu_l2_lf = zeros(length(Karray),length(Narray));

error_h_max_hll     = zeros(length(Karray),length(Narray));
error_hu_max_hll    = zeros(length(Karray),length(Narray));
error_h_l2_hll  = zeros(length(Karray),length(Narray));
error_hu_l2_hll = zeros(length(Karray),length(Narray));

error_h_max_hllc  = zeros(length(Karray),length(Narray));
error_hu_max_hllc = zeros(length(Karray),length(Narray));
error_h_l2_hllc   = zeros(length(Karray),length(Narray));
error_hu_l2_hllc  = zeros(length(Karray),length(Narray));

error_h_max_roe  = zeros(length(Karray),length(Narray));
error_hu_max_roe = zeros(length(Karray),length(Narray));
error_h_l2_roe   = zeros(length(Karray),length(Narray));
error_hu_l2_roe  = zeros(length(Karray),length(Narray));



% k_lines    = [];
alpha = 1;
for K = Karray
    kidx = kidx + 1;
    nidx = 0;
        for N = Narray
        nidx = nidx + 1;
        initialize_solver
        struct_pack; param.x = x;param.b = zeros(size(x));param.bx = zeros(size(x));
        h0  = hfun(x,0);
        hu0 = hufun(x,0);
        htrue  = hfun(x,FinalTime);
        hutrue = hufun(x,FinalTime);
        param.flux = @LaxFswe;
        [h_lf,hu_lf] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D_MMS,param,maps);
        param.flux = @HLLswe;
        [h_hll,hu_hll] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D_MMS,param,maps);
        param.flux = @HLLCswe;
        [h_hllc,hu_hllc] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D_MMS,param,maps);
        param.flux = @Roe_swe;
        [h_roe,hu_roe] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D_MMS,param,maps);

        eh_lf    = h_lf-htrue;
        ehu_lf   = hu_lf-hutrue;
        eh_hll   = h_hll-htrue;
        ehu_hll  = hu_hll-hutrue;
        eh_hllc  = h_hllc-htrue;
        ehu_hllc = hu_hllc-hutrue;
        eh_roe   = h_roe-htrue;
        ehu_roe  = hu_roe-hutrue;
        
        
        errh_lf    = zeros(1,K);
        errhu_lf   = zeros(1,K);
        errh_hll   = zeros(1,K);
        errhu_hll  = zeros(1,K);
        errh_hllc  = zeros(1,K);
        errhu_hllc = zeros(1,K);
        errh_roe   = zeros(1,K);
        errhu_roe  = zeros(1,K);
        for i = 1:K
        errh_lf(i)    = (h_lf(:,i)-htrue(:,i))'*Mk*(h_lf(:,i)-htrue(:,i));
        errhu_lf(i)   = (hu_lf(:,i)-hutrue(:,i))'*Mk*(hu_lf(:,i)-hutrue(:,i));
        
        errh_hll(i)   = (h_hll(:,i)-htrue(:,i))'*Mk*(h_hll(:,i)-htrue(:,i));
        errhu_hll(i)  = (hu_hll(:,i)-hutrue(:,i))'*Mk*(hu_hll(:,i)-hutrue(:,i));
        
        errh_hllc(i)  = (h_hllc(:,i)-htrue(:,i))'*Mk*(h_hllc(:,i)-htrue(:,i));
        errhu_hllc(i) = (hu_hllc(:,i)-hutrue(:,i))'*Mk*(hu_hllc(:,i)-hutrue(:,i));
        
        errh_roe(i)  = (h_roe(:,i)-htrue(:,i))'*Mk*(h_roe(:,i)-htrue(:,i));
        errhu_roe(i) = (hu_roe(:,i)-hutrue(:,i))'*Mk*(hu_roe(:,i)-hutrue(:,i));
        
        end
        error_h_l2_lf(kidx,nidx)   = sqrt(sum(errh_lf));
        error_hu_l2_lf(kidx,nidx)  = sqrt(sum(errhu_lf));
        error_h_max_lf(kidx,nidx)  = max(max(abs(eh_lf)));
        error_hu_max_lf(kidx,nidx) = max(max(abs(ehu_lf)));
        
        error_h_l2_hll(kidx,nidx)   = sqrt(sum(errh_hll));
        error_hu_l2_hll(kidx,nidx)  = sqrt(sum(errhu_hll));
        error_h_max_hll(kidx,nidx)  = max(max(abs(eh_hll)));
        error_hu_max_hll(kidx,nidx) = max(max(abs(ehu_hll)));  
        
        error_h_l2_hllc(kidx,nidx)   = sqrt(sum(errh_hllc));
        error_hu_l2_hllc(kidx,nidx)  = sqrt(sum(errhu_hllc));
        error_h_max_hllc(kidx,nidx)  = max(max(abs(eh_hllc)));
        error_hu_max_hllc(kidx,nidx) = max(max(abs(ehu_hllc)));
        
        error_h_l2_roe(kidx,nidx)   = sqrt(sum(errh_roe));
        error_hu_l2_roe(kidx,nidx)  = sqrt(sum(errhu_roe));
        error_h_max_roe(kidx,nidx)  = max(max(abs(eh_roe)));
        error_hu_max_roe(kidx,nidx) = max(max(abs(ehu_roe)));     
        end
end
toc

k_lines = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
    k_lines(:,i) = (1./Karray').^(Narray(i)+1);
end

%%
% subplot(311)
% htrue  = hfun(x,FinalTime);
% hutrue = hufun(x,FinalTime);
% 
% hold on
% plot(x,h_roe,'k','Linewidth',2)
% plot(x,htrue,'r--','Linewidth',2)
% hold off
% subplot(312)
% hold on
% plot(x,hu_roe,'k','Linewidth',2)
% plot(x,hutrue,'r--','Linewidth',2)
% hold off
% 
% subplot(313)
% hold on
% plot5 = plot(x,abs(h_roe-htrue),'kx-','MarkerSize',3,'DisplayName','h');
% plot6 = plot(x,abs(hu_roe-hutrue),'bo-','MarkerSize',3,'DisplayName','hu');
% hold off


%% subplots
% Max norm h
figure
% h1 = subplot(2,2,1);
t = tiledlayout(2,2,'TileSpacing','tight');
h1 = nexttile;
% h1 = subplot(221)
defplotname = 'slopelimited_convergence';
plot1 = loglog(Karray,error_h_max_lf,'o','MarkerSize',4,'Linewidth',1.5);

k_lines_shift = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
k_lines_shift(:,i) = k_lines(:,i) * (error_h_max_lf(1,i)/k_lines(1,i));
end
hold on
%  b = plot(Karray(i),k_lines(i,:),'--');
% figure
for i = 1:length(plot1)
    leg = ['K','\^',num2str(-(Narray(i)+1))];
    plot2(i) = plot(Karray,k_lines_shift(:,i),'--','Color',plot1(i).Color,'DisplayName',leg,'Linewidth',1.5);
    plot3(i) = loglog(Karray,error_h_max_hll(:,i),'S','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
    plot4(i) = loglog(Karray,error_h_max_hllc(:,i),'^','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
    plot5(i) = loglog(Karray,error_h_max_roe(:,i),'+','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
end
% LFleg = plot(nan, nan, 'ko','DisplayName','LF','MarkerSize',4,'Linewidth',1.5);
% HLLleg = plot(nan, nan, 'ks','DisplayName','HLL','MarkerSize',4,'Linewidth',1.5);
% HLLCleg = plot(nan, nan, 'k^','DisplayName','HLLC','MarkerSize',4,'Linewidth',1.5);
% Roeleg = plot(nan, nan, 'k+','DisplayName','Roe','MarkerSize',4,'Linewidth',1.5);



% legend([LFleg,HLLleg,HLLCleg,Roeleg,plot2],'Location','sw')
% xlim(get(gca,'Xlim').*[1,1.1])
% ylim(get(gca,'Ylim').*[1,2])
xlim([1e0,0.33*1e3])
xticklabels([])
ylim([1e-15,1e2])
% pos = [1e2,5e-15, 1.5e2, 1e-11];

% abc = gca;
% pos = [plot5(4).XData(end)/1.2,plot5(4).YData(end)/20, plot5(4).XData(end)/2, plot5(4).YData(end)*20];
% rectangle('Position',pos);


hold off
% xlabel('$\mathrm K$','interpreter','latex')
% ylabel('$||\mathbf{h}-\mathbf{h}_{\mathrm {true}}||_{\infty}$','interpreter','latex')
ylabel('$||\textrm{error}||_{\infty}$','interpreter','latex','fontsize',14)
title('Free surface ($h$)','interpreter','latex','fontsize',14)
grid minor
h2 = nexttile;
% h2 = subplot(222);

% Max norm hu
% defplotname = 'max_hu';

plot1 = loglog(Karray,error_hu_max_lf,'o','MarkerSize',4,'Linewidth',1.5);

k_lines_shift = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
k_lines_shift(:,i) = k_lines(:,i) * (error_hu_max_lf(1,i)/k_lines(1,i));
end
hold on
%  b = plot(Karray(i),k_lines(i,:),'--');
% figure
for i = 1:length(plot1)
    leg = ['K','\^',num2str(-(Narray(i)+1))];
    plot2(i) = plot(Karray,k_lines_shift(:,i),'--','Color',plot1(i).Color,'DisplayName',leg,'Linewidth',1.5);
    plot3(i) = loglog(Karray,error_hu_max_hll(:,i),'S','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
    plot4(i) = loglog(Karray,error_hu_max_hllc(:,i),'^','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
    plot5(i) = loglog(Karray,error_hu_max_roe(:,i),'+','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
end
% LFleg = plot(nan, nan, 'ko','DisplayName','LF','MarkerSize',4,'Linewidth',1.5);
% HLLleg = plot(nan, nan, 'ks','DisplayName','HLL','MarkerSize',4,'Linewidth',1.5);
% HLLCleg = plot(nan, nan, 'k^','DisplayName','HLLC','MarkerSize',4,'Linewidth',1.5);
% Roeleg = plot(nan, nan, 'k+','DisplayName','Roe','MarkerSize',4,'Linewidth',1.5);

title('Flow rate ($hu$)','interpreter','latex','fontsize',14)


% legend([LFleg,HLLleg,HLLCleg,Roeleg,plot2],'Location','sw')
% xlim(get(gca,'Xlim').*[1,1.1])
% ylim(get(gca,'Ylim').*[1,2])
xlim([1e0,0.33*1e3])
ylim([1e-15,1e2])
xticklabels([])
% yticklabels([])
% pos = [plot5(4).XData(end)/1.2,plot5(4).YData(end)/20, plot5(4).XData(end)/2, plot5(4).YData(end)*20];
% rectangle('Position',pos);
hold off
% xlabel('$\mathrm K$','interpreter','latex')
% ylabel('$||\mathbf{hu}-\mathbf{hu}_{\mathrm {true}}||_{\infty}$','interpreter','latex')

grid minor
h3 = nexttile;
% h3 = subplot(223);

% L2 Norm h

plot1 = loglog(Karray,error_h_l2_lf,'o','MarkerSize',4,'Linewidth',1.5);

k_lines_shift = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
k_lines_shift(:,i) = k_lines(:,i) * (error_h_l2_lf(1,i)/k_lines(1,i));
end
hold on
%  b = plot(Karray(i),k_lines(i,:),'--');
% figure
for i = 1:length(plot1)
    leg = ['K','\^',num2str(-(Narray(i)+1))];
    plot2(i) = plot(Karray,k_lines_shift(:,i),'--','Color',plot1(i).Color,'DisplayName',leg,'Linewidth',1.5);
    plot3(i) = loglog(Karray,error_h_l2_hll(:,i),'S','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
    plot4(i) = loglog(Karray,error_h_l2_hllc(:,i),'^','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
    plot5(i) = loglog(Karray,error_h_l2_roe(:,i),'+','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
end
% LFleg = plot(nan, nan, 'ko','DisplayName','LF','MarkerSize',4,'Linewidth',1.5);
% HLLleg = plot(nan, nan, 'ks','DisplayName','HLL','MarkerSize',4,'Linewidth',1.5);
% HLLCleg = plot(nan, nan, 'k^','DisplayName','HLLC','MarkerSize',4,'Linewidth',1.5);
% Roeleg = plot(nan, nan, 'k+','DisplayName','Roe','MarkerSize',4,'Linewidth',1.5);



% legend([LFleg,HLLleg,HLLCleg,Roeleg,plot2],'Location','sw')
% xlim(get(gca,'Xlim').*[1,1.1])
% ylim(get(gca,'Ylim').*[1,2])

xlim([1e0,0.33*1e3])
ylim([1e-15,1e2])
% pos = [plot5(4).XData(end)/1.2,plot5(4).YData(end)/20, plot5(4).XData(end)/2, plot5(4).YData(end)*20];
% rectangle('Position',pos);
hold off
xlabel('Elements ($\mathrm K$)','interpreter','latex')
% ylabel('$||\mathbf{h}-\mathbf{h}_{\mathrm {true}}||_{\mathrm L^2}$','interpreter','latex')
ylabel('$||\textrm{error}||_{\mathrm L^2}$','interpreter','latex','fontsize',14)


grid minor
h4 = nexttile;

% h4 = subplot(224)

% L2 norm hu
% h4 = subplot(224);

plot1 = loglog(Karray,error_hu_l2_lf,'o','MarkerSize',4,'Linewidth',1.5);

k_lines_shift = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
k_lines_shift(:,i) = k_lines(:,i) * (error_hu_l2_lf(1,i)/k_lines(1,i));
end
hold on
%  b = plot(Karray(i),k_lines(i,:),'--');
% figure
for i = 1:length(plot1)
    leg = ['K','\^',num2str(-(Narray(i)+1))];
    plot2(i) = plot(Karray,k_lines_shift(:,i),'--','Color',plot1(i).Color,'DisplayName',leg,'Linewidth',1.5);
    plot3(i) = loglog(Karray,error_hu_l2_hll(:,i),'S','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
    plot4(i) = loglog(Karray,error_hu_l2_hllc(:,i),'^','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
    plot5(i) = loglog(Karray,error_hu_l2_roe(:,i),'+','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
end
LFleg = plot(nan, nan, 'ko','DisplayName','LF','MarkerSize',4,'Linewidth',1.5);
HLLleg = plot(nan, nan, 'ks','DisplayName','HLL','MarkerSize',4,'Linewidth',1.5);
HLLCleg = plot(nan, nan, 'k^','DisplayName','HLLC','MarkerSize',4,'Linewidth',1.5);
Roeleg = plot(nan, nan, 'k+','DisplayName','Roe','MarkerSize',4,'Linewidth',1.5);

xlim([1e0,0.33*1e3])
ylim([1e-15,1e2])
% yticklabels([])

% xlim(get(gca,'Xlim').*[1,2])
% ylim(get(gca,'Ylim').*[2,3])
pos = [plot5(4).XData(end)/1.2,plot5(4).YData(end)/20, plot5(4).XData(end)/2, plot5(4).YData(end)*20];

hold off
xlabel('Elements ($\mathrm K$)','interpreter','latex')
% ylabel('$||\mathbf{hu}-\mathbf{hu}_{\mathrm {true}}||_{\mathrm L^2}$','interpreter','latex')
lgd = legend([LFleg,HLLleg,HLLCleg,Roeleg,plot2]);
lgd.Layout.Tile = 'west';
grid minor

% rectangle('Position',pos);
% axes('position',[0.225 0.15 0.05 0.1])
% box on
% loglog(Karray(end-2:end),k_lines_shift(end-2:end,4),'--','Color',plot1(i).Color,'Linewidth',1.5)
% hold on
% set(gca,'xtick',[]);set(gca,'xticklabel',[]);set(gca,'ytick',[]);set(gca,'yticklabel',[])
% xlim([plot5(4).XData(end)/1.2,plot5(4).XData(end)/2+plot5(4).XData(end)/1.2])
% ylim([plot5(4).YData(end)/20,plot5(4).YData(end)/20 + plot5(4).YData(end)*20])
% loglog(Karray(end),error_hu_l2_hll(end,end),'S','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
% loglog(Karray(end),error_hu_l2_hllc(end,end),'^','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
% loglog(Karray(end),error_hu_l2_roe(end,end),'+','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
% 
% hold off
set(gcf, 'Position',  [300, 300, 1100, 500])
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,[plotpath,defplotname],'epsc');
saveas(gcf,[plotpath,defplotname],'png');

% 
% 
% %% Max norm h
% figure
% defplotname = 'max_h';
% plot1 = loglog(Karray,error_h_max_lf,'o','MarkerSize',4,'Linewidth',1.5);
% 
% k_lines_shift = zeros(length(Karray),length(Narray));
% for i = 1:length(Narray)
% k_lines_shift(:,i) = k_lines(:,i) * (error_h_max_lf(1,i)/k_lines(1,i));
% end
% hold on
% %  b = plot(Karray(i),k_lines(i,:),'--');
% % figure
% for i = 1:length(plot1)
%     leg = ['K','\^',num2str(-(Narray(i)+1))];
%     plot2(i) = plot(Karray,k_lines_shift(:,i),'--','Color',plot1(i).Color,'DisplayName',leg,'Linewidth',1.5);
%     plot3(i) = loglog(Karray,error_h_max_hll(:,i),'S','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
%     plot4(i) = loglog(Karray,error_h_max_hllc(:,i),'^','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
%     plot5(i) = loglog(Karray,error_h_max_roe(:,i),'+','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
% end
% LFleg = plot(nan, nan, 'ko','DisplayName','LF','MarkerSize',4,'Linewidth',1.5);
% HLLleg = plot(nan, nan, 'ks','DisplayName','HLL','MarkerSize',4,'Linewidth',1.5);
% HLLCleg = plot(nan, nan, 'k^','DisplayName','HLLC','MarkerSize',4,'Linewidth',1.5);
% Roeleg = plot(nan, nan, 'k+','DisplayName','Roe','MarkerSize',4,'Linewidth',1.5);
% 
% 
% 
% legend([LFleg,HLLleg,HLLCleg,Roeleg,plot2],'Location','sw')
% xlim(get(gca,'Xlim').*[1,1.1])
% ylim(get(gca,'Ylim').*[1,2])
% 
% hold off
% xlabel('$\mathrm K$','interpreter','latex')
% ylabel('$||\mathbf{h}-\mathbf{h}_{\mathrm {true}}||_{\infty}$','interpreter','latex')
% 
% grid minor
% 
% saveas(gcf,[plotpath,defplotname],'epsc');
% saveas(gcf,[plotpath,defplotname],'png');
% 
% 
% %% Max norm hu
% figure
% defplotname = 'max_hu';
% 
% plot1 = loglog(Karray,error_hu_max_lf,'o','MarkerSize',4,'Linewidth',1.5);
% 
% k_lines_shift = zeros(length(Karray),length(Narray));
% for i = 1:length(Narray)
% k_lines_shift(:,i) = k_lines(:,i) * (error_hu_max_lf(1,i)/k_lines(1,i));
% end
% hold on
% %  b = plot(Karray(i),k_lines(i,:),'--');
% % figure
% for i = 1:length(plot1)
%     leg = ['K','\^',num2str(-(Narray(i)+1))];
%     plot2(i) = plot(Karray,k_lines_shift(:,i),'--','Color',plot1(i).Color,'DisplayName',leg,'Linewidth',1.5);
%     plot3(i) = loglog(Karray,error_hu_max_hll(:,i),'S','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
%     plot4(i) = loglog(Karray,error_hu_max_hllc(:,i),'^','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
%     plot5(i) = loglog(Karray,error_hu_max_roe(:,i),'+','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
% end
% LFleg = plot(nan, nan, 'ko','DisplayName','LF','MarkerSize',4,'Linewidth',1.5);
% HLLleg = plot(nan, nan, 'ks','DisplayName','HLL','MarkerSize',4,'Linewidth',1.5);
% HLLCleg = plot(nan, nan, 'k^','DisplayName','HLLC','MarkerSize',4,'Linewidth',1.5);
% Roeleg = plot(nan, nan, 'k+','DisplayName','Roe','MarkerSize',4,'Linewidth',1.5);
% 
% 
% 
% legend([LFleg,HLLleg,HLLCleg,Roeleg,plot2],'Location','sw')
% xlim(get(gca,'Xlim').*[1,1.1])
% ylim(get(gca,'Ylim').*[1,2])
% 
% 
% hold off
% xlabel('$\mathrm K$','interpreter','latex')
% ylabel('$||\mathbf{hu}-\mathbf{hu}_{\mathrm {true}}||_{\infty}$','interpreter','latex')
% 
% grid minor
% saveas(gcf,[plotpath,defplotname],'epsc');
% saveas(gcf,[plotpath,defplotname],'png');
% %% L2 Norm h
% figure
% defplotname = 'L2_h';
% 
% plot1 = loglog(Karray,error_h_l2_lf,'o','MarkerSize',4,'Linewidth',1.5);
% 
% k_lines_shift = zeros(length(Karray),length(Narray));
% for i = 1:length(Narray)
% k_lines_shift(:,i) = k_lines(:,i) * (error_h_l2_lf(1,i)/k_lines(1,i));
% end
% hold on
% %  b = plot(Karray(i),k_lines(i,:),'--');
% % figure
% for i = 1:length(plot1)
%     leg = ['K','\^',num2str(-(Narray(i)+1))];
%     plot2(i) = plot(Karray,k_lines_shift(:,i),'--','Color',plot1(i).Color,'DisplayName',leg,'Linewidth',1.5);
%     plot3(i) = loglog(Karray,error_h_l2_hll(:,i),'S','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
%     plot4(i) = loglog(Karray,error_h_l2_hllc(:,i),'^','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
%     plot5(i) = loglog(Karray,error_h_l2_roe(:,i),'+','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
% end
% LFleg = plot(nan, nan, 'ko','DisplayName','LF','MarkerSize',4,'Linewidth',1.5);
% HLLleg = plot(nan, nan, 'ks','DisplayName','HLL','MarkerSize',4,'Linewidth',1.5);
% HLLCleg = plot(nan, nan, 'k^','DisplayName','HLLC','MarkerSize',4,'Linewidth',1.5);
% Roeleg = plot(nan, nan, 'k+','DisplayName','Roe','MarkerSize',4,'Linewidth',1.5);
% 
% 
% 
% legend([LFleg,HLLleg,HLLCleg,Roeleg,plot2],'Location','sw')
% xlim(get(gca,'Xlim').*[1,1.1])
% ylim(get(gca,'Ylim').*[1,2])
% 
% hold off
% xlabel('$\mathrm K$','interpreter','latex')
% ylabel('$||\mathbf{h}-\mathbf{h}_{\mathrm {true}}||_{\mathrm L^2}$','interpreter','latex')
% 
% grid minor
% saveas(gcf,[plotpath,defplotname],'epsc');
% saveas(gcf,[plotpath,defplotname],'png');
% %% L2 norm hu
% figure
% % subplot(2,2,1)
% defplotname = 'L2_hu';
% 
% plot1 = loglog(Karray,error_hu_l2_lf,'o','MarkerSize',4,'Linewidth',1.5);
% 
% k_lines_shift = zeros(length(Karray),length(Narray));
% for i = 1:length(Narray)
% k_lines_shift(:,i) = k_lines(:,i) * (error_hu_l2_lf(1,i)/k_lines(1,i));
% end
% hold on
% %  b = plot(Karray(i),k_lines(i,:),'--');
% % figure
% for i = 1:length(plot1)
%     leg = ['K','\^',num2str(-(Narray(i)+1))];
%     plot2(i) = plot(Karray,k_lines_shift(:,i),'--','Color',plot1(i).Color,'DisplayName',leg,'Linewidth',1.5);
%     plot3(i) = loglog(Karray,error_hu_l2_hll(:,i),'S','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
%     plot4(i) = loglog(Karray,error_hu_l2_hllc(:,i),'^','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
%     plot5(i) = loglog(Karray,error_hu_l2_roe(:,i),'+','Color',plot1(i).Color,'MarkerSize',4,'Linewidth',1.5);
% end
% LFleg = plot(nan, nan, 'ko','DisplayName','LF','MarkerSize',4,'Linewidth',1.5);
% HLLleg = plot(nan, nan, 'ks','DisplayName','HLL','MarkerSize',4,'Linewidth',1.5);
% HLLCleg = plot(nan, nan, 'k^','DisplayName','HLLC','MarkerSize',4,'Linewidth',1.5);
% Roeleg = plot(nan, nan, 'k+','DisplayName','Roe','MarkerSize',4,'Linewidth',1.5);
% 
% 
% 
% legend([LFleg,HLLleg,HLLCleg,Roeleg,plot2],'Location','sw')
% xlim(get(gca,'Xlim').*[1,1.1])
% ylim(get(gca,'Ylim').*[1,2])
% 
% 
% hold off
% xlabel('$\mathrm K$','interpreter','latex')
% ylabel('$||\mathbf{hu}-\mathbf{hu}_{\mathrm {true}}||_{\mathrm L^2}$','interpreter','latex')
% 
% grid minor
% saveas(gcf,[plotpath,defplotname],'epsc');
% saveas(gcf,[plotpath,defplotname],'png');