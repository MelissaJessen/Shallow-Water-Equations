clear,clc,close all
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\analytical_solver');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\tests')
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin')
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\nonlin\LaiKhan\';
defplotname = 'LaiKhanTest1';

% Order of polymomials used for approximation 
NODETOL = 1e-10;

periodic = false;
BC = 1;
flux = @Roe_swe;
param.SL = 0;
dt_factor = 0.1;
param.SLtype = 2;
param.Limiter = @SlopeLim;
N = 2;
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


%%
[h,hu,ht,hut,tt] = NSWE1D(h,hu,x,FinalTime,dt_factor,@NSWE_RHS1D_MMS,param,maps);


% h = h+b;


%%
subplot(311)
htrue  = hfun(x,FinalTime);
hutrue = hufun(x,FinalTime);

hold on
plot(x,h,'k','Linewidth',2)
plot(x,htrue,'r--','Linewidth',2)
hold off
subplot(312)
hold on
plot(x,hu,'k','Linewidth',2)
plot(x,hutrue,'r--','Linewidth',2)
hold off

subplot(313)
hold on
plot5 = plot(x,abs(h-htrue),'kx-','MarkerSize',3,'DisplayName','h');
plot6 = plot(x,abs(hu-hutrue),'bo-','MarkerSize',3,'DisplayName','hu');
hold off
% saveas(gcf,PlotName(plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat','.eps'),'epsc');
% 
% saveas(gcf,PlotName(plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat','.png'),'png');


% 1/N*(h(end,:)-htrue(end,:))*M*(h(end,:)-htrue(end,:))'
%%

% (h(:,end)-htrue(:,end))'*M*(h(:,end)-htrue(:,end))
%%
SWE2gif(x,ht,hut,b,tt,plotpath,defplotname,FinalTime,N,flux,param.SL,['D','O'],'flat',100,[],[],false);

%% GIF
clear mov;
figure(2)
printcount = 0;
filename  = [plotpath,'LaiKhantest2.gif'];
frameidx = 0;
for i = 1:length(tt)
    
    printcount = printcount + 1;
    if (printcount == 10||i==length(tt))
        frameidx = frameidx + 1;
%         title(['T = ',num2str(tt(i))])

        subplot(2,1,1)
        plot(x,squeeze(ht(i,:,:)))
        hold on
        plot(x,b,'k--')
        hold off
        title('Free surface height')

        xlim([xmin,xmax])
%         ylim([-0.5,5])

        subplot(2,1,2)
        plot(x,squeeze(hut(i,:,:)))
        title('hu')

%         ylim([-3,8])
        xlim([xmin,xmax])
        set(gcf,'color','w')

        drawnow
        mov(frameidx) = getframe(gcf);
        printcount = 0;
        pause(0.01)
    end
end
% movie2gif(mov, filename, 'LoopCount', 3, 'DelayTime', 0)


%% max and L2 norms
% nk = [1,2,3,4,5,6,7,8,9];


% nk = [1,2,3,4,5,6,7,8,9];
% Narray = [1,2,3,4];

nk = [1,2,3,4,5,6,7];
Narray = [1,2,3];
Karray = 2.^nk;
% 
% k_line = 1./(Karray.^Narray);
tic
nidx       = 0;
kidx       = 0;
error_h_max     = zeros(length(Karray),length(Narray));
error_hu_max    = zeros(length(Karray),length(Narray));
error_h_L2  = zeros(length(Karray),length(Narray));
error_hu_L2 = zeros(length(Karray),length(Narray));

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
        [h,hu] = NSWE1D(h0,hu0,x,FinalTime,dt_factor,@NSWE_RHS1D_MMS,param,maps);

        eh = h-htrue;
        ehu = hu-hutrue;
        errh = zeros(1,K);
        errhu = zeros(1,K);
        for i = 1:K
        errh(i)  = (h(:,i)-htrue(:,i))'*Mk*(h(:,i)-htrue(:,i));
        errhu(i) = (hu(:,i)-hutrue(:,i))'*Mk*(hu(:,i)-hutrue(:,i));
        end
        error_h_L2(kidx,nidx)  = sqrt(sum(errh));
        error_hu_L2(kidx,nidx) = sqrt(sum(errhu));
        error_h_max(kidx,nidx)     = max(max(abs(eh)));
        error_hu_max(kidx,nidx)    = max(max(abs(ehu)));

%         k_lines = [k_lines,1./Karray'.^(N+1)];         
        end
end
toc
%%
k_lines = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
    k_lines(:,i) = (1./Karray').^(Narray(i)+1);
end

% 1/N*(h(end,:)-htrue(end,:))*M*(h(end,:)-htrue(end,:))'
%% Max norm h


plot1 = loglog(Karray,error_h_max,'o','MarkerSize',4,'Linewidth',1.5);

k_lines_shift = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
k_lines_shift(:,i) = k_lines(:,i) * (error_h_max(1,i)/k_lines(1,i));
end
hold on
%  b = plot(Karray(i),k_lines(i,:),'--');
% figure
for i = 1:length(plot1)
    leg = ['K','\^',num2str(-(Narray(i)+1))];
    b(i) = plot(Karray,k_lines_shift(:,i),'--','Color',plot1(i).Color,'DisplayName',leg,'Linewidth',1.5);

end

legend(b)

hold off
xlabel('K','interpreter','latex')
ylabel('$||\mathbf{h}-\mathbf{h}_{true}||_{\infty}$','interpreter','latex')
% for l = 1:length(a)
%     set(b(l),{'Color'},num2cell(a(l).Color,2))
% end

grid minor

%% max norm HU

plot1 = loglog(Karray,error_hu_max,'o','MarkerSize',4,'Linewidth',1.5);

% k_lines_shift = zeros(size(k_lines));
k_lines_shift = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
k_lines_shift(:,i) = k_lines(:,i) * (error_hu_max(1,i)/k_lines(1,i));
end
hold on
%  b = plot(Karray(i),k_lines(i,:),'--');
% figure
for i = 1:length(plot1)
    leg = ['K','\^',num2str(-(Narray(i)+1))];
    b(i) = plot(Karray,k_lines_shift(:,i),'--','Color',plot1(i).Color,'DisplayName',leg,'Linewidth',1.5);

end

% legend(b)

hold off
xlabel('K','interpreter','latex')
ylabel('$||\mathbf{hu}-\mathbf{hu}_{true}||_{\infty}$','interpreter','latex')
% for l = 1:length(a)
%     set(b(l),{'Color'},num2cell(a(l).Color,2))
% end

grid minor

%% L2 Norm h

plot1 = loglog(Karray,error_h_L2,'o','MarkerSize',4,'Linewidth',1.5);
k_lines_shift = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
k_lines_shift(:,i) = k_lines(:,i) * (error_h_L2(1,i)/k_lines(1,i));
end
hold on
%  b = plot(Karray(i),k_lines(i,:),'--');
% figure
for i = 1:length(plot1)
    leg = ['K','\^',num2str(-(Narray(i)+1))];
    b(i) = plot(Karray,k_lines_shift(:,i),'--','Color',plot1(i).Color,'DisplayName',leg,'Linewidth',1.5);

end

% legend(b)

hold off
xlabel('K','interpreter','latex')
ylabel('$||\mathbf{U}-\mathbf{U}_{true}||_{L^2}$','interpreter','latex')
% for l = 1:length(a)
%     set(b(l),{'Color'},num2cell(a(l).Color,2))
% end

drawnow
grid minor


