clear,clc,close all
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\Linear';
% Order of polymomials used for approximation 

NODETOL = 1e-10;
FinalTime = 2;
periodic = 1;

xmin = 0.0;
xmax = 2*pi;
A = 1;
T = 2*pi;
g = 9.81;
omega = 2*pi/FinalTime;
k = 2*pi/(abs(xmax-xmin));
h0 = omega^2/(k^2*g);
etasol = @(x,t,A,omega,k) A*cos(omega*t-k*x); 
usol = @(x,t,A,omega,k,h0) omega/(k*h0)*A*cos(omega*t-k*x);
D = [sqrt(g*h0),0;0,-sqrt(g*h0)];
param.g      = g;
param.A      = A;
param.h0     = h0;
param.omega  = omega;
param.k      = k;
param.etasol = etasol;
param.usol   = usol;
param.D      = D;
param.xmin   = xmin;
param.xmax   = xmax;



% Polynomial konvergence
nk = [1,2,3,4,5,6];

Karray = 2.^nk;
Narray = [1,2,3];
% 
% k_line = 1./(Karray.^Narray);

nidx       = 0;
kidx       = 0;
errorKCper = zeros(length(Karray),length(Narray));
alpha = 1;

error_u_max     = zeros(length(Karray),length(Narray));
error_eta_max    = zeros(length(Karray),length(Narray));
error_u_l2     = zeros(length(Karray),length(Narray));
error_eta_l2    = zeros(length(Karray),length(Narray));
for K = Karray
    kidx = kidx + 1;
    nidx = 0;
        for N = Narray
        nidx = nidx + 1;
        initialize_solver
        struct_pack
        eta0 = etasol(x,0,A,omega,k);
        u0   = usol  (x,0,A,omega,k,h0);
        [u,eta] = LinSWE1D(u0,eta0,x,FinalTime,1,param,maps);
        utrue = usol(x,FinalTime,A,omega,k,h0);
        etatrue = etasol(x,FinalTime,A,omega,k);
        e_u = u-utrue;
        e_eta = eta-etatrue;

        err_u_loc = zeros(1,length(K));
        err_eta_loc = zeros(1,length(K));
        
            for i = 1:K
                 err_u_loc(i)    = e_u(:,i)'*Mk*e_u(:,i);
                 err_eta_loc(i)   = e_eta(:,i)'*Mk*e_eta(:,i);
            end
        error_u_max(kidx,nidx) = max(max(abs(e_u)));
        error_eta_max(kidx,nidx) = max(max(abs(e_eta)));
        error_u_l2(kidx,nidx) = sqrt(sum(err_u_loc));
        error_eta_l2(kidx,nidx) = sqrt(sum(err_eta_loc));

        end
%     
end

k_lines = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
    k_lines(:,i) = (1./Karray').^(Narray(i)+1);
end
%% plot

set(0,'defaultfigurecolor',[1 1 1])
set(groot,'DefaultAxesTickLabelInterpreter','Tex');

t = tiledlayout(2,2,'TileSpacing','tight');

% u max

k_lines_shift = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
k_lines_shift(:,i) = k_lines(:,i) * (error_u_max(1,i)/k_lines(1,i));
end

tilearray(1) = nexttile;
plot1 = loglog(Karray,k_lines_shift,'--','Linewidth',1.5);

hold on
for i = 1:length(Narray)
    loglog(Karray,error_u_max(:,i),'^','Color',plot1(i).Color,'Linewidth',1.5);
end

xticklabels([])

ylim([1e-8,1e1])
grid minor
ylabel('$||\textrm{error}||_{\infty}$','interpreter','latex','fontsize',14)
title('$u(\textrm{m/s})$','interpreter','latex')
% eta max
k_lines_shift = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
k_lines_shift(:,i) = k_lines(:,i) * (error_eta_max(1,i)/k_lines(1,i));
end

tilearray(2) = nexttile;
plot1 = loglog(Karray,k_lines_shift,'--','Linewidth',1.5);

hold on
for i = 1:length(Narray)
    loglog(Karray,error_eta_max(:,i),'^','Color',plot1(i).Color,'Linewidth',1.5);
end

ylim([1e-8,1e1])
grid minor
% ylabel('$||\textrm{error}_{h}||_{\infty}$','interpreter','latex','fontsize',14)
title('$h(\textrm{m})$','interpreter','latex')
xticklabels([])
yticklabels([])


% u l2
k_lines_shift = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
k_lines_shift(:,i) = k_lines(:,i) * (error_u_l2(1,i)/k_lines(1,i));
end

tilearray(3) = nexttile;
plot1 = loglog(Karray,k_lines_shift,'--','Linewidth',1.5);

hold on
for i = 1:length(Narray)
    loglog(Karray,error_u_l2(:,i),'^','Color',plot1(i).Color,'Linewidth',1.5);
end

ylim([1e-8,1e1])
grid minor
ylabel('$||\textrm{error}||_{L^2}$','interpreter','latex','fontsize',14)
xlabel('$K$','interpreter','latex')


k_lines_shift = zeros(length(Karray),length(Narray));
for i = 1:length(Narray)
k_lines_shift(:,i) = k_lines(:,i) * (error_eta_l2(1,i)/k_lines(1,i));
end

tilearray(4) = nexttile;
plot1 = loglog(Karray,error_eta_l2,'^','Linewidth',1.5);

hold on
for i = 1:length(Narray)
    leg = ['K','\^',num2str(-(Narray(i)+1))];
    plot2(i) = loglog(Karray,k_lines_shift(:,i),'--','Linewidth',1.5,'Color',plot1(i).Color,'Linewidth',1.5,'DisplayName',leg);

end
ylim([1e-8,1e1])
grid minor
yticklabels([])
xlabel('$K$','interpreter','latex')
legend(plot2)


set(gcf, 'Position',  [300, 300, 1100, 500])
set(gca,'LooseInset',get(gca,'TightInset'));

% saveas(gcf,[plotpath,defplotname],'epsc');
% saveas(gcf,[plotpath,defplotname],'png');

%%
% 
% a = loglog(Karray,errorKUper,'k--^');
% for i = 1:length(Narray)
% a(i).MarkerFaceColor = a(i).Color;
% klegends{i} = ['K^{(',num2str(-Narray(i)), ')}'];
% end
% 
% 
% % colors = cell2mat(get(a,'Color'))
% colors = {[1,0,0],[0,1,0],[0,0,1]};
% hold on
% for i = 1:length(Narray)
% % b = loglog(Karray,errorKCper(:,i),'-o','Color',colors(i,:));
% b = loglog(Karray,errorKCper(:,i),'k--o');
% 
% b.MarkerFaceColor = b.Color;
% % c.MarkerFaceColor = c.Color;
% end
% 
% for i = 1:length(Narray)
%     
%     c = loglog(Karray,k_lines(:,i),'-','Color',squeeze(colors{i}));
% end
% 
% hold off
% legend('Upwind','','','Central','','',klegends{1},klegends{2},klegends{3})
% grid('minor')
% xlabel('K')
% ylabel('E = ||q-q_{true}||_\infty')
% title('Periodic Domain')
% 
% 
% %% Non-periodic K convergence
% 
% periodic = 0;
% factor = 0.5;
% % Polynomial konvergence
% 
% nidx = 0;
% errorKC = zeros(1,length(Karray));
% alpha = 1;
% for K = Karray
%     nidx = nidx + 1;
%     initialize_solver
%     struct_pack
%     eta0 = etasol(x,0,A,omega,k);
%     u0   = usol  (x,0,A,omega,k,h0);
%     [u,eta] = LinSWE1D(u0,eta0,x,FinalTime,1,param,maps);
%     utrue = usol(x,FinalTime,A,omega,k,h0);
%     etatrue = etasol(x,FinalTime,A,omega,k);
%     q = [u,eta];
%     qtrue = [utrue,etatrue];
% 
%     errorKC(nidx) = max(max(abs(qtrue-q)));
% %     errorP(nidx) = max(max(abs(u-utrue)));
% %     
% end
% 
% errorKU = zeros(1,length(Karray));
% alpha = 0;  
% nidx = 0;
% for K = Karray
%     nidx = nidx + 1;
%     initialize_solver
%     struct_pack
%     eta0 = etasol(x,0,A,omega,k);
%     u0   = usol  (x,0,A,omega,k,h0);
%     [u,eta] = LinSWE1D(u0,eta0,x,FinalTime,1,param,maps);
%     utrue = usol(x,FinalTime,A,omega,k,h0);
%     etatrue = etasol(x,FinalTime,A,omega,k);
%     q = [u,eta];
%     qtrue = [utrue,etatrue];
% 
%     errorKU(nidx) = max(max(abs(qtrue-q)));
% %     errorP(nidx) = max(max(abs(u-utrue)));
% %     
% end
% %%
% figure
% loglog(Karray,errorKU,'k-*')
% hold on
% loglog(Karray,errorKC,'k-o')
% loglog(Karray,k_lines,'k--')
% hold off
% legend('Upwind','Central')
% grid('minor')
% xlabel('K')
% ylabel('E = ||q-q_{true}||_\infty')
% title('Non-periodic Domain')
% %%
% figure
% loglog(Karray,errorKU,'k-*')
% hold on
% loglog(Karray,errorKC,'k-o')
% loglog(Karray,errorKUper,'r-*')
% loglog(Karray,errorKCper,'r-o')
% loglog(Karray,k_line,'b--')
% hold off
% grid('minor')
% legend('Upwind','Central','Upwind Periodic','Central Periodic','1/K^N')
% 
% xlabel('K')
% ylabel('E = ||q-q_{true}||_\infty')
% 
% 
% % %%
% % 
% % errorKUper = zeros(length(Karray),length(Narray));
% % alpha = 0;
% % kidx  = 0;
% % nidx  = 0;
% % for K = Karray
% %     kidx = kidx + 1;
% %     nidx = 0;
% %         for N = Narray
% %         nidx = nidx + 1;
% %         initialize_solver
% %         struct_pack
% %         eta0 = etasol(x,0,A,omega,k);
% %         u0   = usol  (x,0,A,omega,k,h0);
% %         [u,eta] = LinSWE1D(u0,eta0,x,FinalTime,1,param,maps);
% %         utrue = usol(x,FinalTime,A,omega,k,h0);
% %         etatrue = etasol(x,FinalTime,A,omega,k);
% %         q = [u,eta];
% %         qtrue = [utrue,etatrue];
% %         errorKUper(kidx,nidx) = max(max(abs(qtrue-q)));
% %         end     
% % end
