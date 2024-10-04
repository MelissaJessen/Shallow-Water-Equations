
clear,clc,close all
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\utilities');
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\Linear\';

% Order of polymomials used for approximation 

NODETOL = 1e-10;
alpha = 0;
periodic = 0;

N = 6;
K = 5;
FinalTime = 20;


xmin = 0.0;
xmax = 2*pi;
A = 1;
g = 9.81;
T = 4;
omega = 2*pi/T;
k = 2*pi/(abs(xmax-xmin));
h0 = omega^2/(k^2*g);
etasol = @(x,t,A,omega,k) A*cos(omega*t-k*x); 
usol = @(x,t,A,omega,k,h0) omega/(k*h0)*A*cos(omega*t-k*x);
D = [sqrt(g*h0),0;0,-sqrt(g*h0)];

% Run initialization script
initialize_solver




eta0 = etasol(x,0,A,omega,k);
u0   = usol  (x,0,A,omega,k,h0);
%eta0 = 0*x;
%u0 = 0*x;
% Pack structs
struct_pack
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
%%
[u,eta,ut,etat,tt] = LinSWE1D(u0,eta0,x,FinalTime,1,param,maps);

%%
periodic
subplot(2,1,1)
a = plot(x,etasol(x,FinalTime,A,omega,k));

colors = cell2mat(get(a,'Color'));
hold on
for i = 1:size(x,2)
%     b = plot(x(:,i),eta(:,i),'ko','Color',colors(i,:),'MarkerSize',3);
    b = plot(x(:,i),eta(:,i),'kx');
    b.MarkerFaceColor = [0,0,0];
end
hold off
title('Free surface \eta')
grid('minor');
xlim([xmin,xmax])
ylim([1.1*min(min(eta)),1.1*max(max(eta))])

subplot(2,1,2)
plot(x,usol(x,FinalTime,A,omega,k,h0));
hold on
for i = 1:size(x,2)
%     b = plot(x(:,i),u(:,i),'o','Color',colors(i,:));
    b = plot(x(:,i),u(:,i),'kx');

    b.MarkerFaceColor = [0,0,0];
end
hold off
xlim([xmin,xmax])
ylim([1.1*min(min(u)),1.1*max(max(u))])

title('Velocity u')
grid('minor');



%% Time plots
close all
clear mov;
figure
printcount = 0;
filename  = [plotpath,'linswe.gif'];
figure
mov_count = 0;
umin = min(ut(:));
umax = max(ut(:));
hmin = min(etat(:));
hmax = max(etat(:));
for i = 1:length(tt)
    
    printcount = printcount + 1;
    if printcount == 20
        mov_count = mov_count + 1;
        title(['T = ',num2str(tt(i))])
        subplot(2,1,1)
        plot(x,squeeze(etat(i,:,:)))
        subtitle('Free surface elevation')
        xlim([xmin,xmax])
        ylim([hmin,hmax]*1.1)
        grid minor



        %     
    %     etasol = @(x,t,A,omega,k) A*cos(omega*t-k*x); 
    % usol = @(x,t,A,omega,k,h0) omega/(k*h0)*A*cos(omega*t-k*x);

        subplot(2,1,2)
%         subtitle('Velocity')
        plot(x,squeeze(ut(i,:,:)))
        subtitle('Velocity')

        xlim([xmin,xmax])
        ylim([umin,umax]*1.1)
        grid minor

        set(gcf,'color','w')

        drawnow
%         if mov_count ==1
%             gif(filename)
%         else
%             gif
%         end

        mov(mov_count) = getframe(gcf);
        printcount = 0;
        pause(0.01)
        
    end
end

% movie2gif(mov,[plotpath,'line_swe_reflection2-',date(),'.gif'], 'LoopCount', 0, 'DelayTime', 0)

% %% Time plots + true sol
% close all
% clear mov;
% printcount = 0;
% filename  = [plotpath,'linswe.gif'];
% figure
% mov_count = 0;
% for i = 1:length(tt)
%     
%     printcount = printcount + 1;
%     if printcount == 10
%         mov_count = mov_count + 1;
%         title(['T = ',num2str(tt(i))])
%         subplot(2,1,1)
% 
%         plot(x,etasol(x,tt(i),A,omega,k));
%         title('Free surface elevation')
% 
%         hold on
%         plot(x,squeeze(etat(i,:,:)),'kx')
%         hold off
%         xlim([xmin,xmax])
%         ylim([-1.1,1.1])
% 
% 
% 
%         %     
%     %     etasol = @(x,t,A,omega,k) A*cos(omega*t-k*x); 
%     % usol = @(x,t,A,omega,k,h0) omega/(k*h0)*A*cos(omega*t-k*x);
% 
%         subplot(2,1,2)
%         plot(x,usol(x,tt(i),A,omega,k,h0));
%         subtitle('Velocity')
%         hold on
%         plot(x,squeeze(ut(i,:,:)),'kx')
%         hold off
%         ylim([-7,7])
%         xlim([xmin,xmax])
% 
%         drawnow
% %         if mov_count ==1
% %             gif(filename)
% %         else
% %             gif
% %         end
% 
%         mov(mov_count) = getframe(gcf);
%         printcount = 0;
%         pause(0.01)
%         
%     end
% end
% % movie2gif(mov, filename, 'LoopCount', 0, 'DelayTime', 0)






