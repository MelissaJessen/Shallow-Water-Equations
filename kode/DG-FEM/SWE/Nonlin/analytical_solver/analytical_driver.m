clear,clc,close all

addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\utilities');
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\Analytical\';
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
NODETOL = 1e-10;
alpha = 0;
periodic = 0;

N = 1;
K = 100;


xmin = 0.0;
xmax = 50;
chalen = xmax-xmin; % Channel length

A = 1;
g = 9.81;
initialize_solver
x = linspace(xmin,xmax,400);


toggle_gif = 0;

% hl   = 1;  % left height
% ul   = -5;
% hr   = 1; % right height
% ur   = 5;
% gate = 25;

hl = 1;
hr = 0;
ul = 0;
ur = 0;
gate = 25;
FinalTime = 2.5;

% d in fortran is h in my code
h          = zeros(size(x(:)));
h(x(:)<=gate) = hl;
gate_idx = sum(x(:)<=gate);
h(x(:)>gate)  = hr;
u          = zeros(size(x(:)));
u(x(:)<=gate) = ul;
u(x(:)>gate)  = ur;
hu         = u.*h;
timout     = FinalTime;
dt = FinalTime/100;
% dt = 0.1;

% TIME = linspace(0,FinalTime);
% TIME = dt:dt:FinalTime;
TIME = FinalTime;


nrtol  = 1e-10;
niter = 1e3;
figure(1)
subplot(211)
plot(x(:),h)
xline(gate,'k--')

cl = sqrt(g*hl);
cr = sqrt(g*hr);

alpha = 0.0*pi;

mcells = length(x(:));

pack_analytical_parameters


dcrit = (ur - ul) - 2.0*(cl + cr);

if(hl<=0 || hr<=0 || dcrit>=0)
filename = [plotpath,'x0=',num2str(gate),'_T=',num2str(FinalTime),'_hL=',num2str(hl),'_hR=',num2str(hr),'_uL=',num2str(ul),'_uR=',num2str(ur),'_alpha=',num2str(alpha),'dry_analytical.gif'];
'Computing dry-bed solution'

else
filename = [plotpath,'x0=',num2str(gate),'_T=',num2str(FinalTime),'_hL=',num2str(hl),'_hR=',num2str(hr),'_uL=',num2str(ul),'_uR=',num2str(ur),'_alpha=',num2str(alpha),'wet_analytical.gif'];
'Computing wet-bed solution'

end

clear mov
for i = 1:length(TIME)
   analytical_param.timout = TIME(i); 

    if(hl<=0 || hr<=0 || dcrit>=0)
        [h,u] = drybed(h,u,analytical_param);
%         plot_title = ['Dry',' x_0=',num2str(gate),' T=',num2str(FinalTime),' h_L=',num2str(hl),' h_R=',num2str(hr),' u_L=',num2str(ul),' u_R=',num2str(ur), ' \alpha = ',num2str(alpha),', t = ',num2str(TIME(i))];

    else
        [h,u,converged_flag] = wetbed(h,u,niter,analytical_param);
%         plot_title = ['Wet',' x_0=',num2str(gate),' T=',num2str(FinalTime),' h_L=',num2str(hl),' h_R=',num2str(hr),' u_L=',num2str(ul),' u_R=',num2str(ur), ' \alpha = ',num2str(alpha),', t = ',num2str(TIME(i))];
    end

    u = u + alpha*g*TIME(i);
    x = x + 0.5*g*alpha*TIME(i)^2;

    figure(2)
%     sgtitle(plot_title)
    subplot(211)
    plot(x(:),h)
    xline(gate,'k--')
%     ylim([0,4])
    ylabel('h')

    subplot(212)
    plot(x(:),u)
    xline(gate,'k--')
%     ylim([0,12])
    ylabel('u')
    xlabel('x')
    set(gcf,'color','w')
    drawnow
    
    mov(i) = getframe(gcf);
end

movie2gif(mov, filename, 'LoopCount', 3, 'DelayTime', 0)
