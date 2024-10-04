clear,clc,close all

set(0,'defaultfigurecolor',[1 1 1])
set(groot,'DefaultAxesTickLabelInterpreter','Tex');
set(groot, 'DefaultTextInterpreter', 'Tex')
set(groot, 'DefaultLegendInterpreter', 'Tex')
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\analytical_solver');

datapath = 'C:\Users\Bruger\Documents\GitHub\Thesis\data\analytical_data\';

FinalTime = 2.5;

xmin = 0;
xmax = 50;
A = 1;
g = 9.81;

gate = 20;
x = linspace(xmin,xmax,1000)';
dt = 0.01;
t = 0:dt:FinalTime;

hl = 1;
hr = 0.1;
hlt = hl;
hrt = hr;
ul = 2.5;
ur = 0;  
ult = ul;
urt = ur;
h  = zeros(size(x));
u = zeros(size(x));

h(x<=gate) = hl;
h(x>gate) = hr;


u(x<=gate) = ul;
u(x>gate) = ur;
hu = u.*h;

% Pack analytical parameters

nrtol  = 1e-10;
niter = 1e3;
cl = sqrt(g*hl);
cr = sqrt(g*hr);
mcells = length(x(:));
chalen = xmax;
timout = FinalTime;

analytical_param.chalen = chalen;
analytical_param.mcells = mcells;
analytical_param.cl     = cl;
analytical_param.cr     = cr;
analytical_param.dr     = hrt;
analytical_param.dl     = hlt;
analytical_param.ul     = ult;
analytical_param.ur     = urt;

analytical_param.gate   = gate;
analytical_param.g      = g;

analytical_param.tol    = nrtol;
dcrit = (ur - ul) - 2.0*(cl + cr);
h0 = h;
hu0 = hu;

htime = zeros(length(x),length(t));
utime = zeros(length(x),length(t));


for i = 1:length(t)
    analytical_param.timout = t(i);

    if(hl<=0 || hr<=0 || dcrit>=0)
    %     disp(['Dry',' BC: ',num2str(BC),' flux: ',flux])
        [htrue,utrue] = drybed(h,u,analytical_param);
    %     plot_title = ['Dry',' x_0=',num2str(gate),' T=',num2str(FinalTime),' h_L=',num2str(hl),' h_R=',num2str(hr),' u_L=',num2str(ul),' u_R=',num2str(ur), ' \alpha = ',num2str(alpha),', t = ',num2str(FinalTime)];

    else
    %     disp(['Wet',' BC: ',num2str(BC),' flux: ',flux])
        [htrue,utrue,converged_flag] = wetbed(h,u,niter,analytical_param);
    %     plot_title = ['Wet',' x_0=',num2str(gate),' T=',num2str(FinalTime),' h_L=',num2str(hl),' h_R=',num2str(hr),' u_L=',num2str(ul),' u_R=',num2str(ur), ' \alpha = ',num2str(alpha),', t = ',num2str(FinalTime)];
    end
    htime(:,i) = htrue;
    utime(:,i) = utrue;

end
htime = htime';
utime = utime';


[X,T] = meshgrid(x,t);


contourf(X,T,htime)
colorbar

figure
plot(X(end,:),utime(end,:).*htime(end,:))

dat_filename = 'analytical1';
% 
filename = fullfile(datapath,dat_filename);
save(filename,'X','T','htime','utime')

%% GIF
