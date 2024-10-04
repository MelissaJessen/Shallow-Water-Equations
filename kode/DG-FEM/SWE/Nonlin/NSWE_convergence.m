clear,clc,close all
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\utilities');
addpath('C:\Users\Bruger\Documents\GitHub\Thesis\matlab\DG-FEM\SWE\Nonlin\analytical_solver')
plotpath = 'C:\Users\Bruger\OneDrive - Danmarks Tekniske Universitet\Thesis\Report\fig\DG-FEM\SWE\nonlin\';


NODETOL = 1e-10;
alpha = 0;
periodic = 0;

N = 2;
K = 100;
FinalTime = 2.5;


xmin = 0.0;
xmax = 50;
A = 1;
g = 9.81;
gate = 20;


initialize_solver
hl = 1;
hr = 1;
ul = -1;
ur = 1;
h  = zeros(size(x));
u = zeros(size(x));

h(x<=gate) = hl;
h(x>gate) = hr;


u(x<=gate) = ul;
u(x>gate) = ur;

hu = u.*h;

% Pack structs
struct_pack
param.g = g;
param.xmin   = xmin;
param.xmax   = xmax;


% Pack analytical parameters

nrtol  = 1e-10;
niter = 1e3;
cl = sqrt(g*hl);
cr = sqrt(g*hr);
alpha = 0.0*pi;
mcells = length(x(:));
chalen = xmax;
timout = FinalTime;

pack_analytical_parameters
dcrit = (ur - ul) - 2.0*(cl + cr);



if(hl<=0 || hr<=0 || dcrit>=0)
    'Dry'
    [htrue,utrue] = drybed(h,u,analytical_param);
    plot_title = ['Dry',' x_0=',num2str(gate),' T=',num2str(FinalTime),' h_L=',num2str(hl),' h_R=',num2str(hr),' u_L=',num2str(ul),' u_R=',num2str(ur), ' \alpha = ',num2str(alpha),', t = ',num2str(FinalTime)];

else
    'Wet'
    [htrue,utrue,converged_flag] = wetbed(h,u,niter,analytical_param);
    plot_title = ['Wet',' x_0=',num2str(gate),' T=',num2str(FinalTime),' h_L=',num2str(hl),' h_R=',num2str(hr),' u_L=',num2str(ul),' u_R=',num2str(ur), ' \alpha = ',num2str(alpha),', t = ',num2str(FinalTime)];
end

tic
[h,hu,ht,hut,tt] = NSWE1D(h,hu,x,FinalTime,1,param,maps);
toc

figure(1)
subplot(211)
hold on
plot(x,h)
plot(x,htrue,'k--')
hold off
subplot(212)
hold on
plot(x,hu./h,'k--')
plot(x,utrue)


%% K and P convergence

% Polynomial konvergence
nk = [6,7,8];
Karray = 2.^nk;
Narray = [1,2];
% 
% k_line = 1./(Karray.^Narray);

nidx       = 0;
kidx       = 0;
error = zeros(length(Karray),length(Narray));
k_lines    = [];
alpha = 1;
for K = Karray
    kidx = kidx + 1;
    nidx = 0;
    for N = Narray
        nidx = nidx + 1;

        initialize_solver
        h  = zeros(size(x));
        u = zeros(size(x));

        h(x<=gate) = hl;
        h(x>gate) = hr;


        u(x<=gate) = ul;
        u(x>gate) = ur;
        hu = h.*u;
        mcells = length(x(:));

        struct_pack
        pack_analytical_parameters

        [h,hu] = NSWE1D(h,hu,x,FinalTime,1,param,maps);
        [htrue,utrue,converged_flag] = wetbed(h,u,niter,analytical_param);

        q = [u,hu./u];
        qtrue = [htrue,utrue];
        error(kidx,nidx) = norm(qtrue-q,2)/(mcells);

        k_lines = [k_lines,1./Karray'.^(N+1)];         
    end
%     
end
%%
figure
loglog(Karray,error,'k-*')
hold on
loglog(Karray,error,'k-o')
loglog(Karray,k_lines,'k--')
hold off
legend('Upwind','Central')
grid('minor')
xlabel('K')
ylabel('E = ||q-q_{true}||_\infty')
title('Non-periodic Domain')
%%

figure(1)
subplot(211)
hold on
plot(x,h)
plot(x,htrue,'k--') 
hold off
subplot(212)
hold on
plot(x,hu./h,'k--')
plot(x,utrue)

