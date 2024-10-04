% Semi-Impicit scheme for the coupled 1D Shallow Water and 2D Richards 
% equations using the nested Newton-type method of Casulli&Zanolli

clear;
clc
global alpha thetas thetar n m Ks psic 
global dx dz dt K Kx Kz Hx Nx Nz dTheta12 g gamma

% Physical model parameters in SI units 
day     = 24*3600;    
g       = 9.81;         % m/s^2 
Ks      = 0.062/day;    % [m/s] 
thetas  = 0.41;         % [-] saturated water content 
thetar  = 0.095;        % [-] residual water content 
n       = 1.31;         % [-] parameter n 
m       = 1 - 1/n;      % [-] parameter m 
alpha   = 1.9;          % [m^(-1)] 
gamma   = 1e-2;         % friction coefficient 
psic    =-1/alpha*( (n-1)/n )^(1/n); % critical value of psi where the maximum of the C = dTheta/dPsi is located 
% Domain 
xL   = -1000; % left
xR   = +1000; % right 
zL   = -2;    % bottom 
zR   = 0;     % surface
Nx   = 100;   % horizontal 
Nz   = 20;    % vertical 
dx   = (xR-xL)/Nx; % x mesh spacing 
dz   = (zR-zL)/Nz; % y mesh spacing 
x    = linspace(xL+dx/2,xR-dx/2,Nx); 
z    = linspace(zL+dz/2,zR-dz/2,Nz); 
tend = 30000;   % set the final simulation time 
t    = 0;     % initial time 

% set the initial conditions
psi = zeros(Nz+1,Nx);
psi(1:Nz,:)     = -1;   % pressure head in the subsurface
psi(Nz+1,40:60) = 0.5;  % eta in the free surface layer (SWE layer)

u    = zeros(1,Nx+1);       % velocity at the x faces of the cells, only in the free-surface layer
bath = zeros(1,Nx);         % bathymetry 

K  = zeros(Nz+1,Nx);    % hydraulic conductivity at the cell centersin
Kx = zeros(Nz+1,Nx+1);  % hydraulic conductivity at the cell interfaces in x-direction
Kz = zeros(Nz+2,Nx);    % hydraulic conductivity at the cell interfaces in z-direction
Hx = zeros(1,Nx+1);     % free surface height at the cell interfaces
rhs= zeros(Nz+1,Nx);    % right hand side of the mildly nonlinear system

theta    = zeros(Nz+1,Nx);  % volume fraction of the fluid
dTheta12 = zeros(Nz+1,Nx);  %

% time loop
for nt = 1:10000000
    dt = 60;
    if (t +dt > tend)
        dt = tend - t;
    end
    if (t >= tend)
        break;
    end

    % *** Explitit step
    ustar = u;

    % *** Implicit step (Nested Newton mehtod)
    [theta,rhs] = LinearPartCoeff(psi,theta,bath,ustar);

    tol = 1e-7;
    maxNewton = 100;
    % ---------------- Outer iterations
    psi(1:Nz,1:Nx) = min(psic,psi(1:Nz,1:Nx));  % initial guess for the Outer iterations
    for iOut = 1:maxNewton
        [ResOut,NormResOut] = ResidualOuter(psi,bath,rhs);
        disp(strcat(' Outer iteration = ', int2str(iOut), ', outres = ', num2str(NormResOut))); 
        if ( NormResOut < tol )
            break;
        end
        % ---------------- Inner iterations
        psi_alpha = psi;    % save the value of psi at the current outer iteration (alpha in the w1_d5.pdf)
        % psi = max(psic,psi);  % initial guess for the inner iterations as suggested by the paper by Casulli&Zanolli
        for iInn = 1:maxNewton
            [ResInn,NormResInn] = ResidualInner(psi,psi_alpha,bath,rhs);
            disp(strcat('  - Inner iteration = ',int2str(iInn), ', inres = ', num2str(NormResInn))) 
            if ( NormResInn < tol )
                break;
            end
            [dpsi,CGerr,CGk] = CGsolver(ResInn,@MatVecProd_psi);
            psi = psi - dpsi;
        end
    end
    % Now, update the velocity
    u(1   ) = 0;    % wall BC on the left
    u(2:Nx) = ( ustar(2:Nx) - g*dt/dx*( psi(Nz+1,2:Nx) - psi(Nz+1,1:Nx-1) ) )./( 1 + dt*gamma./Hx(2:Nx) );
    u(Nx+1) = 0;    % wall BC on the right

    % time update
    t = t + dt;

    % postprocessing and plotting part
    H   = theta(Nz+1,:);    % the last z-layer of theta is H = max(0, eta-b)
    eta = H + bath;
    subplot(3,1,1)
    plot(x,eta,'o-')
    axis([xL xR -0.1 0.7])
    xlabel('$x$','Interpreter','latex')
    ylabel('\eta')
    title(strcat('$t = $',num2str(t)),"Interpreter","latex") 

    subplot(3,1,2) 
    surf(x,z,psi(1:Nz,:),'EdgeColor','none','FaceColor','interp')   
    xlabel('$x$',Interpreter='latex')
    ylabel('$z$',Interpreter='latex')
    zlabel('\psi')
    view([0 90])
    title('\psi')
    colorbar;
    
    subplot(3,1,3) 
    surf(x,z,Theta(psi(1:Nz,:)),'EdgeColor','none','FaceColor','interp')   
    xlabel('$x$','Interpreter','latex')
    ylabel('$z$','Interpreter','latex')
    zlabel('\theta(\psi)')
    view([0 90])
    title('\theta(\psi)') 
    colorbar;
    
    pause(0.0001) 


end

















