% Godunov-type Finite volume method for the 1D Shallow Water Equation
clear;
clc;

global Nx g;

% Physical parameters
xL =-1;
xR = 1;
t = 0;
tend = 1;
g = 9.81;   % acceleration due to gravity

% numerical parameters
Nx = 100;           % the number of the cell-centers
dx = (xR - xL)/Nx;  % distance between the mesh points
CFL= 0.9;           % CFL number <= 1
x  = linspace(xL+dx/2,xR-dx/2,Nx);

% Initial conditoins
q = zeros(2,Nx);
hL = 1;   
hR = 0.1;
q(1,1:Nx/2)    = hL;
q(1,Nx/2+1:Nx) = hR;

% Define parameters for the Gaussian
h0 = 1;          % Amplitude of the Gaussian
x0 = (xL + xR)/2; % Mean (center) of the Gaussian, centered in the domain
sigma = 0.1;

% Create the Gaussian blob for h(x)
h = h0 * exp(-(x - x0).^2 / (2 * sigma^2));
q(1,:) = h;

qnew = q;

% plot the initial data
subplot(2,1,1)
plot(x,q(1,:),'-o')
ylabel('h')
xlabel('x')

subplot(2,1,2)
plot(x,q(2,:)./q(1,:),'-o')
ylabel('u')
xlabel('x')

% compute the total amount of water content in the domain (integral of h
% between xL and xR):
H0 = sum(dx*q(1,:));

% time loop
for n = 1:1000000
    amax = max(max(abs(Lambda(q))));
    dt   = CFL*dx/amax; 
    if (t + dt > tend)
        dt = tend - t;
    end
    if (t >= tend)
        break;  % stop the time loop
    end
    
    % space loop
    dtdx = dt/dx;
    for i = 1:Nx
        if    ( i == 1  )
            q_ghost    = q(:,i);    % "create" a ghost state outside of the domain
            q_ghost(2) =-q(2,i);    % invert the sing of its momentum/velocity
            qm = ExactRP(q_ghost,q(:,i  ),0);  % compute the exact sol of the RP at x_{i-1/2}
            fm = f(qm);
            qp = ExactRP(q(:,i  ),q(:,i+1),0);  % compute the exact sol of the RP at x_{i+1/2}
            fp = f(qp);
            % stadrand Finite-Volume update:
            qnew(:,i) = q(:,i) - dtdx*( fp - fm );
        elseif( i == Nx )
            q_ghost    = q(:,i);    % "create" a ghost state outside of the domain
            q_ghost(2) =-q(2,i);    % invert the sing of its momentum/velocity
            qm = ExactRP(q(:,i-1),q(:,i  ),0);  % compute the exact sol of the RP at x_{i-1/2}
            fm = f(qm);
            qp = ExactRP(q(:,i  ),q_ghost,0);  % compute the exact sol of the RP at x_{i+1/2}
            fp = f(qp);
            % stadrand Finite-Volume update:
            qnew(:,i) = q(:,i) - dtdx*( fp - fm );
        else
            qm = ExactRP(q(:,i-1),q(:,i  ),0);  % compute the exact sol of the RP at x_{i-1/2}
            fm = f(qm);
            qp = ExactRP(q(:,i  ),q(:,i+1),0);  % compute the exact sol of the RP at x_{i+1/2}
            fp = f(qp);
            % stadrand Finite-Volume update:
            qnew(:,i) = q(:,i) - dtdx*( fp - fm );
        end

    end

    % update time and solution
    t = t + dt;
    q = qnew;

    % plot the snumerical solution
    subplot(2,1,1)
    plot(x,q(1,:),'-o')
    ylabel('h')
    xlabel('x')
    title(strcat('time = ',num2str(t)))
    
    subplot(2,1,2)
    plot(x,q(2,:)./q(1,:),'-o')
    ylabel('u')
    xlabel('x')
    errMassCons = sum(dx*q(1,:)) - H0;  % compare integral of h with the one at t=0
    title(strcat('mass cons error = ',num2str(errMassCons)))

    pause(0.001)
    % for monitoring: Compute the amount of water and compare with H0

end





