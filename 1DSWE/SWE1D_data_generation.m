% Godunov-type Finite volume method for the 1D Shallow Water Equation
clear;
clc;

% Start timing;
tic;

global Nx g;
   
% Physical parameters
xL = 0;
xR = 1;
t = 0;
tend = 1;
g = 9.81;   % acceleration due to gravity

% numerical parameters
Nx = 200;           % the number of the cell-centers
dx = (xR - xL)/Nx;  % distance between the mesh points
CFL= 0.9;           % CFL number <= 1
x  = linspace(xL+dx/2,xR-dx/2,Nx);

% Initial conditions
q = zeros(2,Nx);
u = q(2,:);

% Sigma values
%sigma_vals = linspace(0.1, 1, 10);
sigma_vals = 0.1;

mu_vals = linspace(0.3,0.8,6); % Mean of the Gaussian

% Define folder
folder = 'C:\Users\Matteo\Shallow-Water-Equations\dataFNO';

% Loop over each sigma value
for s = 1:length(sigma_vals)
    sigma = sigma_vals(s); % Standard deviation (width) of the Gaussian

    % Print current sigma value
    fprintf('Current sigma: %.2f\n', sigma);

    % Define parameters for the Gaussian
    h0 = 1;          % Amplitude of the Gaussian
    
    % Loop over each mu value
    for m = 1:length(mu_vals)
        mu = mu_vals(m);

        % Print current mu value
        fprintf('Current mu: %.2f\n', mu);
   
        % Initial conditoins
        q = zeros(2,Nx);
        u = q(2,:);

        % Create the Gaussian blob for h(x)
        h = h0 * exp(-(x - mu).^2 / (2 * sigma^2));
        q(1,:) = h;
        qnew = q;

        % compute the total amount of water content in the domain (integral of h between xL and xR):
        H0 = sum(dx*q(1,:));

        % Initialize storage for h and u at each timestep
        h_all = h;
        u_all = u;
        time_all = [];

        t = 0; % Reset time for each sigma and mu value

        % Time loop
        for n = 1:1000000
            amax = max(max(abs(Lambda(q))));
            dt   = CFL*dx/amax; 
            if (t + dt > tend)
                dt = tend - t;
            end
            if (t >= tend)
                break;  % stop the time loop
            end
        
            % Space loop
            dtdx = dt/dx;
            for i = 1:Nx
                if    ( i == 1  )
                    q_ghost    = q(:,i);    % "create" a ghost state outside of the domain
                    q_ghost(2) =-q(2,i);    % invert the sing of its momentum/velocity
                    qm = ExactRP(q_ghost,q(:,i  ),0);  % compute the exact sol of the RP at x_{i-1/2}
                    fm = f(qm);
                    qp = ExactRP(q(:,i  ),q(:,i+1),0);  % compute the exact sol of the RP at x_{i+1/2}
                    fp = f(qp);
                    % standard Finite-Volume update:
                    qnew(:,i) = q(:,i) - dtdx*( fp - fm );
                elseif( i == Nx )
                    q_ghost    = q(:,i);    % "create" a ghost state outside of the domain
                    q_ghost(2) =-q(2,i);    % invert the sing of its momentum/velocity
                    qm = ExactRP(q(:,i-1),q(:,i  ),0);  % compute the exact sol of the RP at x_{i-1/2}
                    fm = f(qm);
                    qp = ExactRP(q(:,i  ),q_ghost,0);  % compute the exact sol of the RP at x_{i+1/2}
                    fp = f(qp);
                    % standard Finite-Volume update:
                    qnew(:,i) = q(:,i) - dtdx*( fp - fm );
                else
                    qm = ExactRP(q(:,i-1),q(:,i  ),0);  % compute the exact sol of the RP at x_{i-1/2}
                    fm = f(qm);
                    qp = ExactRP(q(:,i  ),q(:,i+1),0);  % compute the exact sol of the RP at x_{i+1/2}
                    fp = f(qp);
                    % standard Finite-Volume update:
                    qnew(:,i) = q(:,i) - dtdx*( fp - fm );
                end
            end
    
        h = qnew(1,:);
        u = qnew(2,:);
    
        % Append current timestep's h and u to storage
        h_all = [h_all; h];
        u_all = [u_all; u];
        time_all = [time_all; t];
    
        % update time and solution
        t = t + dt;
        q = qnew;

        end

        % Save data for the current sigma and mu value
        filename = ['data_sigma_no',num2str(s),'mu_no',num2str(m)];
        fullpath = fullfile(folder, filename);
        save(fullpath,'time_all','h0','sigma','mu','x','h_all','u_all');
    end
end

% Stop timing
elapsedTime = toc;

% Print the result
fprintf('The program took %.2f seconds to run.\n', elapsedTime);
