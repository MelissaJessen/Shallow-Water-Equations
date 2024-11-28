% Linear Advection on a Great Circle of a Sphere
clear; clc;

% Parameters
L = 2*pi; % Length of the domain (0 to 2*pi in spherical coordinates)
N_theta = 500; % Number of grid points
dtheta = L / N_theta; % Grid spacing
theta = linspace(0, L, N_theta+1); % Grid points
theta = theta(1:end-1); % Remove duplicate point for periodicity
dt = 0.025; % Time step
c = 1; % Wave speed
T = 10; % Total simulation time

% Initial condition: Gaussian pulse
sigma = pi / 16; % Width of the Gaussian
theta0 = pi / 4; % Center of the Gaussian
u = exp(-((theta - theta0) / sigma).^2); % Initial pulse

% Exact solution as a function of time
uexact = @(t) exp(-((mod(theta - c*t, L) - theta0) / sigma).^2);

% Central finite difference operator
D = @(u) (circshift(u, -1) - circshift(u, 1)) / (2 * dtheta);

% ERK4 coefficients
rk4a = [0, 0.5, 0.5, 1]; % Stage weights
rk4b = [1/6, 1/3, 1/3, 1/6]; % Final combination weights

% Time-stepping loop
time = 0;
u_storage = u; % Storage for visualization
t_storage = 0;

while time < T
    % ERK4 Time-stepping
    k1 = -c * D(u);
    k2 = -c * D(u + 0.5 * dt * k1);
    k3 = -c * D(u + 0.5 * dt * k2);
    k4 = -c * D(u + dt * k3);
    
    % Update the solution
    u = u + dt * (rk4b(1) * k1 + rk4b(2) * k2 + rk4b(3) * k3 + rk4b(4) * k4);

    % Enforce periodicity explicitly (optional, usually not required)
    %u(1) = u(end);
    
    % Update time
    time = time + dt;
    
    % Store solution for visualization
    u_storage = [u_storage; u];
    t_storage = [t_storage; time];
end

% Visualization of the pulse propagating on the great circle
figure;
r0 = 1.5; % Base radius of the circle
tim = 0; % Initial time for visualization

for k = 1:size(u_storage, 1)
    % Plot the solution on a circle (spherical projection)
    plot(r0*cos(theta), -r0*sin(theta), 'k'); % Base circle
    hold on;
    plot((r0 + uexact(tim)).*cos(theta), -(r0 + uexact(tim)).*sin(theta), 'b', 'LineWidth', 1.5); % Exact solution
    plot((r0 + u_storage(k, :)).*cos(theta), -(r0 + u_storage(k, :)).*sin(theta), 'r', 'LineWidth', 2); % Numerical solution
    hold off;
    title(['Time: ', num2str(t_storage(k), '%.2f'), ' s']);
    axis equal;
    axis([-2.5, 2.5, -2.5, 2.5]);
    legend('Base Circle', 'Exact Solution', 'Numerical Solution');
    drawnow;
    tim = tim + dt; % Update time for exact solution
end


% % Linear Advection on a Great Circle of a Sphere
% clear; clc;
% 
% % Parameters
% L = 2*pi; % Length of the domain (0 to 2*pi in spherical coordinates)
% N_theta = 500; % Number of grid points
% dtheta = L / N_theta; % Grid spacing
% theta = linspace(0, L, N_theta+1); % Grid points
% theta = theta(1:end-1); % Remove duplicate point for periodicity
% dt = 0.005; % Time step
% c = 1; % Wave speed
% T = 10; % Total simulation time
% 
% % Initial condition: Gaussian pulse
% sigma = pi / 16; % Width of the Gaussian
% theta0 = pi / 4; % Center of the Gaussian
% u = exp(-((theta - theta0) / sigma).^2); % Initial pulse
% 
% % Exact solution as a function of time
% uexact = @(t) exp(-((mod(theta - c*t, L) - theta0) / sigma).^2);
% 
% % Storage for visualization
% u_storage = u;
% t_storage = 0;
% 
% % Time-stepping loop (Finite Difference Scheme)
% time = 0;
% while time < T
%     % Upwind finite difference scheme (periodic boundary)
%     u_new = u - c * dt / dtheta * (u - circshift(u, 1));
% 
%     % Update solution and time
%     u = u_new;
%     time = time + dt;
% 
%     % Store solution for visualization
%     u_storage = [u_storage; u];
%     t_storage = [t_storage; time];
% end
% 
% % Visualization of the pulse propagating on the great circle
% figure;
% r0 = 1.5;
% tim = 0;
% for k = 1:size(u_storage, 1)
%     % Plot the solution on a circle (spherical projection)
%     plot(r0*cos(theta),-r0*sin(theta),'k')
%     hold on
%     plot((r0 + uexact(tim)).*cos(theta),-(r0 + uexact(tim)).*sin(theta),'b')
%     plot((r0 + u_storage(k, :)).*cos(theta),-(r0 + u_storage(k, :)).*sin(theta), 'LineWidth', 2);
%     hold off
%     title(['Time: ', num2str(t_storage(k), '%.2f'), ' s']);
%     %rlim([0 1.1]);
%     axis equal
%     axis([-2.5, 2.5, -2.5, 2.5])
% %    pause(0.025);
% drawnow
% tim = tim + dt;
% end
% 
