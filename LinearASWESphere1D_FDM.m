% Simplified Linearized Shallow Water Equations on a Great Circle
clear; clc;

% Parameters
L = 2*pi; % Length of the domain (0 to 2*pi in spherical coordinates)
N_theta = 500; % Number of grid points
dtheta = L / N_theta; % Grid spacing
theta = linspace(0, L, N_theta+1); % Grid points
theta = theta(1:end-1); % Remove duplicate point for periodicity
dt = 0.0025; % Time step
T = 10; % Total simulation time
a = 1;%6371e3; % Radius of Earth in meters
phi = pi / 6; % Latitude (fixed)
g = 9.81; % Gravitational acceleration
h0 = 10; % Mean layer depth in meters
f = 2 * 7.2921e-5 * sin(phi); % Coriolis parameter

% Initial conditions: Gaussian pulse for h' and zero velocity
sigma = pi / 16; % Width of the Gaussian
theta0 = pi / 4; % Center of the Gaussian
h_prime = exp(-((theta - theta0) / sigma).^2); % Initial perturbation in wave elevation
v = zeros(size(theta)); % Initial velocity

% Exact solution (for comparison, assumes no Coriolis force)
c = sqrt(g * h0); % Wave speed
uexact = @(t) exp(-((mod(theta - c*t, L) - theta0) / sigma).^2);

% Central finite difference operator
D = @(u) (circshift(u, -1) - circshift(u, 1)) / (2 * dtheta);

% Time-stepping loop
time = 0;
h_storage = h_prime; % Storage for wave elevation visualization
t_storage = 0;

while time < T
    % Initialize Runge-Kutta increments
    k_h1 = -h0 / (a * cos(phi)) * D(v);
    k_v1 = -g * D(h_prime) - f * v;

    k_h2 = -h0 / (a * cos(phi)) * D(v + 0.5 * dt * k_v1);
    k_v2 = -g * D(h_prime + 0.5 * dt * k_h1) - f * (v + 0.5 * dt * k_v1);

    k_h3 = -h0 / (a * cos(phi)) * D(v + 0.5 * dt * k_v2);
    k_v3 = -g * D(h_prime + 0.5 * dt * k_h2) - f * (v + 0.5 * dt * k_v2);

    k_h4 = -h0 / (a * cos(phi)) * D(v + dt * k_v3);
    k_v4 = -g * D(h_prime + dt * k_h3) - f * (v + dt * k_v3);

    % Update h' and v
    h_prime = h_prime + dt * (1/6 * k_h1 + 1/3 * k_h2 + 1/3 * k_h3 + 1/6 * k_h4);
    v = v + dt * (1/6 * k_v1 + 1/3 * k_v2 + 1/3 * k_v3 + 1/6 * k_v4);

    % Update time
    time = time + dt;

    % Store wave elevation for visualization
    h_storage = [h_storage; h_prime];
    t_storage = [t_storage; time];
end

% Visualization of the wave elevation propagating on the great circle
figure;
r0 = 1.5; % Base radius of the circle
tim = 0; % Initial time for visualization

for k = 1:size(h_storage, 1)
    % Plot the wave elevation on a circle (spherical projection)
    plot(r0*cos(theta), -r0*sin(theta), 'k'); % Base circle
    hold on;
    %plot((r0 + uexact(tim)).*cos(theta), -(r0 + uexact(tim)).*sin(theta), 'b', 'LineWidth', 1.5); % Exact solution
    plot((r0 + h_storage(k, :)).*cos(theta), -(r0 + h_storage(k, :)).*sin(theta), 'r', 'LineWidth', 2); % Numerical solution
    hold off;
    title(['Time: ', num2str(t_storage(k), '%.2f'), ' s']);
    axis equal;
    axis([-2.5, 2.5, -2.5, 2.5]);
    %legend('Base Circle', 'Exact Solution', 'Numerical Solution');
    drawnow;
    tim = tim + dt; % Update time for exact solution
end
