% Simplified Linearized Shallow Water Equations on a Great Circle (FVM)
clear; clc;

% Parameters
L = 2*pi; % Length of the domain (0 to 2*pi in spherical coordinates)
N_theta = 500; % Number of control volumes
dtheta = L / N_theta; % Width of each control volume
theta = linspace(0, L, N_theta+1); % Cell edges
theta_center = theta(1:end-1) + dtheta/2; % Cell centers
dt = 0.0025; % Time step
T = 0.1; % Total simulation time
a = 1;%6371e3; % Radius of Earth in meters
phi = pi / 10; % Latitude (fixed)
g = 9.81; % Gravitational acceleration
h0 = 10; % Mean layer depth in meters
f = 2 * 7.2921e-5 * sin(phi); % Coriolis parameter

% Initial conditions: Gaussian pulse for h' and zero velocity
sigma_vals = [pi/8, pi/16, pi/32]; % Width of the Gaussian
sigma_no = 3;
sigma = sigma_vals(sigma_no);
%sigma = pi / 16; 
theta0 = pi / 4; % Center of the Gaussian
h_prime = exp(-((theta_center - theta0) / sigma).^2); % Initial perturbation in wave elevation
v = zeros(size(theta_center)); % Initial velocity

% Exact solution (for comparison, assumes no Coriolis force)
c = sqrt(g * h0); % Wave speed
uexact = @(t) exp(-((mod(theta_center - c*t, L) - theta0) / sigma).^2);
h_exact = @(t) exp(-((mod(theta_center - c*t, L) - theta0) / sigma).^2); % Exact wave profile

% Flux computation function 
compute_flux = @(q_left, q_right) 0.5 * (q_left + q_right); % Upwind or central flux

% Time-stepping loop
time = 0;
h_storage = h_prime; % Storage for wave elevation visualization
t_storage = 0;

while time < T
    % Compute fluxes for h' and v
    h_flux = compute_flux(h_prime, circshift(h_prime, -1));
    v_flux = compute_flux(v, circshift(v, -1));
    
    % Compute finite volume updates (flux differences)
    k_h1 = -h0 / (a * cos(phi)) * (v_flux - circshift(v_flux, 1)) / dtheta;
    k_v1 = -g / (a * cos(phi)) * (h_flux - circshift(h_flux, 1)) / dtheta;

    % Stage 2
    h_flux = compute_flux(h_prime + 0.5 * dt * k_h1, circshift(h_prime + 0.5 * dt * k_h1, -1));
    v_flux = compute_flux(v + 0.5 * dt * k_v1, circshift(v + 0.5 * dt * k_v1, -1));
    k_h2 = -h0 / (a * cos(phi)) * (v_flux - circshift(v_flux, 1)) / dtheta;
    k_v2 = -g / (a * cos(phi)) * (h_flux - circshift(h_flux, 1)) / dtheta;
    
    % Stage 3
    h_flux = compute_flux(h_prime + 0.5 * dt * k_h2, circshift(h_prime + 0.5 * dt * k_h2, -1));
    v_flux = compute_flux(v + 0.5 * dt * k_v2, circshift(v + 0.5 * dt * k_v2, -1));
    k_h3 = -h0 / (a * cos(phi)) * (v_flux - circshift(v_flux, 1)) / dtheta;
    k_v3 = -g / (a * cos(phi))* (h_flux - circshift(h_flux, 1)) / dtheta;
    
    % Stage 4
    h_flux = compute_flux(h_prime + dt * k_h3, circshift(h_prime + dt * k_h3, -1));
    v_flux = compute_flux(v + dt * k_v3, circshift(v + dt * k_v3, -1));
    k_h4 = -h0 / (a * cos(phi)) * (v_flux - circshift(v_flux, 1)) / dtheta;
    k_v4 = -g / (a*cos(phi)) * (h_flux - circshift(h_flux, 1)) / dtheta;

    % Update h' and v
    h_prime = h_prime + dt * (1/6 * k_h1 + 1/3 * k_h2 + 1/3 * k_h3 + 1/6 * k_h4);
    v = v + dt * (1/6 * k_v1 + 1/3 * k_v2 + 1/3 * k_v3 + 1/6 * k_v4);

    % Update time
    time = time + dt;

    % Store wave elevation for visualization
    h_storage = [h_storage; h_prime];
    t_storage = [t_storage; time];
end

% Visualization
figure;
r0 = 1.5; % Base radius of the circle
for k = 1:size(h_storage, 1)
    h_exact_k = h_exact(t_storage(k)); % Compute exact solution at time t_storage(k)

    plot(r0*cos(theta_center), -r0*sin(theta_center), 'k'); % Base circle
    hold on;

    plot((r0 + h_storage(k, :)).*cos(theta_center), -(r0 + h_storage(k, :)).*sin(theta_center), 'r', 'LineWidth', 2); % Numerical solution
    hold off;
    title(['Time: ', num2str(t_storage(k), '%.2f'), ' s']);
    axis equal;
    axis([-2.5, 2.5, -2.5, 2.5]);
    drawnow;
end

%% Save data to mat.file
% Define folder

folder = 'Shallow-Water-Equations\dataFNO';
filename = ['linear-SWE-sphere-1d_sigma=',num2str(sigma_no),',T=',num2str(T)];
fullpath = fullfile(folder, filename);
save(fullpath,'h_storage', 't_storage', 'theta');


%% Plot the initial condition

k = 1;

figure;
plot(r0*cos(theta_center), -r0*sin(theta_center), 'k'); % Base circle
hold on;
plot((r0 + h_storage(k, :)).*cos(theta_center), -(r0 + h_storage(k, :)).*sin(theta_center), 'r', 'LineWidth', 2); % Numerical solution
hold off;
%title(['Inital conditions for sigma at time: ', num2str(t_storage(k), '%.2f'), ' s']);
title(['Initial conditions for \sigma = \pi/32']);
%title(['Inital conditions for sigma : ', num2str(sigma, '%.2f'), ' s']);
axis equal;
axis([-2.5, 2.5, -2.5, 2.5]);
drawnow;

% Save the figure as PDF
%exportgraphics(gcf, 'C:/Users/Matteo/Shallow-Water-Equations/plots/SWE-spherical-1d-initial_conditions_sigma3.pdf', 'ContentType', 'vector');
 