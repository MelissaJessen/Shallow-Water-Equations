% Simplified Linearized Shallow Water Equations on a Sphere (FVM) - 2D Version
clear; clc;

% Parameters
L_theta = 2*pi; % Length of the domain in theta (0 to 2*pi in longitude)
L_phi = pi;     % Length of the domain in phi (0 to pi in latitude)
N_theta = 200;  % Number of control volumes in theta
N_phi = 200;    % Number of control volumes in phi
dtheta = L_theta / N_theta; % Width of each control volume in theta
dphi = L_phi / N_phi;       % Width of each control volume in phi
theta = linspace(0, L_theta, N_theta+1); % Cell edges in theta
phi = linspace(0, L_phi, N_phi+1);       % Cell edges in phi
theta_center = theta(1:end-1) + dtheta/2; % Cell centers in theta
phi_center = phi(1:end-1) + dphi/2;      % Cell centers in phi
dt = 0.0025; % Time step
T = 1; % Total simulation time
a = 1; % Radius of the sphere (for simplicity)
g = 9.81; % Gravitational acceleration
h0 = 10; % Mean layer depth in meters

% Initial conditions: Gaussian pulse for h' and zero velocity
sigma_vals = [pi/8, pi/16, pi/32]; % Width of the Gaussian
sigma_no = 2;
sigma = sigma_vals(sigma_no);
theta0 = pi / 4; % Center of the Gaussian
phi0 = pi / 2;   % Center of the Gaussian
h_prime = exp(- ( (theta_center - theta0).^2/(2*sigma.^2) ) + ( (phi_center - phi0).^2/(2*sigma.^2))  ); % Gaussian 2D
%h_prime = exp(-((theta_center - theta0) / sigma).^2 - ((phi_center - phi0) / sigma).^2); % Initial perturbation in wave elevation
v_theta = zeros(size(theta_center)); % Initial velocity in theta direction
v_phi = zeros(size(phi_center)); % Initial velocity in phi direction

% Prepare the spherical coordinates for plotting
[Theta, Phi] = meshgrid(theta_center, phi_center); % 2D grid of theta and phi
X = a * sin(Phi) .* cos(Theta); % X-coordinates on the sphere
Y = a * sin(Phi) .* sin(Theta); % Y-coordinates on the sphere
Z = a * cos(Phi);               % Z-coordinates on the sphere

% Since h_prime is 2D (phi x theta), we replicate it to match the meshgrid
h_prime = exp(-((Theta - theta0).^2 / (2 * sigma.^2)) - ((Phi - phi0).^2 / (2 * sigma.^2)));

% Define a custom colormap as an Nx3 matrix (RGB values)
custom_cmap = [0.5 0.7 0.89;  % Blue
               0 0 1;  % Mostly blue
               0 0 0.5]; % Red

% % Plot the data on the sphere
% figure;
% surf(X, Y, Z, h_prime, 'EdgeColor', 'none'); % Map h' values as color
% colormap(winter); % Set colormap
% colormap(interp1([1, 2, 3], custom_cmap, linspace(1, 3, 256)));
% colorbar; % Add colorbar for reference
% title("Initial water height h (m) on a Sphere");
% axis equal;
% xlabel('X'); ylabel('Y'); zlabel('Z');
% view(3); % 3D view


%
% Exact solution (for comparison)
c = sqrt(g * h0); % Wave speed
uexact = @(t) exp(-((mod(theta_center - c*t, L_theta) - theta0) / sigma).^2 - ((mod(phi_center - c*t, L_phi) - phi0) / sigma).^2);
h_exact = @(t) exp(-((mod(theta_center - c*t, L_theta) - theta0) / sigma).^2 - ((mod(phi_center - c*t, L_phi) - phi0) / sigma).^2); % Exact wave profile

% Flux computation function
compute_flux = @(q_left, q_right) 0.5 * (q_left + q_right); % Central flux

% Time-stepping loop
time = 0;
%h_storage = h_prime; % Storage for wave elevation visualization
t_storage = 0;

% Initalize empty storage
h_storage = [];

while time < T
    % Compute fluxes for h', v_theta, and v_phi
    h_flux_theta = compute_flux(h_prime, circshift(h_prime, -1, 2)); % Flux in the theta direction
    h_flux_phi = compute_flux(h_prime, circshift(h_prime, -1, 1)); % Flux in the phi direction
    
    v_flux_theta = compute_flux(v_theta, circshift(v_theta, -1, 2)); % Flux for v_theta
    v_flux_phi = compute_flux(v_phi, circshift(v_phi, -1, 1)); % Flux for v_phi
    
    % Compute finite volume updates (flux differences)
    k_h1 = -h0 / a * (v_flux_theta - circshift(v_flux_theta, 1, 2)) / dtheta - h0 / a * (v_flux_phi - circshift(v_flux_phi, 1, 1)) / dphi;
    k_v_theta1 = -g * (h_flux_theta - circshift(h_flux_theta, 1, 2)) / dtheta;
    k_v_phi1 = -g * (h_flux_phi - circshift(h_flux_phi, 1, 1)) / dphi;

    % Stage 2
    h_flux_theta = compute_flux(h_prime + 0.5 * dt * k_h1, circshift(h_prime + 0.5 * dt * k_h1, -1, 2));
    h_flux_phi = compute_flux(h_prime + 0.5 * dt * k_h1, circshift(h_prime + 0.5 * dt * k_h1, -1, 1));
    
    v_flux_theta = compute_flux(v_theta + 0.5 * dt * k_v_theta1, circshift(v_theta + 0.5 * dt * k_v_theta1, -1, 2));
    v_flux_phi = compute_flux(v_phi + 0.5 * dt * k_v_phi1, circshift(v_phi + 0.5 * dt * k_v_phi1, -1, 1));
    
    k_h2 = -h0 / a * (v_flux_theta - circshift(v_flux_theta, 1, 2)) / dtheta - h0 / a * (v_flux_phi - circshift(v_flux_phi, 1, 1)) / dphi;
    k_v_theta2 = -g * (h_flux_theta - circshift(h_flux_theta, 1, 2)) / dtheta;
    k_v_phi2 = -g * (h_flux_phi - circshift(h_flux_phi, 1, 1)) / dphi;

    % Stage 3
    h_flux_theta = compute_flux(h_prime + 0.5 * dt * k_h2, circshift(h_prime + 0.5 * dt * k_h2, -1, 2));
    h_flux_phi = compute_flux(h_prime + 0.5 * dt * k_h2, circshift(h_prime + 0.5 * dt * k_h2, -1, 1));
    
    v_flux_theta = compute_flux(v_theta + 0.5 * dt * k_v_theta2, circshift(v_theta + 0.5 * dt * k_v_theta2, -1, 2));
    v_flux_phi = compute_flux(v_phi + 0.5 * dt * k_v_phi2, circshift(v_phi + 0.5 * dt * k_v_phi2, -1, 1));
    
    k_h3 = -h0 / a * (v_flux_theta - circshift(v_flux_theta, 1, 2)) / dtheta - h0 / a * (v_flux_phi - circshift(v_flux_phi, 1, 1)) / dphi;
    k_v_theta3 = -g * (h_flux_theta - circshift(h_flux_theta, 1, 2)) / dtheta;
    k_v_phi3 = -g * (h_flux_phi - circshift(h_flux_phi, 1, 1)) / dphi;

    % Stage 4
    h_flux_theta = compute_flux(h_prime + dt * k_h3, circshift(h_prime + dt * k_h3, -1, 2));
    h_flux_phi = compute_flux(h_prime + dt * k_h3, circshift(h_prime + dt * k_h3, -1, 1));
    
    v_flux_theta = compute_flux(v_theta + dt * k_v_theta3, circshift(v_theta + dt * k_v_theta3, -1, 2));
    v_flux_phi = compute_flux(v_phi + dt * k_v_phi3, circshift(v_phi + dt * k_v_phi3, -1, 1));
    
    k_h4 = -h0 / a * (v_flux_theta - circshift(v_flux_theta, 1, 2)) / dtheta - h0 / a * (v_flux_phi - circshift(v_flux_phi, 1, 1)) / dphi;
    k_v_theta4 = -g * (h_flux_theta - circshift(h_flux_theta, 1, 2)) / dtheta;
    k_v_phi4 = -g * (h_flux_phi - circshift(h_flux_phi, 1, 1)) / dphi;

    % Update h', v_theta, v_phi
    h_prime = h_prime + dt * (1/6 * k_h1 + 1/3 * k_h2 + 1/3 * k_h3 + 1/6 * k_h4);
    v_theta = v_theta + dt * (1/6 * k_v_theta1 + 1/3 * k_v_theta2 + 1/3 * k_v_theta3 + 1/6 * k_v_theta4);
    v_phi = v_phi + dt * (1/6 * k_v_phi1 + 1/3 * k_v_phi2 + 1/3 * k_v_phi3 + 1/6 * k_v_phi4);

    % Update time
    time = time + dt;
    
    % Store wave elevation for visualization
    %h_storage = [h_storage; h_prime];
    h_storage(:,:,end+1) = h_prime; % Append h_prime to the 3rd dimension
    t_storage = [t_storage; time];
end


% Visualization
figure;
for k = 1:size(t_storage, 1)
    surf(X, Y, Z, h_storage(:, :, k), "EdgeColor", "none")
    colormap(winter);
    title(['Time: ', num2str(t_storage(k), '%.2f'), ' s']);
    axis equal;
    view([120, -10]) % View
    drawnow;
end

