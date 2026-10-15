% 2D Uniform Plane
clear; clc;
figure(1);
clf;

f = 1e9;                % Frequency 1 GHz
c = 3e8;                % Speed of Light (m/s)
lambda = c / f;         % Wavelength (0.3m)
k = (2 * pi) / lambda;  % Wave Number

Nx = 8;                 % Length of x-axis
Ny = 8;                 % Length of y-axis

dx = lambda/2;
dy = lambda/2;

% Axes for centering in origin 
x_axis = (-(Nx-1)/2:(Nx-1)/2) * dx;
y_axis = (-(Ny-1)/2:(Ny-1)/2) * dy;

% Creates the coordinate system
[x_grid,y_grid] = meshgrid(x_axis,y_axis);
x_coord = x_grid(:);
y_coord = y_grid(:);
z_coord = zeros(size(x_coord));


subplot(1,3,1);
% Plots plane of elements
scatter3(x_coord,y_coord,z_coord,'filled');
hold on;

theta_steer = deg2rad(0);      % Azimuth
phi_steer = deg2rad(90);        % Elevation

% Computes vetcor pointed to desired theat and phi 
u = ones(size(x_coord))*cos(phi_steer)*cos(theta_steer);
v = ones(size(y_coord))*cos(phi_steer)*sin(theta_steer);
w = ones(size(z_coord))*sin(phi_steer);

% Plots vector to each element
quiver3(x_coord,y_coord,z_coord,u,v,w,0.5,'r');
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
axis tight;
zlim([0 1]);
title(['2D Uniform Plane ' num2str(Nx) 'x' num2str(Ny) ...
    ' dx = \lambda/2, dy = \lambda/2, \lambda = ' num2str(lambda) ]);

% Weights appling a phase shift 
weights = exp(1i * k * (x_coord * cos(phi_steer) * cos(theta_steer) + ...
                             y_coord * cos(phi_steer) * sin(theta_steer)));

% Simulate Radiation Pattern
theta_sweep = linspace(-pi, pi, 180);
phi_sweep = linspace(0, pi/2, 90);    % 0 is horizon, pi/2 is zenith
[Theta, Phi] = meshgrid(theta_sweep, phi_sweep);

AF = zeros(size(Theta));   

for i = 1:numel(Theta)
    a_theta_phi = exp(1i * k * (x_coord * cos(Phi(i)) * cos(Theta(i)) + ...
                                y_coord * cos(Phi(i)) * sin(Theta(i))));
                                     
    AF(i) = abs(sum(weights .* conj(a_theta_phi)));
end 

[~, peak_idx] = max(AF(:));

peak_az = rad2deg(Theta(peak_idx));
peak_el = rad2deg(Phi(peak_idx));

fprintf('Peak direction: Az = %.1f°, El = %.1f°\n', ...
    peak_az, peak_el);
    
% Plotting Circular Polar Heatmap 
subplot(1,3,2);

% 3D Plotting (Spherical to Cartesian mapping)
% AF_normalized = AF / max(AF(:));
% AF_db = 20 * log10(AF_normalized);
% AF_db(AF_db < -30) = -30; 
AF_db = 20 * log10(AF);
AF_db(AF_db < 0) = 0;  

% Map Elevation (Phi) to Radius. 
% Zenith (90 degrees) becomes the center (R = 0).
% Horizon (0 degrees) becomes the outer edge (R = 90).
% R_polar = 90 - rad2deg(Phi); 
R_polar = 90 - rad2deg(Phi); 

% Map Azimuth (Theta) to the circular angle
X_polar = R_polar .* cos(Theta);
Y_polar = R_polar .* sin(Theta);

% Plot the flat heatmap
surf(X_polar, Y_polar, AF_db, 'EdgeColor', 'none'); 
colormap('jet'); 
colorbar;
hold on;

% Draw Custom Polar Grid
% Hovers lines slightly above the plot so they don't get buried
Z_grid = max(AF_db(:)) + 1;

% Draw Elevation Rings (e.g., at 0, 30, and 60 degrees elevation)
theta_circle = linspace(0, 2*pi, 100);
elevations_to_mark = [0, 30, 60]; 
for el = elevations_to_mark
    r = 90 - el;    % Convert elevation to radius
    plot3(r*cos(theta_circle), r*sin(theta_circle), ...
        Z_grid*ones(size(theta_circle)), 'w--', 'LineWidth', 1, ... 
        'Color', [1 1 1 0.5]);
    text(r, 0, Z_grid, sprintf(' %d\\circ El', el), 'Color', 'white', ...
        'VerticalAlignment', 'bottom');
end

% 2. Draw Azimuth Spokes (every 45 degrees)
azimuths_to_mark = 0:45:315;
for az = azimuths_to_mark
    az_rad = deg2rad(az);
    plot3([0, 90*cos(az_rad)], [0, 90*sin(az_rad)], [Z_grid, Z_grid], ...
        'w--', 'LineWidth', 1, 'Color', [1 1 1 0.5]);
    text(98*cos(az_rad), 98*sin(az_rad), Z_grid, ...
        sprintf('%d\\circ Az', az), 'Color', 'white', ...
        'HorizontalAlignment', 'center');
end

axis equal off; 
t = title({['Polar Heatmap: Gain (Peak: ' num2str(max(AF_db(:)),'%.2f') ...
    'db)'], ['Az = ' num2str(rad2deg(theta_steer)) ', El = ' ...
    num2str(rad2deg(phi_steer))]});
t.Position = [0, 125, 0]; 
hold off;
view(2);

subplot(1,3,3);


% Convert spherical pattern to Cartesian for the surf() plot
% X_plot = AF_normalized .* cos(Phi) .* cos(Theta);
% Y_plot = AF_normalized .* cos(Phi) .* sin(Theta);
% Z_plot = AF_normalized .* sin(Phi);

X_plot = AF_db .* cos(Phi) .* cos(Theta);
Y_plot = AF_db .* cos(Phi) .* sin(Theta);
Z_plot = AF_db .* sin(Phi);

surf(X_plot, Y_plot, Z_plot, AF_db, 'EdgeColor', 'none');
colormap('jet');
axis tight;
%daspect([1 1 10])
colorbar;
title('Beam of Array');
xlabel('X-axis'); 
ylabel('Y-axis'); 
zlabel('Z-axis');
