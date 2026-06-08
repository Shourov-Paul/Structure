clear all
clc

% =========================================================================
% 1. Load or Generate the Dataset
% =========================================================================
% In your actual work, you would load your FEM data here:
% data = readmatrix('mode_profile.csv');
% x = unique(data(:,1)); y = unique(data(:,2)); 
% Z = reshape(data(:,3), [length(y), length(x)]);

% FOR THIS EXAMPLE: Generating a dummy dataset mimicking a triangular mode
x = linspace(6, 24, 200); % x coordinates in um
y = linspace(6, 24, 200); % y coordinates in um
[X, Y] = meshgrid(x, y);

% Simulating a triangular-ish intensity profile at the center (15, 15)
R = sqrt((X - 15).^2 + (Y - 15).^2);
Angle = atan2(Y - 15, X - 15);
% Adding a 3-fold symmetry to mimic the triangular core
Z = exp(-(R.^2)/8) .* (1 + 0.15*cos(3*Angle)); 
Z(Z < 0) = 0; % Ensure no negative intensity
Z = Z / max(Z(:)); % Normalize max to 1

% =========================================================================
% 2. Plotting the 2D Mode Profile
% =========================================================================
figure(1)
clf;

% Use imagesc for high-performance 2D heatmaps
% Note: imagesc flips the Y-axis by default, so we set 'YDir' to 'normal'
h = imagesc(x, y, Z);
set(gca, 'YDir', 'normal');

% Apply the 'jet' colormap to match the reference image
colormap('jet');

% Add and format the colorbar
cb = colorbar;
ylabel(cb, 'Intensity (a.u.)', 'FontSize', 16, 'FontWeight', 'bold');
% Set colorbar limits if you want to match the 0 to 0.3 scale from the paper
% caxis([0 0.3]); 

% =========================================================================
% 3. Formatting the Axes (matching your preferred style)
% =========================================================================
xlabel('x [\mum]', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('y [\mum]', 'FontSize', 20, 'FontWeight', 'bold');

% Ensure the aspect ratio is perfectly square (so the core isn't stretched)
axis image; 

% Adjust limits to match the specific bounding box of the core
xlim([9, 21]);
ylim([9, 21]);
xticks([9, 12, 15, 18, 21]);
yticks([9, 12, 15, 18, 21]);

set(gca, 'FontSize', 16, 'LineWidth', 1.5);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
box on;

% Add inner text label for the side length
text(15, 11, 'a = 8 \mum', 'FontSize', 18, 'FontWeight', 'bold', ...
    'FontAngle', 'italic', 'Color', 'white', 'HorizontalAlignment', 'center');

% Save the plot
output_image_path = fullfile(fileparts(mfilename('fullpath')), 'Mode_Profile.png');
saveas(gcf, output_image_path);
disp('Mode profile plot saved successfully.');