clc
clear all
close all
format short

% Robust CSV file path loading
% csv_path = fullfile(fileparts(mfilename('fullpath')), 'All_n_Confloss_R2.5(um).csv');
csv_path = fullfile(fileparts(mfilename('fullpath')), 'All_n_Confloss_R2.5(um)_rearranged.csv');
c = readcell(csv_path);

% Extract headers and data rows
headers = c(1, :);
data = c(2:end, :);

% Set up figure
figure(1)
clf;
hold on;

% Pre-allocate handles for 3 plots instead of 4
handles = zeros(1, 3);
plot_colors = {'k-', 'r--', 'b-.'}; % Black, Red, Blue styles (removed Green)

del_lambda = 100*1e-12;
lambda = (1000*1e-12:del_lambda:12000*1e-12);
lambda_grid = (lambda * 1e9)';

% Store data for inset plot
all_lambda = cell(1, 3);
all_loss = cell(1, 3);

% Loop limited to 3 columns
for col = 1:3
    vals = data(:, col);
    
    % Find non-missing values
    mask = ~cellfun(@(x) any(ismissing(x)), vals);
    neff_col = cell2mat(vals(mask));
    lambda_col = lambda_grid(mask);
    
    % Calculate Confinement Loss (dB/cm)
    loss_col = -8.66 * (2*pi ./ lambda_col) .* imag(neff_col) * 1e4;
    
    % Store values
    all_lambda{col} = lambda_col;
    all_loss{col} = loss_col;
    
    % Plot
    handles(col) = plot(lambda_col, loss_col, plot_colors{col}, 'linewidth', 2.5);
end

% Draw a rectangle on the main plot to indicate the zoom-in region
rectangle('Position', [11.0, 12, 1.0, 15], 'LineStyle', ':', 'LineWidth', 1.5, 'EdgeColor', [0.5, 0.5, 0.5]);

% Plot labels and limits
xlabel('Wavelength [\mum]', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Confinement Loss [dB/cm]', 'FontSize', 20, 'FontWeight', 'bold');
xlim([8, 12]);
ylim([0, 35]); % Fits the curves beautifully
set(gca, 'FontSize', 16, 'LineWidth', 1.5);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
box on;
grid off;

% Add legend in southwest for only the first 3 columns
legend(handles, 'H = 0.05 \mum', 'H = 0.1 \mum', 'H = 0.15 \mum', ...
    'FontSize', 18, 'FontWeight', 'bold', 'FontAngle', 'italic', 'Location', 'southwest');
legend('boxoff');

% Text annotation at appropriate scale (Y = 32) in the empty top-left region
text(8.2, 32, 'R = 2.5 \mum', 'FontSize', 18, 'FontWeight', 'bold', 'FontAngle', 'italic');
hold off;

% --- Create Inset Plot ---
% Define position for inset axes: [left, bottom, width, height]
inset_ax = axes('Position', [0.26, 0.50, 0.32, 0.35]);
box on;
hold on;

% Plot same data on inset axes
for col = 1:3
    plot(all_lambda{col}, all_loss{col}, plot_colors{col}, 'linewidth', 2.0);
end

% Set limits and formatting for inset
xlim(inset_ax, [11.0, 12.05]);
ylim(inset_ax, [12, 27]);
set(inset_ax, 'FontSize', 12, 'LineWidth', 1.2, 'YTick', [12, 15, 18, 21, 24, 27]);
set(inset_ax, 'XMinorTick', 'on', 'YMinorTick', 'on');
title(inset_ax, 'Zoom-in (Top Portion)', 'FontSize', 12, 'FontWeight', 'bold');
hold off;

% Save the plot to Conf_Loss directory
output_image_path = fullfile(fileparts(mfilename('fullpath')), 'loss_r=2.5.png');
saveas(gcf, output_image_path);
disp('Done');