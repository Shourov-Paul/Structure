clear all
clc

N_betas = 10;
N_pump = 31;
c = 3e-7;
del_lambda = 100*1e-12;
lambda = (1000*1e-12:del_lambda:12000*1e-12);
N = length(lambda);

csv_path = fullfile(fileparts(mfilename('fullpath')), 'All_n_data.csv');
imported_data = importdata(csv_path);

if isstruct(imported_data)
    n_eff_data = imported_data.data;
else
    n_eff_data = imported_data;
end

% Compute dispersion for each of the 3 columns
D_all = zeros(N, 3);
for i = 1:3
    n_eff_col = n_eff_data(1:N, i);
    [betas_col, D_col] = dispersion(lambda, n_eff_col, N, del_lambda, N_pump);
    D_all(:, i) = D_col;
    
    if i == 3
        betas = betas_col;
        D = D_col;
    end
end

lambda = lambda * 1e3; % Convert to nm

% Set up figure data
plot_lambda = (lambda(2:N-1) / 1e-6)'; % Wavelength in um (ensuring column vector)
D_cropped = D_all(2:N-1, :);          % Matching rows for dispersion

% -------------------------------------------------------------------------
% TERMINAL DISPLAY: Dispersion at Wavelength = 4 um
% -------------------------------------------------------------------------
target_wavelength = 4.0; 
radius_labels = {'R = 2.00 um', 'R = 2.25 um', 'R = 2.50 um'};

fprintf('\n=======================================\n');
fprintf(' Dispersion Values at Wavelength = %.2f um\n', target_wavelength);
fprintf('=======================================\n');

for i = 1:3
    % Interpolate to find the exact value at 4 um
    D_at_4um = interp1(plot_lambda, D_cropped(:, i), target_wavelength, 'linear');
    % Changed %.4f to %.2f to display only two decimal places
    fprintf('%s : D = %.2f ps/nm/km\n', radius_labels{i}, D_at_4um);
end
fprintf('=======================================\n\n');
% -------------------------------------------------------------------------

figure(1)
clf;
hold on;

% 1. Plot the 3 dispersion lines
h1 = plot(plot_lambda, D_cropped(:, 1), 'k-',  'linewidth', 2.5); % Black solid
h2 = plot(plot_lambda, D_cropped(:, 2), 'r--', 'linewidth', 2.5); % Red dashed
h3 = plot(plot_lambda, D_cropped(:, 3), 'b-.', 'linewidth', 2.5); % Blue dash-dot

% 2. Plot reference axes lines
plot([1.5 12], [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'linewidth', 1.5);
plot([4 4], [-500 200], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'linewidth', 1.5);

xlabel('Wavelength [\mum]', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('D [ps/nm/km]', 'FontSize', 20, 'FontWeight', 'bold');
xlim([2.5, 12]);
ylim([-100, 50]);
set(gca, 'FontSize', 16, 'LineWidth', 1.5);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
box on;
grid off;

% 3. Legend and Text
legend([h1, h2, h3], 'R = 2 \mum', 'R = 2.25 \mum', 'R = 2.5 \mum', ...
    'FontSize', 18, 'FontWeight', 'bold', 'FontAngle', 'italic', 'Location', 'southeast');
legend('boxoff');

text(2.6, 40, 'H = 0.05 \mum', 'FontSize', 18, 'FontWeight', 'bold', 'FontAngle', 'italic');
hold off;

% Save the plot
output_image_path = fullfile(fileparts(mfilename('fullpath')), 'test_plot.png');
saveas(gcf, output_image_path);
disp('Plot saved successfully.');