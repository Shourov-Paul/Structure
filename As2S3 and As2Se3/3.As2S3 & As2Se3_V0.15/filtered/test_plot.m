clear all
clc
lambdapump = 3000;
N_betas = 10;
N_pump = 6;
c = 3e-7;
del_lambda = 100*1e-12;
lambda = (1000*1e-12:del_lambda:12000*1e-12);
N = length(lambda);
csv_path = fullfile(fileparts(mfilename('fullpath')), 'all_n_H0.15.csv');
imported_data = importdata(csv_path);

if isstruct(imported_data)
    n_eff_data = imported_data.data;
else
    n_eff_data = imported_data;
end

% Compute dispersion for each of the 4 columns
D_all = zeros(N, 4);
for i = 1:4
    n_eff_col = n_eff_data(1:N, i);
    [betas_col, D_col] = dispersion(lambda, n_eff_col, N, del_lambda, N_pump);
    D_all(:, i) = D_col;
    if i == 4
        betas = betas_col;
        D = D_col;
    end
end
lambda = lambda * 1e3; % Convert to nm

% Set up figure
figure(1)
clf;
plot_lambda = lambda(2:N-1) / 1e-6; % Convert to um
hold on;

% 1. Plot the 4 dispersion lines
h1 = plot(plot_lambda, D_all(2:N-1, 1), 'k-',  'linewidth', 2.5); % Black solid for H = 0.05 um
h2 = plot(plot_lambda, D_all(2:N-1, 2), 'r--', 'linewidth', 2.5); % Red dashed for H = 0.1 um
h3 = plot(plot_lambda, D_all(2:N-1, 3), 'b-.', 'linewidth', 2.5); % Blue dash-dot for H = 0.15 um
h4 = plot(plot_lambda, D_all(2:N-1, 4), 'g--', 'linewidth', 2.5); % Green dashed for H = 0.2 um

% 2. Plot the zero line in Gray so it doesn't interfere with the red lines
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

% 3. EXPLICITLY pass the handles [h1, h2, h3, h4] to the legend in southeast
legend([h1, h2, h3, h4], 'R = 2 \mum', 'R = 2.25 \mum', 'R = 2.5 \mum', 'R = 2.75 \mum', ...
    'FontSize', 18, 'FontWeight', 'bold', 'FontAngle', 'italic', 'Location', 'southeast');
legend('boxoff');

text(2.6, 40, 'H = 0.15 \mum', 'FontSize', 18, 'FontWeight', 'bold', 'FontAngle', 'italic');
hold off;

% Save the plot as an image in the code directory
output_image_path = fullfile(fileparts(mfilename('fullpath')), 'test_plot.png');
saveas(gcf, output_image_path);

disp('Done');