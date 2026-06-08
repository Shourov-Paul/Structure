clear all
clc

% =========================================================================
% 1. Generate Lambda Array (Based on your provided parameters)
% =========================================================================
del_lambda = 100 * 1e-9; % 100 nm step
lambda = (1000 * 1e-9 : del_lambda : 12000 * 1e-9)'; % Wavelength in meters (1 um to 12 um)
N = length(lambda);

lambda_um = lambda * 1e6; % Convert to micrometers for plotting

% =========================================================================
% 2. Data Loading (A_eff only, 3 columns for R = 2, 2.25, 2.5 um)
% =========================================================================
csv_path = fullfile(fileparts(mfilename('fullpath')), 'Aeff_data_RV.csv');

if isfile(csv_path)
    if exist('readmatrix', 'file')
        A_eff_m2 = readmatrix(csv_path);
    else
        imported_data = importdata(csv_path);
        if isstruct(imported_data)
            A_eff_m2 = imported_data.data;
        else
            A_eff_m2 = imported_data;
        end
    end
    
    % Convert Aeff from m^2 (values around 10^-11) to um^2 for plotting
    A_eff_um2 = A_eff_m2 * 1e12;
    
    % Data validation: check if lengths match
    if size(A_eff_m2, 1) ~= N
        warning('Data length mismatch! Lambda has %d points, but Aeff data has %d points. Adjusting to match...', N, size(A_eff_m2, 1));
        min_N = min(N, size(A_eff_m2, 1));
        lambda = lambda(1:min_N);
        lambda_um = lambda_um(1:min_N);
        A_eff_m2 = A_eff_m2(1:min_N, :);
        A_eff_um2 = A_eff_um2(1:min_N, :);
    end
else
    error('Could not find Aeff_data_RV.csv. Please ensure it is in the same folder as this script.');
end

% =========================================================================
% 3. Calculate Nonlinear Parameter (Gamma) for each R value
% =========================================================================
% Using the reference n2 value from the paper:
n2_base = 7.33e-18; 
n2 = n2_base / 3; % [m^2/W]

% Calculate Gamma formula (for all 3 columns)
gamma = (2 * pi * n2) ./ (lambda .* A_eff_m2); % [m^-1 W^-1]

% =========================================================================
% 4. Plotting (Dual Y-Axis for R = 2, 2.25, 2.5 um)
% =========================================================================
figure(1)
clf;
hold on;

% Define beautiful, distinct colors for each R value
color_R20  = [0.0, 0.4, 0.8];  % Cobalt Blue for R = 2 um
color_R225 = [0.0, 0.6, 0.5];  % Teal for R = 2.25 um
color_R25  = [0.8, 0.1, 0.2];  % Crimson Red for R = 2.5 um

% --- Left Y-Axis: Gamma ---
yyaxis left
h_gamma1 = plot(lambda_um, gamma(:, 1), 'LineStyle', '--', 'Color', color_R20, 'LineWidth', 2.5);
h_gamma2 = plot(lambda_um, gamma(:, 2), 'LineStyle', '--', 'Color', color_R225, 'LineWidth', 2.5);
h_gamma3 = plot(lambda_um, gamma(:, 3), 'LineStyle', '--', 'Color', color_R25, 'LineWidth', 2.5);
ylabel('\gamma [m^{-1}W^{-1}]', 'FontSize', 20, 'FontWeight', 'bold');

ax = gca;
ax.YAxis(1).Color = 'k'; % Force left axis text to black

% --- Right Y-Axis: Aeff ---
yyaxis right
h_Aeff1 = plot(lambda_um, A_eff_um2(:, 1), 'LineStyle', '-', 'Color', color_R20, 'LineWidth', 2.5);
h_Aeff2 = plot(lambda_um, A_eff_um2(:, 2), 'LineStyle', '-', 'Color', color_R225, 'LineWidth', 2.5);
h_Aeff3 = plot(lambda_um, A_eff_um2(:, 3), 'LineStyle', '-', 'Color', color_R25, 'LineWidth', 2.5);
ylabel('A_{eff} [\mum^2]', 'FontSize', 20, 'FontWeight', 'bold');

ax.YAxis(2).Color = 'k'; % Force right axis text to black

% --- General Formatting ---
xlabel('Wavelength [\mum]', 'FontSize', 20, 'FontWeight', 'bold');
xlim([min(lambda_um), max(lambda_um)]); 

% Customizing box and ticks
set(gca, 'FontSize', 16, 'LineWidth', 1.5);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
box on;
grid off; % Enable grid for better readability

% --- Legend ---
legend([h_Aeff1, h_Aeff2, h_Aeff3], ...
    {'A_{eff}, \gamma, R=2 \mum', 'A_{eff}, \gamma, R=2.25 \mum', 'A_{eff}, \gamma, R=2.5 \mum'}, ...
    'FontSize', 14, 'FontWeight', 'bold', 'Location', 'north');
legend('boxoff');

% --- Bottom Text Annotation ---
text(0.5, 0.6, 'H = 0.1 \mum', 'Units', 'normalized', 'FontSize', 16, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center');

hold off;

% Save the plot
output_image_path = fullfile(fileparts(mfilename('fullpath')), 'Gamma_Aeff_plot_RV.png');
saveas(gcf, output_image_path);
disp('Plot saved successfully.');