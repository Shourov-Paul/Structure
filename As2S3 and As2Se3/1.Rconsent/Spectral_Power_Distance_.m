% Spectral_Power_Distance.m
% This script loads refractive index data from All_n_data.csv in a user-specified 
% folder, computes dispersion coefficients for each column (R = 2.00, 2.25, 2.50, 2.75 um),
% runs the supercontinuum generation simulation, and generates/saves:
% 1) Spectral Power vs Wavelength
% 2) Distance vs Wavelength (Spectral Density)
%
% Files are saved with names based on the input CSV column name in the script's directory.

% Preserve folder_path if running programmatically
if exist('folder_path', 'var')
    temp_folder_path = folder_path;
end
clearvars -except temp_folder_path;
if exist('temp_folder_path', 'var')
    folder_path = temp_folder_path;
    clear temp_folder_path;
end
clc;
close all;

% Add script directory to search path
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);

% ==================================================================================================================
% --- EDITABLE PARAMETERS (Simulation Configuration & Inputs) ------------------------------------------------------
% ==================================================================================================================

% Manually specified constant Aeff values for each structure (in m^2):
% Values corresponding to H = 0.05, 0.1, 0.15 um at R = 2.50 um and wavelength = 4 um
aeff_H05 = 1.36389525850693E-11; % Aeff for H = 0.05 um (Column 1)
aeff_H10 = 1.3682739349945E-11; % Aeff for H = 0.10 um (Column 2)
aeff_H15 = 1.37428672469036E-11; % Aeff for H = 0.15 um (Column 3)

% Primary parameters
lambdapump = 4100;           % pump wavelength [nm]
power = 3000;                % peak power of input [W]
FWHM = 50e-3;                % pulse width in FWHM [ps]

%% Dispersion Engineering Parameters
N_betas = 10;                % number of higher order dispersion terms
c_light = 3e-7;              % Unit of velocity (c) of light is km/ps
del_lambda = 100*1e-12;      % Wavelength step size in km
lambda = (1000*1e-12:del_lambda:12000*1e-12); % Wavelength array in km

% Nonlinearity and waveguide parameters
n2 = 1.1e-17;                % nonlinear refractive index [m^2/W] for GeAsSe
loss = 65;                   % loss [dB/m]
beta_TPA = 0;                % TPA coefficient [m/W]
chirp = 0;
mshape = 0;
dlength = 10e-3;             % waveguide length [m] (10 mm)

% Raman response parameters
fr = 0.013;                  % fractional Raman contribution for GeAsSe
tau1 = 21.3;                 % duration [ps]
tau2 = 195;                  % duration [ps]

% ==================================================================================================================
% --- END OF EDITABLE PARAMETERS -----------------------------------------------------------------------------------
% ==================================================================================================================

% Simulation grids setup
nt = 2^13;                   % number of grid points (FFT)
T = 20;                      % Time window [ps]
dt = T/nt;                   % time step
c_shock = 3e8*1e9/1e12;      % speed of light in nm/ps for shock calculations (300,000 nm/ps)
f0 = c_shock/lambdapump;     % pump frequency [THz]
w0 = 2*pi*f0;                % angular pump frequency
t = linspace(-T/2, T/2, nt); % time grid
V = 2*pi*(-nt/2:nt/2-1)'/(nt*dt); % frequency grid

% Raman response function
RT = ((tau1^2+tau2^2)/(tau1*tau2^2)).*sin(t/tau1).*exp(-t/tau2);
RT(t<0) = 0;                 % heaviside step function
RT = RT/trapz(t,RT);         % normalise RT to unit integral
Tr = 0;
nplot = 240;                 % number of length steps to save field at

% Locate and import All_n_dataf_R2.5(um).csv in the script's directory
csv_path = fullfile(fileparts(mfilename('fullpath')), 'All_n_dataf_R2.5(um).csv');
fprintf('Loading index data from: %s\n', csv_path);
imported_data = importdata(csv_path);

if isstruct(imported_data)
    n_eff_data = imported_data.data;
else
    n_eff_data = imported_data;
end

% Column mapping setup
column_names = {'H0.05', 'H0.1', 'H0.15'};
aeff_values = [aeff_H05, aeff_H10, aeff_H15];

N = length(lambda);

% Initialize summary structure array to store results
summary_data = struct('name', {}, 'aeff', {}, 'npump', {}, 'pump_wl', {}, 'neff', {}, 'dispersion', {}, 'N', {}, 'Lfiss', {});

% Loop through each column of n_eff data
for col_idx = 1:3
    current_name = column_names{col_idx};
    Aeff = aeff_values(col_idx);
    
    fprintf('\n--------------------------------------------------\n');
    fprintf('Processing structure (%d/3): %s\n', col_idx, current_name);
    fprintf('--------------------------------------------------\n');
    fprintf('Using manually specified constant Aeff: %.6e m^2\n', Aeff);
    
    n_eff_col = n_eff_data(1:N, col_idx);
    
    % Get dynamic index for the pump wavelength in the grid (using lambda in km)
    pump_wl_um = lambdapump / 1000;
    % Wavelength list in lambda (km) converted to um for comparison
    lambda_um = lambda * 1e9; 
    [~, N_pump] = min(abs(lambda_um - pump_wl_um));
    fprintf('--> Dynamic N_pump calculated: %d (corresponds to wavelength = %.2f um in lambda grid)\n', N_pump, lambda_um(N_pump));
    
    % Calculate dispersion and betas using dispersion.m
    [betas, D] = dispersion(lambda, n_eff_col, N, del_lambda, N_pump);
    betas = betas(1:N_betas);
    
    % Calculate nonlinear coefficient gamma
    gamma = 2*pi*n2/(1e-9*lambdapump*Aeff) + 1i*beta_TPA/(2*Aeff);
    
    % Input field setup
    T0 = FWHM/(2*acosh(sqrt(2)));
    if mshape == 0
        A_field = sqrt(power)*sech(t/T0).*exp(-1i*chirp*t.^2/(2*T0.^2));
    else
        A_field = exp(-0.5*(1+1i*chirp).*(t/T0).^(2*mshape));
    end
    
    % Run simulation using gnlse.m
    fprintf('Running GNLSE propagation...\n');
    [Z, AT, AW, W] = gnlse(t, V, A_field, w0, gamma, betas, loss, fr, Tr, dlength, RT, nplot);
    
    % Process output
    IW = abs(AW).^2; 
    mIW = max(max(IW));
    lIW = 10*log10(IW/mIW);
    
    % Soliton order and fission length calculation
    LD = T0^2/abs(betas(2)); % dispersion length [m]
    LNL = 1/(abs(gamma)*power); % non-linear length [m]
    soliton_N = sqrt(LD/LNL); % soliton order
    Lfiss_val = LD/soliton_N; % soliton fission length [m]
    
    IT = abs(AT).^2;
    mIT = max(max(IT));
    lIT = 10*log10(IT/mIT);
    
    WL = 2*pi*c_shock./W;
    iis = (WL > 300 & WL < 20000);
    
    % Safe base name for filenames
    safe_name = strrep(current_name, '.', '_');
    
    % 1. Plot Spectral Power vs Wavelength
    fig1 = figure('Name', 'Spectral Power vs Wavelength', 'NumberTitle', 'off');
    plot(WL(iis)/1000, lIW(240, iis), '-b', 'linewidth', 3);
    xlabel('Wavelength [\mum]', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Spectral Power [dB]', 'FontSize', 16, 'FontWeight', 'bold');
    set(gca, 'FontSize', 16, 'LineWidth', 1.5);
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
    xlim([0, 20]);
    ylim([-70, 0]);
    box on;
    grid on;
    
    % Save Plot 1
    output_img1 = fullfile(script_dir, [safe_name, '_Spectral_Power.png']);
    saveas(fig1, output_img1);
    close(fig1);
    
    % 2. Plot Distance vs Wavelength (Spectral Density evolution)
    fig2 = figure('Name', 'Distance vs Wavelength', 'NumberTitle', 'off');
    pcolor(WL(iis)/1000, Z*1000, lIW(:, iis));
    set(gca, 'CLim', [-40 0]);
    xlim([1.1, 20]);
    shading interp;
    colormap(jet.^3);
    xlabel('Wavelength [\mum]', 'FontSize', 26, 'FontWeight', 'bold');
    ylabel('Distance [mm]', 'FontSize', 26, 'FontWeight', 'bold');
    set(gca, 'FontSize', 22, 'LineWidth', 1.5);
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
    box on;
    view(2);
    
    % Save Plot 2
    output_img2 = fullfile(script_dir, [safe_name, '_Distance_vs_Wavelength.png']);
    saveas(fig2, output_img2);
    close(fig2);
    
    % 3. Plot Delay vs Distance (Temporal Density evolution)
    fig3 = figure('Name', 'Delay vs Distance', 'NumberTitle', 'off');
    pcolor(t, Z*1000, lIT);
    set(gca, 'CLim', [-40 0]);
    xlim([-0.5, 2.2]);
    shading interp;
    colormap(jet.^3);
    xlabel('Delay [ps]', 'FontSize', 26, 'FontWeight', 'bold');
    ylabel('Distance [mm]', 'FontSize', 26, 'FontWeight', 'bold');
    set(gca, 'FontSize', 22, 'LineWidth', 1.5);
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
    box on;
    view(2);
    
    % Save Plot 3
    output_img3 = fullfile(script_dir, [safe_name, '_Delay_vs_Distance.png']);
    saveas(fig3, output_img3);
    close(fig3);
    
    fprintf('Successfully generated and saved plots:\n');
    fprintf('  - %s\n', output_img1);
    fprintf('  - %s\n', output_img2);
    fprintf('  - %s\n', output_img3);
    
    % Store values in summary structure
    summary_data(col_idx).name = current_name;
    summary_data(col_idx).aeff = Aeff;
    summary_data(col_idx).npump = N_pump;
    summary_data(col_idx).pump_wl = lambda_um(N_pump);
    summary_data(col_idx).neff = n_eff_col(N_pump);
    summary_data(col_idx).dispersion = D(N_pump);
    summary_data(col_idx).N = soliton_N;
    summary_data(col_idx).Lfiss = Lfiss_val;
end

disp('========================================================================================================');
disp('                                   SUMMARY OF EFFECTIVE VALUE DETAILS');
disp('========================================================================================================');
fprintf('%-10s | %-12s | %-8s | %-12s | %-10s | %-12s | %-12s | %-12s\n', ...
    'Structure', 'Aeff (m^2)', 'N_pump', 'Pump WL (um)', 'n_eff', 'D (ps/nm/km)', 'N', 'Lfiss (m)');
disp('--------------------------------------------------------------------------------------------------------');
for k = 1:3
    fprintf('%-10s | %-12.4e | %-8d | %-12.2f | %-10.6f | %-12.4f | %-12.4f | %-12.4e\n', ...
        summary_data(k).name, ...
        summary_data(k).aeff, ...
        summary_data(k).npump, ...
        summary_data(k).pump_wl, ...
        summary_data(k).neff, ...
        summary_data(k).dispersion, ...
        summary_data(k).N, ...
        summary_data(k).Lfiss);
end
disp('========================================================================================================');

disp('==================================================');
disp('All columns of All_n_dataf_R2.5(um).csv processed successfully.');
disp('==================================================');
