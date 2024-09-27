%
% Generate all possible growth rate
% plotting them, reproducing results in Delbende, Huerre paper
%

type = 1; % LSA

% baseflow params
baseflow_parameters(1) = 2; % q-vortex
baseflow_parameters(2) = 0.8; % delbende paper, q


% azimuthal and axial parameters
perturbation_parameters(1) = 10e4;

% should loop over this
perturbation_parameters(2) = -11; % m
perturbation_parameters(3) = 2.1; % k
perturbation_parameters(4) = 0; % no corolois
perturbation_parameters(5) = 0; % no brunt-vaisala
perturbation_parameters(6) = 1.; % \ni = k

% numerical algorithm
grid_parameters(1) = 130;
grid_parameters(2) = 2; % algebraic mapping
grid_parameters(3) = 100;
grid_parameters(4) = 0; % no angle direction

% Convergence criterion:
% -------------------------
  
numerical_parameters(1) = 1.e-5;
  
% Percentage of coeffs considered in the computation of residuals:
% ----------------------------------------------------------------

numerical_parameters(2) = 10;

% plot params (no plots)
plot_parameters(1) = 0;
  
plot_parameters(2) = 0;
  
plot_parameters(3) = 0;

plot_parameters(4) = 0;

  
files_parameters = [];

% Define the range for perturbation_parameters(3)
% Define the range for perturbation_parameters(3) (axial wavenumber, k)
k_values = 0:0.1:8;  % k values from 0 to 8 with a step of 0.1

% First range up to just before 3
%k_values1 = 0:0.05:2.95;  % Up to 2.95
% Include 3 only once in the second range
%k_values2 = 3:0.1:8;
%k_values = [k_values1, k_values2];


% Define the range for perturbation_parameters(2) (azimuthal wavenumber, m)
m_values = -1.*[4];  % Example range for m values (azimuthal wavenumber)

% Pre-allocate modes_ matrix to store results
modes_ = NaN(length(m_values), length(k_values));  % Initialize a 2D array for modes

% Loop over m values (azimuthal wavenumber)
for m_idx = 1:length(m_values)
    perturbation_parameters(2) = m_values(m_idx);  % Set azimuthal wavenumber (m)

    fprintf('Mode: %d\n', perturbation_parameters(2));
    
    % Loop over k values (axial wavenumber)
    for k_idx = 1:length(k_values)
        perturbation_parameters(3) = k_values(k_idx);  % Set axial wavenumber (k)
        
        % Execute the lisa function for the current combination of m and k
        [eigenvalues, r, u_r, u_theta, u_z, rho, pressure, axial_vorticity] = ...
            lisa(type, baseflow_parameters, perturbation_parameters, grid_parameters, ...
                 numerical_parameters, plot_parameters);
        
        % Extract the complex number with the highest positive imaginary part
        positive_imag_eigenvalues = eigenvalues(imag(eigenvalues) > 0);  % Filter positive imaginary parts
        
        if ~isempty(positive_imag_eigenvalues)
            % Find the eigenvalue with the maximum positive imaginary part
            [~, max_idx] = max(imag(positive_imag_eigenvalues));
            modes_(m_idx, k_idx) = positive_imag_eigenvalues(max_idx);  % Store the eigenvalue

        else
            % If no positive imaginary eigenvalues, store NaN
            modes_(m_idx, k_idx) = NaN;
        end
    end
end

% Create a figure
figure;

% First subplot: Imaginary part
subplot(1,2,1);
hold on;

% Loop through each m value and plot the corresponding imaginary part of modes_
for m_idx = 1:length(m_values)
    plot(k_values, imag(modes_(m_idx, :)), 'LineWidth', 1.5, ...
        'DisplayName', ['m = ', num2str(m_values(m_idx))]);
end

% Label the axes
xlabel('k values (Axial Wavenumber)');
ylabel('Imaginary Part of Modes (Eigenvalues)');
title('Imaginary Part of Modes vs k Values for Different m');

% Add a legend
legend('show');

% Grid for better visualization
grid on;
hold off;

% Second subplot: Real part
subplot(1,2,2);
hold on;

% Loop through each m value and plot the corresponding real part of modes_
for m_idx = 1:length(m_values)
    plot(k_values, real(modes_(m_idx, :)), 'LineWidth', 1.5, ...
        'DisplayName', ['m = ', num2str(m_values(m_idx))]);
end

% Label the axes
xlabel('k values (Axial Wavenumber)');
ylabel('Real Part of Modes (Eigenvalues)');
title('Real Part of Modes vs k Values for Different m');

% Add a legend
legend('show');

% Grid for better visualization
grid on;
hold off;