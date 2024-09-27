%
% Analyze the bifurcation hypothesis for m=-4 
% compared to Delbdende paper
%

type = 1; % LSA

% baseflow params
baseflow_parameters(1) = 2; % q-vortex
baseflow_parameters(2) = 0.8; % delbende paper, q


% azimuthal and axial parameters
perturbation_parameters(1) = 10e4;

% should loop over this
perturbation_parameters(2) = -4; % m
%perturbation_parameters(3) = 2.1; % k
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
%k_values = 3.6;

% Initialize a cell array to store the eigenvalues with positive imaginary parts
modes_ = cell(1, length(k_values));  % Initialize a cell array for modes

% Loop over k values (axial wavenumber)
for k_idx = 1:length(k_values)
    perturbation_parameters(3) = k_values(k_idx);  % Set axial wavenumber (k)
        
    % Execute the lisa function for the current combination of m and k
    [eigenvalues, r, u_r, u_theta, u_z, rho, pressure, axial_vorticity] = ...
        lisa(type, baseflow_parameters, perturbation_parameters, grid_parameters, ...
                 numerical_parameters, plot_parameters);
        
    % Extract the eigenvalues with positive imaginary parts
    positive_imag_eigenvalues = eigenvalues(imag(eigenvalues) > 0);  % Filter positive imaginary parts
        
    if ~isempty(positive_imag_eigenvalues)
       % Store all eigenvalues with positive imaginary part for this k value
       modes_{k_idx} = positive_imag_eigenvalues;  % Store the entire array of positive eigenvalues
    else
       % If no positive imaginary eigenvalues, store NaN
       modes_{k_idx} = NaN;  % Store NaN if no positive modes found
    end
end


% Loop over k values and plot the eigenvalues with positive imaginary part
figure;  % Create a new figure

% First subplot: Imaginary part
subplot(1,2,1);

hold on; % Hold on to plot multiple points in the same figure
xlabel('k values (Axial Wavenumber)');
ylabel('Imaginary Part of Eigenvalues (Growth rate)');
title('Imaginary Part, continous line from unified criteria, n = 1, m = -4, Re = 10e4');

% Iterate over k values to plot each set of eigenvalues
for k_idx = 1:length(k_values)
    if ~isnan(modes_{k_idx})
        % Extract the eigenvalues for this k value
        positive_imag_eigenvalues = modes_{k_idx};
        
        % Plot real part vs imaginary part of the eigenvalues for this k
        plot(k_values(k_idx), imag(positive_imag_eigenvalues), 'o');
    end
end

% Just for now
plot(k_range, real(sigma_values), 'b-', 'LineWidth', 2);

%legend('show');  % Show legend with k values
hold off;

% Second subplot: Real part
subplot(1,2,2);
hold on;

% Iterate over k values to plot each set of eigenvalues
for k_idx = 1:length(k_values)
    if ~isnan(modes_{k_idx})
        % Extract the eigenvalues for this k value
        positive_imag_eigenvalues = modes_{k_idx};
        
        % Plot real part vs imaginary part of the eigenvalues for this k
        plot(k_values(k_idx), real(positive_imag_eigenvalues), 'o');
    end
end

% Just for now
plot(k_range, -1*imag(sigma_values), 'r--', 'LineWidth', 2);

% Label the axes
xlabel('k values (Axial Wavenumber)');
ylabel('Real Part of Eigenvalues (Frequency)');
title('Real Part, dashed line from unified criteria, n = 1, m = -4, Re = 10e4');


