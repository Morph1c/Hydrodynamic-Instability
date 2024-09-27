% Define constants
n = 0;   % Fixed value for n
delta_k = 0.1;  % Step size for k
k_max = 8;  % Maximum value of k
find_initial_guess = 1;  % Start with finding k0 and r0r

% Prompt for base flow selection
disp(' ');
disp(' Available base flows: ');
disp('    1: Lamb-Oseen vortex');
disp('    2: Q-vortex');
disp('	3: Carton & Mc Williams with gaussian axial jet')
disp('	4: Francis Turbine velocity profile')

baseflow_parameters = [];
baseflow_parameters(1) = ask_user('-> base flow', 1);
switch baseflow_parameters(1)
   case 2 % Q-vortex
    baseflow_parameters(2) = ask_user('Swirling q factor', 2);
   case 3 % CMW
    baseflow_parameters(2) = ask_user('alpha', 4);
	baseflow_parameters(3) = ask_user('Wo ', 0.5);
   case 4 % Francis turbine profile
    baseflow_parameters(2) = ask_user('R1', 0.4664);
	baseflow_parameters(3) = ask_user('R2', 0.1305);
    baseflow_parameters(4) = ask_user('Omega0', 0.3176);
	baseflow_parameters(5) = ask_user('Omega1', -0.6288);
    baseflow_parameters(6) = ask_user('Omega2', 2.2545);
	baseflow_parameters(7) = ask_user('U0', 0.3069);
    baseflow_parameters(8) = ask_user('U1', 0.0105);
	baseflow_parameters(9) = ask_user('U2', -0.3188);
end

% Range of m values
m_values = -2:-1:-6;

% Initialize figure for plotting real and imaginary parts
figure;
hold on;

% Loop over m values, calculate sigma for each m, and plot real and imaginary parts
for m = m_values
    % Compute sigma and k_range for the current m
    [k_range, sigma_values] = calculate_sigma_for_m_criteria(m, baseflow_parameters, n, delta_k, k_max, find_initial_guess);
    
    % Plot real part of sigma
    plot(k_range, real(sigma_values), 'DisplayName', ['Real Part of \sigma, m = ', num2str(m)]);
    
    % Plot imaginary part of sigma on the same figure
    plot(k_range, imag(sigma_values), '--', 'DisplayName', ['Imaginary Part of \sigma, m = ', num2str(m)]);
end

% Add labels, title, and legend
xlabel('k');
ylabel('\sigma (Real and Imaginary Parts)');
title('\sigma vs k for different m values');
legend show;
grid on;
