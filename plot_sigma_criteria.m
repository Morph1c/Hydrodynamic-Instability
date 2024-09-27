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
disp('    3: Carton & Mc Williams with gaussian axial jet');
disp('    4: Francis Turbine velocity profile');

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
m_values = -2:-1:-4;

% Preallocate a matrix to store sigma_values for all m values
% Rows are different k, columns are different m values
num_m = length(m_values);
k_range = 0:delta_k:k_max;
num_k = length(k_range);
sigma_matrix = zeros(num_k, num_m);  % Tensor for storing sigma for all m values

% Loop over m values and calculate sigma for each m
for i = 1:num_m
    m = m_values(i);
    
    % Compute sigma and k_range for the current m
    [~, sigma_values] = calculate_sigma_for_m_criteria(m, baseflow_parameters, n, delta_k, k_max, find_initial_guess);
    
    % Store sigma values for the current m in the matrix
    sigma_matrix(:, i) = sigma_values;
end

% Create figure for subplots
figure;

% Create first subplot for real part of sigma
subplot(2, 1, 1);  % 2 rows, 1 column, first subplot
hold on;
for i = 1:num_m
    % Plot real part of sigma for each m
    plot(k_range, real(sigma_matrix(:, i)), 'DisplayName', ['Real Part of \sigma, m = ', num2str(m_values(i))]);
end

% Add labels, title, and legend for real part
xlabel('k');
ylabel('Real Part of \sigma');
title('Real Part of \sigma vs k for different m values');
legend show;
grid on;

% Create second subplot for imaginary part of sigma
subplot(2, 1, 2);  % 2 rows, 1 column, second subplot
hold on;
for i = 1:num_m
    % Plot imaginary part of sigma for each m
    plot(k_range, imag(sigma_matrix(:, i)), '--', 'DisplayName', ['Imaginary Part of \sigma, m = ', num2str(m_values(i))]);
end

% Add labels, title, and legend for imaginary part
xlabel('k');
ylabel('Imaginary Part of \sigma');
title('Imaginary Part of \sigma vs k for different m values');
legend show;
grid on;

