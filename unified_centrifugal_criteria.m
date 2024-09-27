function [sigma, r0, k0, r0r] = unified_centrifugal_criteria(baseflow_parameters, m, k, n, find_initial_guess, r_guess)
    % Define symbolic variable r
    syms r_sym

    % Define base flow profiles for different types
    type_flow = baseflow_parameters(1);
    switch type_flow
        case 1   % Lamb-Oseen vortex
            Omega_sym = (1 - exp(-r_sym.^2)) ./ r_sym.^2;
            W_sym = 0;  % Axial velocity is zero in this case

        case 2   % Q-vortex
            q = baseflow_parameters(2); % Set q manually or using ask_user function
            Omega_sym = q * (1 - exp(-r_sym.^2)) ./ r_sym.^2;
            W_sym = exp(-r_sym.^2);

        case 3   % Carton & McWilliam vortex
            Omega_sym = exp(-r_sym.^2);
            W_sym = 0;  % Axial velocity is zero in this case
        case 4 % Francis turbine velocity profile
            R1 = baseflow_parameters(2);
	        R2 = baseflow_parameters(3);
            Omega0 = baseflow_parameters(4);
	        Omega1 = baseflow_parameters(5);
            Omega2 = baseflow_parameters(6);
	        U0 = baseflow_parameters(7);
            U1 = baseflow_parameters(8);
	        U2 = baseflow_parameters(9); 
            Omega_sym = Omega0 * r_sym + Omega1 * (R1^2 ./ r_sym) .* (1 - exp(-1 .* (r_sym.^2 ./ R1^2 ))) + ...
                + Omega2 * (R2^2 ./ r_sym) .* (1 - exp(-1 .* (r_sym .^2 ./ R2^2 )));
            W_sym = U0 + U1 * exp(-1 .* (r_sym .^2 ./ R1^2 )) + U2 * exp(-1 .* (r_sym .^2 ./ R2^2 )); 
    end

    % Compute symbolic derivatives of Omega and W
    dOmega_sym = diff(Omega_sym, r_sym);     % First derivative
    ddOmega_sym = diff(dOmega_sym, r_sym);   % Second derivative
    dW_sym = diff(W_sym, r_sym);             % First derivative of W
    ddW_sym = diff(dW_sym, r_sym);           % Second derivative of W

    % Convert symbolic expressions to numerical function handles
    Omega = matlabFunction(Omega_sym, 'Vars', r_sym);
    dOmega = matlabFunction(dOmega_sym, 'Vars', r_sym);
    %ddOmega = matlabFunction(ddOmega_sym, 'Vars', r_sym);
    W = matlabFunction(W_sym, 'Vars', r_sym);
    dW = matlabFunction(dW_sym, 'Vars', r_sym);
    %ddW = matlabFunction(ddW_sym, 'Vars', r_sym);

    % Switch between criteria
    % Billant & Gallaire, 2013
    %m = ask_user("Value of m[0]: ", 0);
    %k = ask_user("Value of k[0]: ", 0);
    %n = ask_user("Value of n[0]: ", 0);
    kappa = sqrt(m^2 + k^2); % total wavenumber
            
    % Compute alpha
    alpha = asin(m / kappa); 

    % Define Phi and related expressions
    f_sym = cos(alpha)^2 + sin(alpha)^2 * (1 / r_sym^2);
    Phi_sym = 2 * Omega_sym * (2 * Omega_sym + r_sym * dOmega_sym);
    Vorticity_sym = 2 * Omega_sym + r_sym * dOmega_sym; % check if it's correct
    UpperPhi_sym = cos(alpha) / f_sym * (Phi_sym * cos(alpha) - 2 * sin(alpha) * Omega_sym * (dW_sym / r_sym));
    Lambda_sym = cos(alpha) * W_sym + sin(alpha) * Omega_sym;


    % Compute derivatives symbolically
    dUpperPhi_sym = diff(UpperPhi_sym, r_sym);
    ddUpperPhi_sym = diff(dUpperPhi_sym, r_sym);
    dLambda_sym = diff(Lambda_sym, r_sym);
    ddLambda_sym = diff(dLambda_sym, r_sym);
    dVorticity_sym = diff(Vorticity_sym, r_sym);
            
    % Compute H
    H_sym = 1i * sin(alpha) / r_sym^2 * (r_sym * dVorticity_sym - UpperPhi_sym / Omega_sym) + 1i * cos(alpha) * r_sym * diff(dW_sym / r_sym, r_sym);

    % Convert to numerical function handles
    UpperPhi = matlabFunction(UpperPhi_sym, 'Vars', r_sym);
    dUpperPhi =  matlabFunction(dUpperPhi_sym, 'Vars', r_sym);
    ddUpperPhi =  matlabFunction(ddUpperPhi_sym, 'Vars', r_sym);

    Lambda = matlabFunction(Lambda_sym, 'Vars', r_sym);
    dLambda =  matlabFunction(dLambda_sym, 'Vars', r_sym);
    ddLambda =  matlabFunction(ddLambda_sym, 'Vars', r_sym);

    f = matlabFunction(f_sym, 'Vars', r_sym);
    H = matlabFunction(H_sym, 'Vars', r_sym);
    % Compute first k0 and r0r that are initial wavenumber from which start
    % and initial guess for netwon, once you have this value you move a
    % little bit to right (k0 - delta_k, k0 + delta_k) and use as initial
    % guess r0

    switch find_initial_guess

        % If we just want to find initial guess for Netwon algorithm
        case 1
            % Define Phi and its derivative
            Phi_tilde_sym = 2 * Omega_sym * dOmega_sym * r_sym / (dOmega_sym ^ 2 * r_sym ^ 2 + dW_sym ^ 2) * (dOmega * r_sym * Vorticity_sym + dW_sym^2);
            dPhi_tilde_sym = diff(Phi_tilde_sym, r_sym);
            dPhi_tilde = matlabFunction(dPhi_tilde_sym, 'Vars', r_sym);
            Phi_tilde = matlabFunction(Phi_tilde_sym, 'Vars', r_sym);

            % Initial guess for r0r
            %r0r_guess = 1;  % Starting guess for r0r
            %r0r = fzero(dPhi_tilde, r0r_guess);  % Solve for r0r
            r0r = fminbnd(Phi_tilde, 0, 2);  % Search between 0.001 and 2
            
            %disp("r0r is")
            %disp(r0r)

            % Compute k0 = -m * Omega'(r0r) / W'(r0r)
            k0 = -m * dOmega(r0r) / dW(r0r);

            % Null output
            sigma = 0;
            r0 = 0;
        
        % If we actually want to find sigma using Billant, Gallaire
        % 2013 criteria in the inviscid limit
        case 2
            % Define the function F(r), complex-valued
            F_complex = @(r) dLambda(r) - 1i * dUpperPhi(r) ./ (2 * kappa * sqrt(-UpperPhi(r)) );

            % Initial guess for r
            %r_guess = 1+1i;  % Choose an initial guess

            % Solve the complex system
            options = optimoptions('fsolve', 'Display', 'off');
            r0 = fsolve(F_complex, r_guess, options);

            % Compute sigma zero
            sigma_0 = -1i * kappa * Lambda(r0) + sqrt(-UpperPhi(r0));
            sigma = sigma_0 - (2 * n + 1) / (2 * kappa * sqrt(2 * f(r0))) * sqrt(ddUpperPhi(r0) - (dUpperPhi(r0))^2 / (2 * UpperPhi(r0)) ...
                + 2i * sqrt(-UpperPhi(r0)) * kappa * ddLambda(r0) ) ...
                - H(r0) / (2 * kappa * f(r0));
    end
end
