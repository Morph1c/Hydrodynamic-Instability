function [k_range, sigma_values] = calculate_sigma_for_m_criteria(m, baseflow_parameters, n, delta_k, k_max, find_initial_guess)

    % Define the range of k values
    k_range = 0:delta_k:k_max;
    
    % Initialize arrays to store results
    sigma_values = zeros(size(k_range));
    r0_values = zeros(size(k_range));
    
    % Step 1: Find k0 and r0r using the criteria
    [~, ~, k0, r0r] = unified_centrifugal_criteria(baseflow_parameters, m, 0, n, find_initial_guess, 0);
    
    % Store sigma and r0 for k0
    [sigma, r0] = unified_centrifugal_criteria(baseflow_parameters, m, k0, n, 2, r0r);
    [~, idx_k0] = min(abs(k_range - k0));
    sigma_values(idx_k0) = sigma;
    r0_values(idx_k0) = r0;
    
    % Step 2: Move left from k0 in steps of delta_k
    for idx = idx_k0-1:-1:1
        k = k_range(idx);
        if idx == idx_k0 - 1
            r_guess = r0r;
        else
            r_guess = r0_values(idx + 1);
        end
        try
            [sigma, r0] = unified_centrifugal_criteria(baseflow_parameters, m, k, n, 2, r_guess);
            sigma_values(idx) = sigma;
            r0_values(idx) = r0;
        catch
            sigma_values(idx) = NaN;
            r0_values(idx) = NaN;
        end
    end
    
    % Step 3: Move right from k0 in steps of delta_k
    for idx = idx_k0+1:length(k_range)
        k = k_range(idx);
        r_guess = r0_values(idx - 1);
        try
            [sigma, r0] = unified_centrifugal_criteria(baseflow_parameters, m, k, n, 2, r_guess);
            sigma_values(idx) = sigma;
            r0_values(idx) = r0;
        catch
            sigma_values(idx) = NaN;
            r0_values(idx) = NaN;
        end
    end

    % Step 1: Set sigma_values to 0 where the real part of sigma is zero
    sigma_values(real(sigma_values) <= 0) = 0;
end
