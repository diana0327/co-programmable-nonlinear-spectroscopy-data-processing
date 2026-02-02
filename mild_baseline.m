function [Base, Corrected_Spectrum] = baseline_mild(Spectrum,n)
    % Input validation
    if isempty(Spectrum) || length(Spectrum) < 3
        error('Spectrum must have at least 3 elements.');
    end

    l = length(Spectrum);
    lp = ceil(0.5 * l);

    % Extend the spectrum with padding
    initial_Spectrum = [ones(lp, 1) * Spectrum(1); Spectrum(:); ones(lp, 1) * Spectrum(end)];

    % Initialize variables
    S = initial_Spectrum;
%     n = 33; % Start with a larger initial value for milder removal
    max_iterations = 20; % Limit the maximum number of iterations
    area_tolerance = 1e-3; % Tolerance for stopping based on area change
    flag1 = 0;

    A = []; % Array to store area under the curve
    i_min = 1; % Default value for i_min to avoid errors

    for iter = 1:max_iterations
        n = n + 2; % Increment n
        i = (n - 1) / 2;

        % Perform peak stripping
        [Baseline, ~] = peak_stripping(S, n);

        A(i) = trapz(S - Baseline);
        Stripped_Spectrum{i} = Baseline;
        S = Baseline;

        % Check for stopping condition based on area change
        if i > 3
            if A(i-1) < A(i-2) && A(i-1) < A(i)
                i_min = i-1;
                flag1 = 1;
            elseif abs(A(i) - A(i-1)) / A(i-1) < area_tolerance
                i_min = i; % Stop if area change is small
                flag1 = 1;
            end
        end

        if flag1 == 1
            break;
        end
    end

    % If no stopping condition was met, use the last iteration
    if flag1 == 0
        i_min = i; % Use the last index if no minimum was found
    end

    % Extract the baseline
    Base = Stripped_Spectrum{i_min};

    % Validate sizes before subtraction
    if length(Base) ~= length(initial_Spectrum)
        error('Baseline length does not match initial_Spectrum length.');
    end

    % Correct the spectrum
    Corrected_Spectrum = initial_Spectrum - Base;

    % Extract the central portion
    Corrected_Spectrum = Corrected_Spectrum(lp+1:lp+l);
    Base = Base(lp+1:lp+l);
%end compared to function [Base, Corrected_Spectrum] = baseline(Spectrum)
    % Input validation
    if isempty(Spectrum) || length(Spectrum) < 3
        error('Spectrum must have at least 3 elements.');
    end

    l = length(Spectrum);
    lp = ceil(0.5 * l);

    % Extend the spectrum with padding
    initial_Spectrum = [ones(lp, 1) * Spectrum(1); Spectrum(:); ones(lp, 1) * Spectrum(end)];

    % Initialize variables
    S = initial_Spectrum;
    n = 1;
    flag1 = 0;

    while flag1 == 0
        n = n + 2;  % Increment n
        i = (n - 1) / 2;

        % Perform peak stripping
        [Baseline, stripping] = peak_stripping(S, n);

        A(i) = trapz(S - Baseline);
        Stripped_Spectrum{i} = Baseline;
        S = Baseline;

        % Check for minimum condition
        if i > 3
            if A(i-1) < A(i-2) && A(i-1) < A(i)
                i_min = i-1;
                flag1 = 1;
            end
        end
    end

    % Extract the baseline
    Base = Stripped_Spectrum{i_min};

    % Validate sizes before subtraction
    if length(Base) ~= length(initial_Spectrum)
        error('Baseline length does not match initial_Spectrum length.');
    end

    % Correct the spectrum
    Corrected_Spectrum = initial_Spectrum - Base;

    % Extract the central portion
    Corrected_Spectrum = Corrected_Spectrum(lp+1:lp+l);
    Base = Base(lp+1:lp+l);
end