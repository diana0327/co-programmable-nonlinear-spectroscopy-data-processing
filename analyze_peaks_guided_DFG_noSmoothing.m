function DFG_results = analyze_peaks_guided_DFG_noSmoothing(Ycompend_subtracted, X_cm, powerstep, DFG_range, power785, first_DFG_peak)

set(groot, 'defaultAxesFontName', 'Helvetica', 'defaultAxesFontSize', 8);
% Map the global first_DFG_peak to local DFG_range
first_DFG_peak_in_range = first_DFG_peak - (DFG_range(1) - 1);

last_peak_idx = first_DFG_peak_in_range;  % <-- Now it's correct local index

% Extract DFG sub-region
y_DFG = Ycompend_subtracted(DFG_range, :);
x_DFG = X_cm(DFG_range)';

numSamplesDFG = size(y_DFG, 2);

% Define color gradients
colorNames = {'skyBlue', 'green', 'pink2', 'purp', 'darkblue', 'lightblue', 'grayblue', 'bluepurp'};
DFG_colors = createColorGradient(colorNames, numSamplesDFG);

% Create figures
figure_DFG_raw = figure; set(figure_DFG_raw, 'Units', 'centimeters', 'Position', [15, 10, 6, 5]);
figure_DFG_fit = figure; set(figure_DFG_fit, 'Units', 'centimeters', 'Position', [15, 10, 6, 5]);

% Parameters
peak_window_pts = 3;          % Search up to 50 pts to the right
min_prominence = 1e-5;       % Minimum peak height
fit_half_window_pts = 30;      % Points to left and right for fitting

% Initialize output
DFG_results = struct('x', {}, 'y', {}, 'xpeak', {}, 'ypeak', {}, 'xFit', {}, 'yFit', {}, 'integral', {});

% Initialize starting point
% last_peak_idx = first_DFG_peak;  % first_DFG_peak is a data point index!

for i = 1:numSamplesDFG
    %% Process DFG (move to larger index --> larger wavelength)
    start_idx = last_peak_idx;

    % Plot raw signal
    figure(figure_DFG_raw); hold on;
    plot(x_DFG, y_DFG(:, i), 'Color', DFG_colors(i, :), 'LineWidth', 0.8);

    % Search in index window (move RIGHT)
    idx_window = max(1, start_idx-1) : min(length(x_DFG), start_idx + peak_window_pts);

%     if isempty(idx_window) || length(idx_window) < 5
%         warning('No points to search for sample %d. Skipping.', i);
%         continue;
%     end

    [pks, locs_rel] = findpeaks(y_DFG(idx_window, i), 'MinPeakProminence', min_prominence);
    locs_abs = idx_window(1) - 1 + locs_rel;  % Correct absolute index

    if isempty(pks)
        warning('No real peak found for sample %d.', i);
        continue;
    end

    % Take the strongest peak
    [~, strongest_idx] = max(pks);
    peak_idx_found = locs_abs(strongest_idx);

    last_peak_idx = peak_idx_found; % Update for next peak search

    % Fit around the found peak
    left = max(1, peak_idx_found - fit_half_window_pts);
    right = min(length(x_DFG), peak_idx_found + fit_half_window_pts);

    x_data = x_DFG(left:right);
    y_data = y_DFG(left:right, i) / powerstep(i, 2) / power785;

    [xmaxwave, maxiy_fit, xFit, yFit] = gaussianfit_stable(x_data, y_data);

    integral_fit = trapz(xFit, yFit);

    % Save results
    DFG_results{i} = struct('x', x_data, 'y', y_data, 'xpeak', xmaxwave, 'ypeak', maxiy_fit, 'xFit', xFit, 'yFit', yFit, 'integral', integral_fit);

    % Plot fit
    figure(figure_DFG_fit); hold on;
    plot(xFit, yFit, 'Color', DFG_colors(i, :), 'LineWidth', 0.8);
end
end
