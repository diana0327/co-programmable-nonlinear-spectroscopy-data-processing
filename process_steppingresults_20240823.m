%%% program for processing stepping CW mode results


clearvars -except powerprofile Powerstep  OPTIR; close all; clc
cd('C:\Users\zhxie\Pictures\images 20250605\chi3 all molecules\cnbpt chi3')
[file,path] = uigetfile('*.sif',...
   'Select One or More Files', ...
   'MultiSelect', 'on');
filenamespectrum = fullfile(path, file);
[Xnm,Y]=sireadkinetic(filenamespectrum);
% for i=1:7
% [Xnm,Y(:,i)]=sireadkinetic(filenamespectrum{i});
% end

% Y= reshape(data(:,2), [2048,length(data)/2048]); % Ym is 2048 x 100 matrix
% 
% % Extract the first row (unique entries, just the first 2048)
% Xnm = data( 1:2048,:); % x is a 1 x 2048 row vector
%% 


npnumber=3;
% [Xref,Yref]=sireadkinetic(filename);
% figure;hold on
% plot(Xnm(1001:end)+6.026,r(1001:end))
% plot(Xref(1001:end),Yref(1001:end,1))
% xlim([850 880])
% 
% figure;hold on
% plot(Xnm(1:1000)+4.4460,r(1:1000))
% plot(Xref(1:1000),Yref(1:1000))

% Y=y;

% X=double(convert_nm_to_cm(7.7968e+02,Xnm));
%% 

[file,path] = uigetfile('*.txt',...
   'Select One or More Files', ...
   'MultiSelect', 'on');
filename = fullfile(path, file);

Powerstep = readmatrix(filename);

% figure;hold on;
% for i=1:7
% powerstep = readmatrix(filename{i});
% plot(powerstep(:,2)),'DisplayName', filename{i})
% end
% legend

%% 
clear powerstep
powerstep=[Powerstep(1:end,1),Powerstep(1:end,2)*10^6*0.4];%W to uW
% powerstep=[Powerstep(1:5:end,1),Powerstep(1:5:end,2)*10^6];%W to uW

power785=1*10^3/3; %MW to UW
% powerprofile=p_power(900,1150);

figure;plot(powerstep(:,1),powerstep(:,2),'-*')
%% 

figure;plot(Y(:,100:200))
% figure;plot(Y)
%% 
region=1:size(Y,2)
% numwrong=[309,317,537,700,703,707,715,788,938,951,1032,1097,1246,1248,1490,1580,1757,1759]';
% numwrong=[246,339,352,487,644,1434]';
% numwrong=[45,258,1067,1690]';
% numwrong=[58,185,916,1437,1442,1443,1444]';
% numwrong=[46,64,132,184,198,221,423,602,610,616,618,624,726,748,790,808,832,859,872,948,972,989,999,1005,1011,1182,1251,1252,1276,1278,1307,1317,1331,1416,1417,1427,1428,1432,1437,1441,1490,1561,1591,1627,1645,1673,1680,1871,1907,1949,1979,1996]';
% numwrong=[702,703,707,715,788 939,951,969,1032, 1097 1246,1248,1490,1505,1580,1757,1972];
% numwrong=[537,1926,309,1839,252,387,317,359,708];
numwrong=[1916];
% numwrong=[1248,938,1032,1580,969,309,537,1097,1757,1759,1826,1839,1842,317,359,387,702,707,715,788,951,1246,1490,1971,708];
for i=1:length(numwrong);
[inten,location]=max(Y(numwrong(i),region));

Y(numwrong(i)-6:numwrong(i)+5,location)=Y(numwrong(i)-6:numwrong(i)+5,location-1);
% i
end
%% 

%%%%%%adjust the center
[M,peakAnti]=max(Y(199:554,1));
[M,peakS]=max(Y(1448:1933,1));
% Adjust indices to be absolute, not relative to the sliced ranges
% Initial peak positions (you can adjust these as needed)
peakAnti = peakAnti + 198;
peakS = peakS + 1447;
%% find the shift should be applied to the nm to cm conversion based on the peaks find on two sides. 

% Define the range for the random offset search
offset_min = -0.25;
offset_max = 5;

% Set a small tolerance for matching wavenumbers
tolerance = 1e-4;

% Set the number of random iterations
num_iterations = 10000;

% Set a small tolerance for matching wavenumbers
tolerance = 1e-4;

% Set the number of random iterations
num_iterations = 10000;

% Initialize variables to store the best offset and smallest difference
best_offset = NaN;
min_difference = Inf;

% Perform the random search
for i = 1:num_iterations
    % Generate a random offset within the specified range
    offset = offset_min + (offset_max - offset_min) * rand();
    
    % Convert wavelength in nm to cm and add the random offset
    lambcm = double(convert_nm_to_cm(7.80e+02, Xnm + offset));

    % Calculate the difference between the rounded values at peak indices
    difference = abs(-round(lambcm(peakAnti), 4) - round(lambcm(peakS), 4));

    % Check if this difference is the smallest found so far
    if difference < min_difference
        min_difference = difference;
        best_offset = offset;
        
        % If the difference is within the tolerance, we can stop searching
        if difference < tolerance
            break;
        end
    end
end
offset=best_offset
offset=0
% Use the best offset found
X_cm = double(convert_nm_to_cm(7.80e+02, Xnm + offset));
%% 
% figure;plot(X_cm,Y(:,170:200))
figure;plot(X_cm,Y)
%% no use unless useful
Ysignal=Y;
Y=Ysignal(:,20:end);
powerraw=powerstep;
powerstep=powerraw(20:end,:);
%% %% rebuild raman from the last one and the one closest to the last one: original one
clear i Ycompend Y_baseremoved scaling_factor Base powerstep_onlyvalid
close all
figure; plot(Y)
numColumns = size(Y, 2);  % Get the number of columns in Y

% Initialize matrices to store results
Y_baseremoved = zeros(size(Y));
Ycompend = zeros(size(Y));
scaling_factor = zeros(1, numColumns);  % Store scaling factors
% Y=Y(:, end:-1:1);
% Indices for valid columns
valid_indices = [];

% Indices for different regions
idx1 = 683:835;
% idx2 = 1003:1300;
idx2 = 1090:1151; %np11
% idx3 = 1553:1613;
% idx3 = 1576:1613; %NP18
idx3 = 1152:1411;
% idx4 = 1730:1810;
% idx4 = 1944:2027;
% idx4 = 1970:2027;
% idx4 = 954:962;
idx4 = 683:835;
compnum = 1; % Reference column number
n = 1;

% Thresholds for scaling factor
scaling_threshold_high = 1.9;%1.5;    % Upper limit
scaling_threshold_low = 0.3;     % Lower limit normally 0.7
Y_baseremoved=Y-mean(Y(1:100,1));
% Loop through each column of Y
for i = 1:760%numColumns
    % Perform baseline removal using the primary method
    [Base, Y_baseremoved(:,i)] = mild_baseline(double(Y(:,i)'),49);
    
    % Refine the Y data using the specified column as reference
    [Ycompend(:,i), scaling_factor(i)] = refine_y_up_multispecs(Y_baseremoved(:,i), Y_baseremoved(:, compnum), idx1, idx2, idx3, idx4);
    
    % Check if the scaling factor is within the acceptable range
    if scaling_factor(i) > scaling_threshold_high || scaling_factor(i) < scaling_threshold_low
        fprintf('Scaling factor %.2f for column %d is out of range. Applying mild baseline removal...\n', scaling_factor(i), i);
        
        % Apply a milder baseline removal
        [Base, Y_baseremoved(:,i)] = mild_baseline(double(Y(:,i)'),49);
        % Recompute the scaling factor
        [Ycompend(:,i), scaling_factor(i)] = refine_y_up_multispecs(Y_baseremoved(:,i), Y_baseremoved(:, compnum), idx1, idx2, idx3, idx4);
        
        % Check scaling factor again after mild baseline removal
        if scaling_factor(i) > scaling_threshold_high || scaling_factor(i) < scaling_threshold_low
            fprintf('Skipping column %d as scaling factor %.2f is still invalid after fallback.\n', i, scaling_factor(i));
            
            % Plot original and compensated spectra for debugging
            figure;
            subplot(2,1,1); hold on;
            plot(Y(:,compnum), 'DisplayName', 'Column Y (Reference)');
            plot(Y(:,i), 'DisplayName', sprintf('Column %d Y', i));
            title(sprintf('Column %d: Scaling Factor = %.2f', i, scaling_factor(i)));
            legend show;
            
            subplot(2,1,2); hold on;
            plot(Ycompend(:,compnum), 'DisplayName', 'Column Y (Reference)');
            plot(Ycompend(:,i), 'DisplayName', sprintf('Column %d Y', i));
            title('After Compensation');
            legend show;
            
            continue; % Skip this column if scaling factor is still invalid
        end
    end
    
    % Additional check: exclude spectra where power(:,2) < 1000
    if powerstep(i,2) < 2000%5000%14000%17560%17460
        fprintf('Excluding column %d as power(:,2) = %.2f is less than 14000.\n', i, powerstep(i,2));
        continue; % Skip this column
    end
    
    % If scaling factor is valid and power(:,2) >= 1000, add the index to valid_indices
    valid_indices(end + 1) = i;
end

% Modify powerstep variable according to valid columns
powerstep_onlyvalid = powerstep(valid_indices, :);

% powerstep_onlyvalid=powerstep;
% Ycompend=Y;

%% Plotting the results
range_plotted = 7:5:length(valid_indices) - 5;
figure;
subplot(1,2,1); plot(X_cm(1350:1510), Y(1350:1510, valid_indices(range_plotted)));
title('Stokes Before Raman Fluctuation Compensation');
subplot(1,2,2); plot(X_cm(1350:1510), Ycompend(1350:1510, valid_indices(range_plotted)));
title('Stokes After Raman Fluctuation Compensation');

figure;
subplot(1,2,1); plot(X_cm(325:440), Y(325:440, valid_indices(range_plotted)));
title('Anti-Stokes Before Raman Fluctuation Compensation');
subplot(1,2,2); plot(X_cm(325:440), Ycompend(325:440, valid_indices(range_plotted)));
title('Anti-Stokes After Raman Fluctuation Compensation');

figure;
subplot(2,1,1); plot(X_cm, Y(:, valid_indices));
title('Before Raman Fluctuation Compensation');
xlim([-1300 1700]); ylim([400 2200]);
xlim([0 1000]); ylim([550 550+310]);
subplot(2,1,2); plot(X_cm, Ycompend(:, valid_indices));
title('After Raman Fluctuation Compensation');
% xlim([-1300 1700]); ylim([-100 1800]);
xlim([0 1000]); ylim([-10 300]);

% figure;
% subplot(1,2,1);plot(X_cm,Y);title('before Raman fluct compen');ylim([550 580]);xlim([-170 200])
% subplot(1,2,2);plot(X_cm,Ycompend);title('after Raman fluct compen');ylim([-10 20]);xlim([-170 200])

%% %% rebuild raman from the last one and the one cloest to the last one
% Concatenate all indices into one vector
% Define your indices
idx5 = 1152:1411;
% all_indices = [idx1, idx2, idx3, idx4, idx5];
all_indices = [idx1, idx2, idx3, idx4];
% Assign weights
weight_idx345 = 3;  % For example, give idx5 twice the weight
weight_other = 1; % Weight for the other indices

% Calculate the mean differences for each region separately
mean_diff_other = mean(abs(Ycompend([idx1, idx2], :) - Ycompend([idx1, idx2], compnum)), 1);
mean_diff_idx345 = mean(abs(Ycompend([idx3, idx4,idx5], :) - Ycompend([idx3, idx4,idx5], compnum)), 1);

% Combine them with different weights
mean_diff = (weight_other * mean_diff_other + weight_idx345 * mean_diff_idx345) / (weight_other + weight_idx345);

% Initialize variables
found_valid_location = false;
% Iterate until a valid location is found
while ~found_valid_location
    % Find the index of the minimum mean difference
    [Mvalue, location] = min(mean_diff);
    
    % Check if the condition abs(compnum - location) > 70 is met and
    % the location is within valid_indices
    if abs(compnum - location) > 500 && ismember(location, valid_indices)
        % If conditions are met, mark as found and proceed
        found_valid_location = true;
    else
        % Exclude the current location by setting it to a high value
        mean_diff(location) = inf;
    end
end

% Display the results
disp(['Valid location found at index: ', num2str(location)]);
disp(['Difference from compnum: ', num2str(abs(compnum - location))]);
disp(['Minimum mean difference value: ', num2str(Mvalue)]);


%% compare to the last spectrum
% Example data (assuming Ycompend is your data matrix)
% Ycompend is an N x M matrix where each column represents a curve

% Define the columns for comparison
figure; hold on;plot(Ycompend(:,1),'DisplayName', 'fill the region in the section below');plot(Ycompend(:,location),'DisplayName', 'select the peak of blue curve');legend

% figure; hold on;plot(X_cm,Ycompend(:,2),'DisplayName', 'to be replaced');plot(X_cm,Ycompend(:,location),'DisplayName', 'parts osed for replaced');legend

% %% if compnum=148
% 
% end_column = size(Y, 2); % The column with the reference curve
% 
% % Define the regions to compare
% % regions = [320, 360; 1490, 1520]; %for NP13
% % regions = [390, 440; 1370, 1420]; %for NP19 
% regions = [400, 440; 1380, 1420]; %for Np20
% 
% % Initialize variables for storing indices
% start_indices = zeros(1, 2);
% end_indices = zeros(1, 2);
% 
% % Loop through the regions to find the biggest differences
% for i = 1:size(regions, 1)
%     % Extract the curves from Ycompend
%     curve1 = Ycompend(regions(i, 1):regions(i, 2), location);
%     curve2 = Ycompend(regions(i, 1):regions(i, 2), end_column);
% 
%     % Find the biggest difference in the current region
%     [~, start_indices(i), end_indices(i)] = find_biggest_difference(curve1, curve2, 15);
% end
% 
% % Build the RamanF curve
% RamanF = [...
%     Ycompend(1:start_indices(1)-1+regions(1, 1), end_column)', ...
%     Ycompend(start_indices(1)+regions(1, 1):end_indices(1)+regions(1, 1), location)', ...
%     Ycompend(end_indices(1)+1+regions(1, 1):start_indices(2)-1+regions(2, 1), end_column)', ...
%     Ycompend(start_indices(2)+regions(2, 1):end_indices(2)+regions(2, 1), location)', ...
%     Ycompend(end_indices(2)+1+regions(2, 1):end, end_column)' ...
% ];
% 
% % Plot the results
% figure;
% plot(RamanF);
% hold on;
% plot(Ycompend(:, end_column));
% legend('RamanF', 'Original Reference Curve');
% xlabel('Data Points');
% ylabel('Intensity');
% title('Comparison of RamanF and Original Reference Curve');
% hold off;


%% 

figure;plot(Ycompend(:,1:end))



%% compare to the first spectrum
% Example data (assuming Ycompend is your data matrix)
% Ycompend is an N x M matrix where each column represents a curve

% Define the columns for comparison

 % The column with the reference curve

% Define the regions to compare
% regions = [380, 430; 1370, 1400]; % the NP 13 
% regions = [385, 430; 1375, 1408]; % the NP 13 NP17
% regions = [385, 430; 1375, 1417]; % the NP 18
% regions = [394, 445; 1370, 1415];% the NP 19
% regions = [394, 445; 1381, 1417];% the NP 20 21 22
% regions = [533-15, 533+15; 1447-15, 1447+15];% the NP 20 21 22 ne results
peak1=531;peak2=1445;
regions = [peak1-15, peak1+15; peak2-15, peak2+15];% the NP 20 21 22 ne results
% regions = [498, 506; 1427, 1465];
% Initialize variables for storing indices
start_indices = zeros(1, 2);
end_indices = zeros(1, 2);

% Loop through the regions to find the biggest differences
for i = 1:size(regions, 1)
    % Extract the curves from Ycompend
    curve1 = Ycompend(regions(i, 1):regions(i, 2), location);
    curve2 = Ycompend(regions(i, 1):regions(i, 2), 1);

    % Find the biggest difference in the current region
    [~, start_indices(i), end_indices(i)] = find_biggest_difference(curve1, curve2, 20);
end

% Build the RamanF curve
RamanF = [...
    Ycompend(1:start_indices(1)-1+regions(1, 1), 1)', ...
    Ycompend(start_indices(1)+regions(1, 1):end_indices(1)+regions(1, 1), location)', ...
    Ycompend(end_indices(1)+1+regions(1, 1):start_indices(2)-1+regions(2, 1), 1)', ...
    Ycompend(start_indices(2)+regions(2, 1):end_indices(2)+regions(2, 1), location)', ...
    Ycompend(end_indices(2)+1+regions(2, 1):end, 1)' ...
];

% Plot the results
figure;

plot(Ycompend(:, 1));
hold on;
plot(RamanF);
legend( 'Original Reference Curve','RamanF');
xlabel('Data Points');
ylabel('Intensity');
title('Comparison of RamanF and Original Reference Curve');
hold off;
% 
% Ramanforwenwrong = [...
%     Y(1:start_indices(1)-1+regions(1, 1), 1)', ...
%     Y(start_indices(1)+regions(1, 1):end_indices(1)+regions(1, 1), location)', ...
%     Y(end_indices(1)+1+regions(1, 1):start_indices(2)-1+regions(2, 1), 1)', ...
%     Y(start_indices(2)+regions(2, 1):end_indices(2)+regions(2, 1), location)', ...
%     Y(end_indices(2)+1+regions(2, 1):end, 1)' ...
% ];
% 
% for i=1:length(n)
% Y_wenwrong(:,i)=Y(:,n(i))-Ramanforwenwrong';
% [Base, Y_baseremoved(:,n(i))] = baseline(double(Y_wenwrong(:,i)'));
%     
%     
% % Refine the Y data using the specified column as reference
% [Ycompend(:,n(i)), scaling_factor(:,n(i))] = refine_y_up_multispecs(Y_baseremoved(:,n(i)), Y_baseremoved(:, compnum), idx1, idx2, idx3, idx4);
% end
%% suntract Raman from here on only valid spectrums are processed, if the plotting looks wrong at the stoke side, redo the raman
% Ensure RamanF is a column vector for correct broadcasting
clear Ycompend_subtracted
RamanF = RamanF(:);
% RamanF=RamanF-562;
% Subtract RamanF from each column of Ycompend
Ycompend_subtracted = bsxfun(@minus, Ycompend(:,valid_indices), RamanF);

% Alternatively, if you are using a recent version of MATLAB, you can do:
% Ycompend_subtracted = Ycompend - RamanF;

% Plot the subtracted result for the first few columns as an example
figure;
plot(Ycompend_subtracted(:,1:5:end)); % Plotting first 5 columns as an example
hold on 
plot(RamanF,'k','LineWidth',2)
xlabel('Data Points');
ylabel('Subtracted Intensity');
title('Ycompend Columns After Subtracting RamanF');
%% 
figure;plot(Ycompend_subtracted(:,540:670))
%% 
%% 
figure;
surf(Ycompend_subtracted)
shading interp;
view(2);
caxis([0 50]);
view(2);
axis tight
%% seperate SFG and DFG sides + find peaks + find maximum + find guassian fit + intergration of guassian fit

plot_with_gradient_and_Raman(X_cm, Ycompend, RamanF,powerstep_onlyvalid,power785)
% Define the filename for saving
filename = 'fullrawplot_smaller';
% Save the figure with high resolution
% print(gcf, filename, '-dpng','-r1200');
print(gcf, filename, '-dsvg');
%% 
figure;plot(Ycompend_subtracted(:,271:end))
%% 
close all
clear DFG_raw  SFG_raw  SFG_fit  DFG_fit SFG_maximums DFG_maximums SFG_integral  DFG_integral
% SFG_range = 305:440;
% DFG_range = 1350:1530;
peak1=531;peak2=1444;
power785=0.33*10^3*0.03; %MW to UW
SFG_range = 1:peak1+10;%11447-10:1950;%1385:1525;
DFG_range = peak2-10:1960;%1385:1525;


% range=1:735
% [DFG_raw, SFG_raw, SFG_fit, DFG_fit,SFG_maximums,DFG_maximums,SFG_integral(range), ~]...
%     = analyze_peaks(Ycompend_subtracted(:,range), X_cm,powerstep_onlyvalid(range,:),SFG_range,DFG_range,power785);




% first_SFG_peak = 533;%501;
first_SFG_peak = peak1;%501;
% first_DFG_peak = 1446;
first_DFG_peak = peak2;
% [~, SFG_results] = analyze_peaks_guided(Ycompend_subtracted, X_cm, powerstep_onlyvalid, SFG_range, DFG_range, power785, first_SFG_peak, first_DFG_peak);
% plot_peaks_guided(SFG_results, DFG_results);
% SFG_results = analyze_peaks_guided_SFG(Ycompend_subtracted, X_cm, powerstep_onlyvalid, SFG_range, power785, first_SFG_peak)
% DFG_results = analyze_peaks_guided_DFG(Ycompend_subtracted(:,1:end), X_cm, powerstep_onlyvalid(1:end,:), DFG_range, power785, first_DFG_peak)
SFG_results1 = analyze_peaks_guided_SFG_noSmoothing(Ycompend_subtracted(:,2:270), X_cm, powerstep_onlyvalid(2:270,:), SFG_range, power785, first_SFG_peak);
peak1=411;SFG_range = 1:peak1+30;first_SFG_peak=peak1;
SFG_results2 = analyze_peaks_guided_SFG_noSmoothing(Ycompend_subtracted(:,271:end), X_cm, powerstep_onlyvalid(271:end,:), SFG_range, power785, first_SFG_peak);
SFG_results=[SFG_results1, SFG_results2];

% [SFG_raw,SFG_fit,DFG_fit,DFG_integral, SFG_integral, SFG_maximums, DFG_maximums] = ...
%     analyze_peakspowerdependent(Ycompend_subtracted, X_cm, powerstep_onlyvalid, SFG_range, DFG_range, power785);

DFG_results = analyze_peaks_guided_DFG_noSmoothing(Ycompend_subtracted(:,1:end), X_cm, powerstep_onlyvalid(1:end,:), DFG_range, power785, first_DFG_peak);
% 
% DFG_results2 = analyze_peaks_guided_DFG_noSmoothing(Ycompend_subtracted(:,232:232+391), X_cm, powerstep_onlyvalid(232:232+391,:), DFG_range, power785, first_DFG_peak);
DFG_results1 = analyze_peaks_guided_DFG_noSmoothing(Ycompend_subtracted(:,1:270), X_cm, powerstep_onlyvalid(1:270,:), DFG_range, power785, first_DFG_peak);
peak2=1610;DFG_range = peak2-30:peak2+200;first_DFG_peak=peak2;
DFG_results2 = analyze_peaks_guided_DFG_noSmoothing(Ycompend_subtracted(:,271:271+272), X_cm, powerstep_onlyvalid(271:271+272,:), DFG_range, power785, first_DFG_peak);
peak2=1782;DFG_range = peak2-5:peak2+200;first_DFG_peak=peak2;
DFG_results3 = analyze_peaks_guided_DFG_noSmoothing(Ycompend_subtracted(:,544:end), X_cm, powerstep_onlyvalid(544:end,:), DFG_range, power785, first_DFG_peak);
peak2=1856;DFG_range = peak2-5:peak2+130;first_DFG_peak=peak2;
DFG_results4 = analyze_peaks_guided_DFG_noSmoothing(Ycompend_subtracted(:,650:end), X_cm, powerstep_onlyvalid(650:end,:), DFG_range, power785, first_DFG_peak);

DFG_results=[DFG_results1, DFG_results2,DFG_results3];
% Assuming you have the figure handles saved when the figures were created
% Calculate ratio of SFG and DFG integrals
%% 
% Calculate ratio of SFG and DFG integrals
% ratio = SFG_integral./ DFG_integral;

% Define colors for the curves
colorNames = {'skyBlue', 'darkblue', 'pink2'};
plot_colors = createColorGradient(colorNames, 3);

% Create a new figure for the ratio plot
figure_ratio = figure; 
% set(figure_ratio, 'Units', 'centimeters', 'Position', [15, 10, 4, 3.5]); % Set figure size
hold on;

% Plot DFG Integral (normalized and scaled)
% h1 = plot(powerstep_onlyvalid(:, 1), DFG_integral/max(DFG_integral), 'Color', plot_colors(1,:), 'LineWidth', 1.5);
% Plot SFG Integral (normalized and shifted)
h2 = plot(powerstep_onlyvalid(35:end, 1), SFG_integral/max(SFG_integral) , 'Color', plot_colors(2,:), 'LineWidth', 1.5); %+ 0.5

h3 = plot(powerstep_onlyvalid1600(:, 1), SFG_integral1600/max(SFG_integral1600), 'Color', plot_colors(3,:), 'LineWidth', 1.5); %+ 0.5
% Plot SFG/DFG ratio (normalized and shifted)
% h3 = plot(powerstep_onlyvalid(:, 1), ratio/max(ratio) + 1, 'Color', plot_colors(3,:), 'LineWidth', 1.5);
% h3 = plot(powerstep_onlyvalid(:, 1), ratio, 'Color', plot_colors(3,:), 'LineWidth', 1.5);
% Plot Raman data (scaled and shifted)
h4 = plot(X_cm, RamanF / 2025 * 5, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
% h3 = plot(r44(:,1), (1-r44(:,4))*0.5+0.51, 'Color', plot_colors(3,:), 'LineWidth', 1.5);
% h5 = plot(r33(:,1), (1-r33(:,2))+0.55, 'Color', [0 0.8 0.2], 'LineWidth', 1.5);

% Adjust axis labels and limits
xlabel('Wavenumber (cm^{-1})');
ylabel('Intensity (a.u.)'); % Since all data is normalized
% xlim([960 1160]);
axis tight
xlim([880 1170]);
xlim([880 1680]);
% ylim([0 1.5]);  % Adjust based on your data range
set(gca, 'XTick', 880:40:1680, 'YTick', linspace(0, 2.5, 6));  % Adjust the Y ticks

% Remove grid lines
grid off;

% Add text labels to the right of each curve to represent the legend
% Adjust the positions of the labels based on the data and your preferences
% text(powerstep(end, 1)-25, DFG_integral(end)/max(DFG_integral)+0.2, 'DFG ', 'Color', plot_colors(1,:), 'FontSize', 8, 'VerticalAlignment', 'middle');
text(powerstep(end, 1) -25, SFG_integral(end)/max(SFG_integral) + 0.7, 'SFG ', 'Color', plot_colors(2,:), 'FontSize', 8, 'VerticalAlignment', 'middle');
% text(powerstep(end, 1) -45, ratio(end)/max(ratio) + 1.15, 'SFG/DFG', 'Color', plot_colors(3,:), 'FontSize', 8, 'VerticalAlignment', 'middle');
text(X_cm(1660) - 20, RamanF(1505)/2025*20+0.1, 'Raman', 'Color', [0.5 0.5 0.5], 'FontSize', 8, 'VerticalAlignment', 'middle');
% text(powerstep(end, 1) -25, SFG_integral(end)/max(SFG_integral) + 0.4, 'FTIR ', 'Color', plot_colors(3,:), 'FontSize', 8, 'VerticalAlignment', 'middle');

% Hold off to finalize plotting
hold off;
title('NTP NP functionlized 1600cm^{-1} sample')
%% Calculate valid integrals and ratio
%% Calculate valid integrals and ratio (same as before)

%% Calculate valid integrals and ratio (with outlier removal option)
n_samples=length(SFG_results);
% Initialization
% n_samples = max(length(SFG_results), length(DFG_results));
SFG_integral = nan(1, n_samples);
DFG_integral = nan(1, n_samples);

% Extract integrals
for i = 1:length(SFG_results)
    if ~isempty(SFG_results{i})
        SFG_integral(i) = SFG_results{i}.integral;
    end
end

for i = 1:length(DFG_results)
    if ~isempty(DFG_results{i})
        DFG_integral(i) = DFG_results{i}.integral;
    end
end

% Find valid integrals
valid_idx = ~isnan(SFG_integral) & ~isnan(DFG_integral);
valid_idx = ~isnan(SFG_integral) & ~isnan(SFG_integral);
SFG_integral_valid = SFG_integral(valid_idx);
DFG_integral_valid = DFG_integral(valid_idx);

% Plot SFG and DFG integrals separately (no x axis, just index)
figure;
subplot(2,1,1);
plot(SFG_integral_valid./DFG_integral_valid, 'o-');
title('SFG Integrals (by index)');
ylabel('Integral');
xlabel('Sample index');

subplot(2,1,2);
plot(DFG_integral_valid, 'o-');
title('DFG Integrals (by index)');
ylabel('Integral');
xlabel('Sample index');

% Ask user to input bad sample indices
disp('If you see bad peaks, input their indices as a vector (e.g., [3 7 12])');
bad_indices = input('Enter bad sample indices to remove (if none, enter []): ');

x_axis_clean = powerstep_onlyvalid(valid_idx);  % powerstep_onlyvalid is your X-axis

if ~isempty(bad_indices)
    good_indices = true(1, length(SFG_integral_valid));
    good_indices(bad_indices) = false;

    SFG_integral_valid_clean = SFG_integral_valid(good_indices);
    DFG_integral_valid_clean = DFG_integral_valid(good_indices);
    ratio_valid_clean = SFG_integral_valid_clean ./ DFG_integral_valid_clean;
    x_axis_clean = x_axis_clean(good_indices);  % Clean x-axis too
else
    % If no bad indices, use as is
    SFG_integral_valid_clean = SFG_integral_valid;
    DFG_integral_valid_clean = DFG_integral_valid;
    ratio_valid_clean = SFG_integral_valid_clean ./ DFG_integral_valid_clean;
end
%% 

% Colors
colorNames = {'skyBlue', 'darkblue', 'pink2'};
plot_colors = createColorGradient(colorNames, 3);

% Create figure
figure_ratio = figure; 
% set(figure_ratio, 'Units', 'centimeters', 'Position', [15, 10, 4, 3.5]);
hold on;

% Plot all data (normalized and shifted)
h1 = plot(x_axis_clean(1:end), DFG_integral_valid_clean(1:end) / max(DFG_integral_valid_clean(1:end))+0.2, 'Color', plot_colors(1,:), 'LineWidth', 1.5);
h2 = plot(x_axis_clean, SFG_integral_valid_clean / max(SFG_integral_valid_clean) + 0.5, 'Color', plot_colors(2,:), 'LineWidth', 1.5);
h3 = plot(x_axis_clean(1:end), ratio_valid_clean(1:end)/max(ratio_valid_clean(1:end))+1.2, 'Color', plot_colors(3,:), 'LineWidth', 1.5);
h4 = plot(X_cm, RamanF / 2025 , 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

% Labels
xlabel('Wavenumber (cm^{-1})');
ylabel('Intensity (a.u.)');
axis tight
xlim([880 1680]);
% ylim([0 2]);
set(gca, 'XTick', 880:40:1680, 'YTick', linspace(0, 2.5, 6));

grid off;

% Text annotations
text(1630-25, DFG_integral_valid_clean(end)/max(DFG_integral_valid_clean)+0.2, 'DFG', 'Color', plot_colors(1,:), 'FontSize', 8, 'VerticalAlignment', 'middle');
text(1630-25, SFG_integral_valid_clean(end)/max(SFG_integral_valid_clean) + 0.5, 'SFG', 'Color', plot_colors(2,:), 'FontSize', 8, 'VerticalAlignment', 'middle');
text(1630-80, ratio_valid_clean(end-170)/max(ratio_valid_clean(end-170)) + 0.8, 'SFG/DFG', 'Color', plot_colors(3,:), 'FontSize', 8, 'VerticalAlignment', 'middle');
text(X_cm(1660) - 20, RamanF(1505)/2025*10+0.1, 'Raman', 'Color', [0.5 0.5 0.5], 'FontSize', 8, 'VerticalAlignment', 'middle');

hold off;
% title('1600 cm^{-1} resonant sample')
title('1080 cm^{-1} resonant BPT sample')
%% 
figure;
subplot(3,1,1)
plot(DFG_integral_valid_clean,'Color', plot_colors(1,:),'LineWidth',1.5);axis tight;
subplot(3,1,2)
plot(SFG_integral_valid_clean,'Color', plot_colors(2,:),'LineWidth',1.5);axis tight;
subplot(3,1,3)
plot(ratio_valid_clean,'Color', plot_colors(3,:),'LineWidth',1.5);axis tight;
%% Plot 1: DFG Integral

figure;
plot(powerstep_onlyvalid(:,1), DFG_integral_valid / max(DFG_integral_valid), 'b-', 'LineWidth', 1.5);
xlabel('Sample Index');
ylabel('Normalized DFG Integral');
title('DFG Integral');
ylim([0 1.2]);
grid on;

%% Plot 2: SFG Integral

figure;
plot(powerstep_onlyvalid(:,1), DFG_integral_valid / max(DFG_integral_valid), 'm-', 'LineWidth', 1.5);
xlabel('Sample Index');
ylabel('Normalized SFG Integral');
title('SFG Integral');
ylim([0 1.2]);
grid on;

%% Plot 3: SFG/DFG Ratio

figure;
plot(powerstep_onlyvalid(:,1), ratio_valid / max(ratio_valid), 'r-', 'LineWidth', 1.5);
xlabel('Sample Index');
ylabel('Normalized SFG/DFG Ratio');
title('SFG/DFG Ratio');
ylim([0 1.2]);
grid on;


%% 
DFG_norm=DFG_integral/max(DFG_integral);
SFG_norm=SFG_integral/max(SFG_integral);
ratio_norm=ratio/max(ratio);
% Combine the parameters into a table with specific column names
data_table = table(DFG_norm, SFG_norm , ratio_norm, powerstep_onlyvalid(:,1), ...
    'VariableNames', {'DFG', 'SFG', 'intensity', 'x'});

% Create a dynamic file name including the NP value
output_filename = sprintf('shorter_NP%d.xlsx', npnumber);

% Write the table to an Excel file
writetable(data_table, output_filename);

disp(['Data saved to ', output_filename]);


%% 

figurePaths = {
    'SFG_raw_plot_smaller', ...
    'SFG_fit_plot_smaller'
%     'DFG_raw_plot_smaller', ...
%     'DFG_fit_plot_smaller',...
%     'DFG_integration_smaller',...
%     'SFG_integration_smaller',...
%     'include all_smaller'
};

% Find all open figures
figureHandles = findobj('Type', 'Figure');

% Ensure the number of paths matches the number of open figures
if length(figureHandles) ~= length(figurePaths)
    error('Number of figure handles does not match number of paths.');
end

% Save each figure
for i = 1:length(figureHandles)
    % Set the figure handle
    figure(figureHandles(i));
    
%     % Customize figure properties if needed
%     set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [15, 10, 7, 5]);
%     set(groot, 'defaultAxesFontName', 'Helvetica', 'defaultAxesFontSize', 8);
    
    % Save figure as TIFF with high resolution
    print(gcf, figurePaths{i}, '-dpng','-r1200'); % Save as TIFF with 300 DPI
% print(gcf, figurePaths{i}, '-dsvg');
end
%% 
figure;
hold on;

% Plot SFG_integral vs. powerstep_onlyvalid(1,:) using a solid line
plot(powerstep_onlyvalid(:,1), SFG_integral, 'LineWidth', 2);

% Plot RamanF on the same figure using dark gray from the color library
darkGray = DcolorLibrary('darkgray');
plot(X_cm, RamanF/max(RamanF)/10, 'Color', darkGray, 'LineWidth', 1.5);

% Customize the axes and labels
xlabel('Power Step');
ylabel('SFG Integral / RamanF');
set(gca, 'FontName', 'Helvetica', 'FontSize', 9);
xlim([960, 1160]);

% Adjust figure size
set(gcf, 'Units', 'centimeters', 'Position', [15, 10, 14, 10]);

hold off;

%% plot results for checking
% close all
clear length
figure;
subplot(2,3,[1,2]);hold on
for i=1:length(DFG_raw)
plot(DFG_raw{i}.x,DFG_raw{i}.y)
end
hold off;
title('DFG Raw Peaks found');
xlabel('wavenumber (cm-)');
ylabel('Intensity (counts/)');
xlim([960 1160])
subplot(2,3,[4,5]);hold on;
for i=1:length(DFG_fit)
plot(DFG_fit{1,i}.xFit,DFG_fit{1,i}.yFit)
end
hold off;
title('DFG guassian fitted Peaks');
xlabel('wavenumber (cm-)');
ylabel('Intensity (counts/uw)');
xlim([960 1160])


subplot(2,3,[3,6]);hold on
for i=1:length(DFG_maximums)
DFGfitmax(i,:)=[DFG_maximums{i}.xmax,DFG_maximums{i}.ymax];
end
plot(powerstep_onlyvalid(:,1),DFG_integral,'DisplayName', 'MIRcat wavenumber vs. peak intergration') %sprintf('DFG Fit %d', i)
% plot(DFGfitmax(:,1),DFGfitmax(:,2),'DisplayName', 'fitted peak wavenumber vs. fitted peak max')
plot(powerstep_onlyvalid(:,1),powerstep_onlyvalid(:,2)/max(powerstep_onlyvalid(:,2)),'DisplayName', 'powermeter') %sprintf('DFG Fit %d', i)
legend
xlabel('wavenumber (cm-)');
ylabel('Intensity  (counts/uw)');
xlim([960 1160])

clear length
figure;
subplot(2,3,[1,2]);hold on
for i=1:length(SFG_raw)
plot(SFG_raw{i}.x,SFG_raw{i}.y)
end
hold off;
title('SFG Raw Peaks found');
xlabel('wavenumber (cm-)');
ylabel('Intensity  (counts/uw)');
xlim([-1160 -950])
subplot(2,3,[4,5]);hold on;
for i=1:length(SFG_fit)
plot(SFG_fit{i}.xFit,SFG_fit{i}.yFit)
end
hold off;
title('SFG guassian fitted Peaks');
xlabel('wavenumber (cm-)');
ylabel('Intensity  (counts/uw)');
xlim([-1160 -950])


subplot(2,3,[3,6]);hold on
for i=1:length(SFG_maximums)
SFGfitmax(i,:)=[-SFG_maximums{i}.xmax,SFG_maximums{i}.ymax];
end
plot(powerstep_onlyvalid(:,1),SFG_integral,'DisplayName', 'MIRcat wavenumber vs. peak intergration') %sprintf('DFG Fit %d', i)
% plot(SFGfitmax(:,1),SFGfitmax(:,2),'DisplayName', 'fitted peak wavenumber vs. fitted peak max')
plot(powerstep_onlyvalid(:,1),powerstep_onlyvalid(:,2)/max(powerstep_onlyvalid(:,2)),'DisplayName', 'powermeter') %sprintf('DFG Fit %d', i)
legend
xlabel('wavenumber (cm-)');
ylabel('Intensity / Integral');
xlim([960 1160])

%% 
% process the noisy results
SFG_denoise0 = wdenoise(SFG_integral,3);
DFG_denoise0 = wdenoise(DFG_integral,2);
ratio_denoise0 = wdenoise(ratio,2);
figure;hold on
plot(powerstep_onlyvalid(:,1),SFG_denoise0,'DisplayName', 'peak intergration desnoized') %sprintf('DFG Fit %d', i)
plot(powerstep_onlyvalid(:,1),SFG_integral,'DisplayName', 'peak intergration')
plot(powerstep_onlyvalid(:,1),powerstep_onlyvalid(:,2)/max(powerstep_onlyvalid(:,2)),'DisplayName', 'powermeter') %sprintf('DFG Fit %d', i)

% plot(X_cm,RamanF/max(RamanF),'DisplayName', 'Raman')
% plot(SFGfitmax(:,1),SFGfitmax(:,2),'DisplayName', 'fitted peak max vs. peak intergration')
legend
xlabel('wavenumber (cm-)');
ylabel('Intensity / Integral');
xlim([960 1160])
title('SFG')

figure;hold on
plot(powerstep_onlyvalid(:,1),DFG_denoise0,'DisplayName', 'peak intergration desnoized') %sprintf('DFG Fit %d', i)
plot(powerstep_onlyvalid(:,1),DFG_integral,'DisplayName', 'peak intergration')
plot(powerstep_onlyvalid(:,1),powerstep_onlyvalid(:,2)/max(powerstep_onlyvalid(:,2)),'DisplayName', 'powermeter') %sprintf('DFG Fit %d', i)

% plot(X_cm,RamanF/max(RamanF),'DisplayName', 'Raman')
% plot(SFGfitmax(:,1),SFGfitmax(:,2),'DisplayName', 'fitted peak max vs. peak intergration')
legend
xlabel('wavenumber (cm-)');
ylabel('Intensity / Integral');
title('DFG')
xlim([960 1160])
%% 
%% 

[filepath, name, ext] = fileparts(filenamespectrum);
% Extract the directory name (e.g., 'NP 11') from the filepath
[filepath_base, last_dir] = fileparts(filepath);

% Combine the directory name, current filename, and additional string to form the new filename


save(fullfile(filepath, sprintf('%s_%s_CWstepping_doubt_newmethod2.mat', last_dir, name)));





