%% Estimate FOV-wise average alpha_2

% alpha_cell=cell(N,1); N is the number of FOVs
% put the alpha_2 data from each FOV to a cell array alpha_cell manually

data_cell=alpha_cell;

% Collect all x-values from all matrices
all_x=cell(length(data_cell),1);


for i = 1:length(data_cell)
    all_x{i} = data_cell{i}(:,1);

end

all_x=cell2mat(all_x);

% Determine the new x-range (you can adjust the step size)
min_x = min(all_x);
max_x = max(all_x);
digits=-(floor(log10(min_x))-3);
all_x=round(all_x,digits);
% step = 0.1; % Choose an appropriate step size based on your data density
step=estimate_step_size(all_x);
new_x = (min_x:step:max_x)';

% Initialize arrays for interpolated values
interp_y = zeros(length(new_x), length(data_cell));
interp_std = zeros(length(new_x), length(data_cell));

for i = 1:length(data_cell)
    current_data = data_cell{i};
    x_orig = current_data(:,1);
    y_orig = current_data(:,2);
    std_orig = current_data(:,3);
    
    % Interpolate y values
    interp_y(:,i) = interp1(x_orig, y_orig, new_x, 'pchip');
    
    % Interpolate standard deviations
    interp_std(:,i) = interp1(x_orig, std_orig, new_x, 'pchip');
end

% Calculate combined mean
combined_mean = mean(interp_y, 2);

% Calculate combined standard error
n = size(interp_y, 2); % number of datasets
between_var = var(interp_y, 0, 2); % variance between datasets
within_var = mean(interp_std.^2, 2); % average of squared std errors
combined_std = sqrt(between_var + within_var);

final_matrix = [new_x, combined_mean, combined_std];

% Remove any NaN values that might occur at the edges
final_matrix = final_matrix(all(~isnan(final_matrix(:,2:3)), 2), :);



function optimal_step = estimate_step_size(all_x)
    
    % Remove duplicates and sort
    unique_x = unique(all_x);
    sorted_x = sort(unique_x);
    
    % Calculate differences between consecutive points
    diffs = diff(sorted_x);
    
    % Handle case where all x-values are identical
    if all(diffs == 0)
        optimal_step = 1; % Default step size when all x-values are identical
        return;
    end
    
    % Calculate robust statistics
    min_diff = min(diffs(diffs > 0)); % Smallest non-zero difference
    median_diff = median(diffs);
    mean_diff = mean(diffs);
    std_diff = std(diffs);
    
    % Calculate density-based factors
    density_factor = 1/(mean_diff + std_diff);
    cluster_factor = sum(diffs < mean_diff/10)/length(diffs); % Measure of clustering
    
    % Calculate base step size
    if cluster_factor > 0.8 % Highly clustered data
        base_step = median_diff;
    else
        base_step = min_diff + 0.5*(median_diff - min_diff)*density_factor;
    end
    
    % Apply reasonable bounds
    data_range = max(all_x) - min(all_x);
    optimal_step = max(base_step, min_diff);
    optimal_step = min(optimal_step, data_range/20); % Prevent too large steps
    
    % Final adjustment based on number of unique points
    n_unique = length(unique_x);
    if n_unique < 10 % Very few unique points
        optimal_step = optimal_step * 0.8; % Use finer steps
    end
end
