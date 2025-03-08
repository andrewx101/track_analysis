clearvars
% change file path here
path='C:\Users\andre\OneDrive\Manuscripts\2021 Alginate\01-05-2022 all data\Alginate(micro)\f=0.158';
file='tr-4.mat';

load(fullfile(path,file)); 

[msd,msd_sp]=calculate_MSD(tr);






function [msd_results, single_particle_msds] = calculate_MSD(tracks)  
    % Input: Nx4 matrix [x, y, step_id, track_id]  
    % Output: MSD matrix [delta_step, mean_sqdisp, count]  
    % Output: cell array containing MSD for each trajectory  

    track_ids = unique(tracks(:, 4));  
    num_tracks = length(track_ids);  
    single_particle_msds = cell(num_tracks, 1);  
    
    % Preallocate arrays for results  
    all_delta_steps_cell = cell(num_tracks, 1);  
    all_sq_disps_cell = cell(num_tracks, 1);  

    parfor k = 1:num_tracks  
        track_id = track_ids(k);  
        track_data = tracks(tracks(:, 4) == track_id, :);  
        
        % Sort by step_id and extract coordinates  
        [~, idx] = sort(track_data(:, 3));  
        x = track_data(idx, 1);  
        y = track_data(idx, 2);  
        steps = track_data(idx, 3);  
        n = length(steps);  
        
        if n < 2  
            continue; % Skip if there's not enough data  
        end  
        
        % Generate all valid step pairs  
        [i, j] = find(triu(true(n), 1)); % Indices for upper triangle  
        dx = x(j) - x(i);  
        dy = y(j) - y(i);  
        delta_steps = steps(j) - steps(i);  
        sq_disps = dx.^2 + dy.^2;  

        % Store calculated results  
        all_delta_steps_cell{k} = delta_steps;  
        all_sq_disps_cell{k} = sq_disps;  

        % Calculate single particle MSD  
        single_particle_msd = zeros(max(steps), 3);  
        for step = 1:max(steps)  
            ind = (delta_steps == step);  
            if any(ind)  
                single_particle_msd(step, 1) = step; % lag time  
                single_particle_msd(step, 2) = mean(sq_disps(ind)); % MSD  
                single_particle_msd(step, 3) = sum(ind); % Number of observations  
            end  
        end  
        single_particle_msds{k} = single_particle_msd;  
    end  
    
    % Concatenate results after parfor loop  
    all_delta_steps = vertcat(all_delta_steps_cell{:});  
    all_sq_disps = vertcat(all_sq_disps_cell{:});  
    
    % Calculate overall statistics  
    [unique_deltas, ~, idx] = unique(all_delta_steps);  
    sum_sqdisp = accumarray(idx, all_sq_disps, [], @sum);  
    counts = accumarray(idx, all_sq_disps, [], @numel);  
    
    % Create output matrix for overall MSD  
    msd_results = [unique_deltas, sum_sqdisp ./ counts, counts];  
    msd_results = sortrows(msd_results, 1);  
end    