function save_results(results, ~)
    % Get the folder where this specific script/function is saved
    [current_folder, ~, ~] = fileparts(mfilename('fullpath'));
    
    % Combine the folder path with your desired filename
    forced_filepath = fullfile(current_folder, 'results.mat');
    
    % Save the results struct
    save(forced_filepath, 'results');
    
    fprintf('\n[SAVED] results struct written to %s\n', forced_filepath);
end