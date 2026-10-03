function results = Wavelet_Spectrum_Shifter(sig_x, t)
    % 1. Initialization & Configuration
    p = get_default_parameters();
    
    % 3. Processing
    results = run_pipeline(sig_x, t, p);

    % 4. Summary display
    print_summary(results);

    % 5. Persist results for use in other scripts (overwrites any
    %    previous pipeline_results.mat - each run starts fresh)
    %save_results(results, 'C:\Users\moham\OneDrive\Desktop\stage dété\codes\CWT_analyserpipeline_results.mat');
end