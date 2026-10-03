classdef FeatureConfig
    % Parameters of pipeline 1 (I1/I2/CSV window features)
    properties
        Fs = 10000
        F0 = 200            % main freq tracked by an FFT Analysis before processing the database 
        H = 5
        WindowLength = 0.2   % seconds (cyclic window)
    end
end
