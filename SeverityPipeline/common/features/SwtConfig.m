classdef SwtConfig
    % Parameters of pipeline 2 (bilateral Morse CWT + synchrosqueezing)
    properties
        Fs = 10000
        TargetFreq = -200
        Gamma = 3
        BetaOpt = 120
        VoicesPerOctave = 14
        ThreshRatio = 1e-4
        SearchRadiusHz = 5
        MaxExpansionHz = 3
        WindowSec = 0.2    % seconds (cyclic window)
    end
end
