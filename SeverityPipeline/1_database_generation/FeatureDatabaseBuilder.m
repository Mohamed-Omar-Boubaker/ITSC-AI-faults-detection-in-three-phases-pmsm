classdef FeatureDatabaseBuilder
    % Orchestrator: read -> resample ONCE -> run every extractor -> merge -> save
    properties (SetAccess = immutable)
        Source
        Extractors
        Writer
        FsTarget
    end
    methods
        function obj = FeatureDatabaseBuilder(source, extractors, writer, fsTarget)
            obj.Source = source; obj.Extractors = extractors;
            obj.Writer = writer; obj.FsTarget = fsTarget;
        end
        function run(obj, startIdx)
            N = obj.Source.count();
            names = {}; featTable = [];
            degree = zeros(N, 1); percent = zeros(N, 1);
            if startIdx > 1, [featTable, names, degree, percent] = obj.Writer.loadExisting(); end  % resume

            for i = startIdx:N
                [raw, degree(i), percent(i)] = obj.Source.read(i);
                cur = raw.resampleTo(obj.FsTarget);       % both pipelines get the same signal
                sets = cellfun(@(e) e.extract(cur), obj.Extractors, 'UniformOutput', false);
                fs = FeatureSet.merge(sets);
                if isempty(featTable)
                    featTable = zeros(N, numel(fs.Values));
                    names = fs.Names;
                end
                featTable(i, :) = fs.Values;
                obj.Writer.write(featTable, names, degree, percent);
                fprintf('sample %d / %d done\n', i, N);
            end
            disp('finished');
        end
    end
end
