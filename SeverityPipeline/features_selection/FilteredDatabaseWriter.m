classdef FilteredDatabaseWriter
    % Same file layout as the Python scripts: feat_table, names (1xN cell), percent (row) [, degree]
    methods (Static)
        function write(outFile, X, names, cols, percent, degree)
            feat_table = X(:, cols); %#ok<NASGU>
            names = names(cols); %#ok<NASGU>
            percent = double(percent(:))'; %#ok<NASGU>
            if isempty(degree)
                save(outFile, 'feat_table', 'names', 'percent');
            else
                degree = degree(:)'; %#ok<NASGU>
                save(outFile, 'feat_table', 'names', 'percent', 'degree');
            end
            fprintf('\nSaved filtered database to: %s\nShape: %d x %d\n', outFile, size(X, 1), numel(cols));
        end
    end
end