function plotDataFields(dataSave, fieldSets, xlabels, ylabels, addons)
    % Create a new figure
    figure;

    % Number of fields to plot
    numFields = length(fieldSets);

    % Create tiled layout for plots
    tiledlayout(numFields, 1);

    % Loop over each field and plot
    for i = 1:numFields
        fieldSet = fieldSets{i};

        % Get current subplot
        ax = nexttile;

        % Check if the fields exist in the dataSave structure
        for f = 1:length(fieldSet)
            if ~isfield(dataSave, fieldSet{f})
                error(['Field "', fieldSet{f}, '" does not exist in dataSave.']);
            end
        end

        % Determine the type of plot (2D or 3D)
        if length(fieldSet) == 2
            xData = dataSave.(fieldSet{1});
            yData = dataSave.(fieldSet{2});

            [dataLength, numScenarios, numMethods] = size(yData);

            % Set colors and line styles for different methods
            lineColors = {'b', 'r'};
            lineStyles = {'-', '--'};

            % Plot each scenario and method
            for j = 1:numScenarios
                for method = 1:numMethods
                    plot(ax, xData(:, j, method), yData(:, j, method), ...
                        'Color', lineColors{method}, 'LineStyle', lineStyles{method}, ...
                        'LineWidth', 2);
                    hold(ax, 'on');
                end
            end

        elseif length(fieldSet) == 3
            xData = dataSave.(fieldSet{1});
            yData = dataSave.(fieldSet{2});
            zData = dataSave.(fieldSet{3});

            [dataLength, numScenarios, numMethods] = size(zData);

            % Set colors and line styles for different methods
            lineColors = {'b', 'r'};
            lineStyles = {'-', '--'};

            % Plot each scenario and method
            for j = 1:numScenarios
                for method = 1:numMethods
                    plot3(ax, xData(:, j, method), yData(:, j, method), zData(:, j, method), ...
                        'Color', lineColors{method}, 'LineStyle', lineStyles{method}, ...
                        'LineWidth', 2);
                    hold(ax, 'on');
                end
            end
        else
            error('Each field set should contain either 2 or 3 field names.');
        end

        % If addons are provided, plot them
        if ~isempty(addons) && length(addons) >= i && ~isempty(addons{i})
            if length(fieldSet) == 2
                plot(ax, addons{i}.t, addons{i}.data, 'k:', 'LineWidth', 2);
            else
                plot3(ax, addons{i}.t, addons{i}.data, addons{i}.zdata, 'k:', 'LineWidth', 2);
            end
        end

        % Set labels
        xlabel(ax, xlabels{i});
        ylabel(ax, ylabels{i});
    end
    hold(ax, 'off');
end


