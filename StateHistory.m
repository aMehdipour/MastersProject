classdef StateHistory < dynamicprops
    properties
        % time                   = []
        % flightPathAngle        = []
        % rangeToGo              = []
        % altitudeUnnormalized   = []
        % velocityUnnormalized   = []
        % rangeToGoUnnormalized  = []
        % normGeocentricDistance = []
        % velocity               = []
        % latitude               = []
        % longitude              = []
        % heading                = []
        % posECEF                = []
        % bankAngle              = []
        % lift                   = []
        % drag                   = []
        % loadFactor             = []
        % heatingRate            = []
        % commandedBankAngle     = []
        % headingError           = []
        % crossrange             = []
    end

    methods
        function obj = StateHistory()
            % Constructor intentionally left blank
        end

        % function obj = appendState(obj, currentTime, state)
        %     obj.time                   = [obj.time; currentTime];
        %     obj.flightPathAngle        = [obj.flightPathAngle; state.flightPathAngle];
        %     obj.rangeToGo              = [obj.rangeToGo; state.rangeToGo];
        %     obj.altitudeUnnormalized   = [obj.altitudeUnnormalized; state.altitudeUnnormalized];
        %     obj.velocityUnnormalized   = [obj.velocityUnnormalized; state.velocityUnnormalized];
        %     obj.normGeocentricDistance = [obj.normGeocentricDistance; state.normGeocentricDistance];
        %     obj.velocity               = [obj.velocity; state.velocity];
        %     obj.latitude               = [obj.latitude; state.latitude];
        %     obj.longitude              = [obj.longitude; state.longitude];
        %     obj.heading                = [obj.heading; state.heading];
        %     obj.posECEF                = [obj.posECEF, state.posECEF];
        %     obj.bankAngle              = [obj.bankAngle; state.bankAngle];
        %     obj.lift                   = [obj.lift; state.lift];
        %     obj.drag                   = [obj.drag; state.drag];
        %     obj.loadFactor             = [obj.loadFactor; state.loadFactor];
        %     obj.heatingRate            = [obj.heatingRate; state.heatingRate];
        %     obj.commandedBankAngle     = [obj.commandedBankAngle; state.commandedBankAngle];
        %     obj.rangeToGoUnnormalized  = [obj.rangeToGoUnnormalized; state.rangeToGoUnnormalized];
        %     obj.headingError           = [obj.headingError; state.headingError];
        %     obj.crossrange             = [obj.crossrange; state.crossrange];
        % end

        function obj = appendState(obj, currentTime, state)
            % Iterate over the fields of the state struct
            fields = fieldnames(state);
            for i = 1:length(fields)
                field = fields{i};
                % Check if the field exists in the StateHistory object
                if ~isprop(obj, field)
                    % Add the field dynamically using addprop
                    addprop(obj, field);
                    obj.(field) = [];
                end
                % Append the value to the corresponding array
                obj.(field) = [obj.(field); state.(field)];
            end

            % Append the current time
            if ~isprop(obj, 'time')
                addprop(obj, 'time');
                obj.time = [];
            end
            obj.time = [obj.time; currentTime];
        end

        function plotState(obj)
            figure('Name', 'Flight Path Angle');
            plot(obj.time, obj.flightPathAngle);
            title('Flight Path Angle vs. Time');
            xlabel('Time (s)');
            ylabel('Angle (rad)');

            figure('Name', 'Range to Go');
            plot(obj.time, obj.rangeToGoUnnormalized);
            title('Range to Go vs. Time');
            xlabel('Time (s)');
            ylabel('Range (m)');

            figure('Name', 'Altitude');
            plot(obj.time, obj.altitudeUnnormalized);
            title('Altitude vs. Time');
            xlabel('Time (s)');
            ylabel('Altitude (m)');

            figure('Name', 'Velocity');
            plot(obj.time, obj.velocityUnnormalized);
            title('Velocity vs. Time');
            xlabel('Time (s)');
            ylabel('Velocity (m/s)');

            figure('Name', 'Norm Geocentric Distance');
            plot(obj.time, obj.normGeocentricDistance);
            title('Norm Geocentric Distance vs. Time');
            xlabel('Time (s)');
            ylabel('Distance (m)');

            figure('Name', 'Normalized Velocity');
            plot(obj.time, obj.velocity);
            title('Normalized Velocity vs. Time');
            xlabel('Time (s)');
            ylabel('Velocity');

            figure('Name', 'Latitude');
            plot(obj.time, rad2deg(obj.latitude));
            title('Latitude vs. Time');
            xlabel('Time (s)');
            ylabel('Latitude (degrees)');

            figure('Name', 'Longitude');
            plot(obj.time, rad2deg(obj.longitude));
            title('Longitude vs. Time');
            xlabel('Time (s)');
            ylabel('Longitude (degrees)');

            figure('Name', 'Heading');
            plot(obj.time, obj.heading);
            title('Heading vs. Time');
            xlabel('Time (s)');
            ylabel('Heading (rad)');

            figure('Name', 'Bank Angle');
            hold on
            plot(obj.time, obj.bankAngle);
            plot(obj.time, obj.commandedBankAngle, '--');
            title('Bank Angle vs. Time');
            xlabel('Time (s)');
            ylabel('Bank Angle (rad)');
            legend('Bank Angle', 'Commanded Bank Angle')

            figure('Name','Normalized Lift And Drag vs Range')
            hold on
            plot(obj.time, obj.drag)
            plot(obj.time, obj.lift)
            plot(obj.time, sqrt(obj.lift.^2 + obj.drag.^2))
            xlabel('Normalized Range to Go')
            ylabel('Lift and Drag Accelerations (g)')
            legend('Lift', 'Drag', 'Total Load')

            figure('Name', 'Heating Rate History')
            plot(obj.time(2:end), obj.heatingRate(2:end))
            xlabel('Normalized Range to Go')
            ylabel('Heating Rate $(W/cm^2)$')

            figure('Name', 'Velocity vs Altitude')
            plot(obj.velocityUnnormalized, obj.altitudeUnnormalized./1e3)
            xlabel('Velocity m/s')
            ylabel('Altitude (km)')

            figure('Name', 'Heading Error')
            plot(obj.time, obj.headingError)
            xlabel('Time (s)')
            ylabel('Heading Error (rad)')

            figure('Name', 'Crossrange')
            plot(obj.time, obj.crossrange)
            xlabel('Time (s)')
            ylabel('Crossrange (m)')

            % figure('Name', 'ECEF Position Components');
            % plot(obj.time, obj.posECEF(1, :), 'DisplayName', 'X');
            % hold on;
            % plot(obj.time, obj.posECEF(2, :), 'DisplayName', 'Y');
            % plot(obj.time, obj.posECEF(3, :), 'DisplayName', 'Z');
            % hold off;
            % legend show;
            % title('ECEF Position Components vs. Time');
            % xlabel('Time (s)');
            % ylabel('Position (m)');
        end
        function plotStateAndSave(obj)
            baseDir = 'Plots';
            % Create a subdirectory name based on the current date and time
            % datetime object with a custom format
            dt = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
            subDir = char(dt);  % Convert datetime to character array suitable for file paths
            fullDir = fullfile(baseDir, subDir);
            if ~exist(fullDir, 'dir')
                mkdir(fullDir);
            end
            
            plotFunctions = {@plotFlightPathAngle, @plotRangeToGo, @plotAltitude, ...
                             @plotVelocity, @plotNormGeocentricDistance, @plotNormalizedVelocity, ...
                             @plotLatitude, @plotLongitude, @plotHeading, @plotBankAngle, ...
                             @plotLiftDrag, @plotHeatingRate, @plotVelocityAltitude, @plotECEFPosition};
            plotNames = {'Flight Path Angle', 'Range to Go', 'Altitude', 'Velocity', ...
                         'Norm Geocentric Distance', 'Normalized Velocity', 'Latitude', ...
                         'Longitude', 'Heading', 'Bank Angle', 'Normalized Lift And Drag vs Range', ...
                         'Heating Rate History', 'Velocity vs Altitude', 'ECEF Position Components'};
                         
            for i = 1:length(plotFunctions)
                fig = figure('Visible', 'off');
                plotFunctions{i}(obj);
                saveas(fig, fullfile(fullDir, [plotNames{i}, '.png']));
                close(fig);
            end
        end
        
        function plotFlightPathAngle(obj)
            plot(obj.time, obj.flightPathAngle);
            title('Flight Path Angle vs. Time');
            xlabel('Time (s)');
            ylabel('Angle (rad)');
        end
        
        function plotRangeToGo(obj)
            plot(obj.time, obj.rangeToGoUnnormalized);
            title('Range to Go vs. Time');
            xlabel('Time (s)');
            ylabel('Range (m)');
        end
        
        function plotAltitude(obj)
            plot(obj.time, obj.altitudeUnnormalized);
            title('Altitude vs. Time');
            xlabel('Time (s)');
            ylabel('Altitude (m)');
        end
        
        function plotVelocity(obj)
            plot(obj.time, obj.velocityUnnormalized);
            title('Velocity vs. Time');
            xlabel('Time (s)');
            ylabel('Velocity (m/s)');
        end
        
        function plotNormGeocentricDistance(obj)
            plot(obj.time, obj.normGeocentricDistance);
            title('Norm Geocentric Distance vs. Time');
            xlabel('Time (s)');
            ylabel('Distance (m)');
        end
        
        function plotNormalizedVelocity(obj)
            plot(obj.time, obj.velocity);
            title('Normalized Velocity vs. Time');
            xlabel('Time (s)');
            ylabel('Velocity');
        end
        
        function plotLatitude(obj)
            plot(obj.time, rad2deg(obj.latitude));
            title('Latitude vs. Time');
            xlabel('Time (s)');
            ylabel('Latitude (degrees)');
        end
        
        function plotLongitude(obj)
            plot(obj.time, rad2deg(obj.longitude));
            title('Longitude vs. Time');
            xlabel('Time (s)');
            ylabel('Longitude (degrees)');
        end
        
        function plotHeading(obj)
            plot(obj.time, obj.heading);
            title('Heading vs. Time');
            xlabel('Time (s)');
            ylabel('Heading (rad)');
        end
        
        function plotBankAngle(obj)
            hold on;
            plot(obj.time, obj.bankAngle);
            plot(obj.time, obj.commandedBankAngle, '--');
            hold off;
            title('Bank Angle vs. Time');
            xlabel('Time (s)');
            ylabel('Bank Angle (rad)');
            legend('Bank Angle', 'Commanded Bank Angle');
        end
        
        function plotLiftDrag(obj)
            hold on;
            plot(obj.time, obj.drag);
            plot(obj.time, obj.lift);
            plot(obj.time, sqrt(obj.lift.^2 + obj.drag.^2));
            hold off;
            title('Normalized Lift and Drag vs. Range');
            xlabel('Normalized Range to Go');
            ylabel('Lift and Drag Accelerations (g)');
            legend('Lift', 'Drag', 'Total Load');
        end
        
        function plotHeatingRate(obj)
            plot(obj.time(2:end), obj.heatingRate(2:end));
            title('Heating Rate History');
            xlabel('Normalized Range to Go');
            ylabel('Heating Rate ($W/cm^2$)');
        end
        
        function plotVelocityAltitude(obj)
            plot(obj.velocityUnnormalized, obj.altitudeUnnormalized./1e3);
            title('Velocity vs. Altitude');
            xlabel('Velocity (m/s)');
            ylabel('Altitude (km)');
        end
        
        function plotECEFPosition(obj)
            % hold on;
            % plot(obj.time, obj.posECEF(1, :), 'DisplayName', 'X');
            % plot(obj.time, obj.posECEF(2, :), 'DisplayName', 'Y');
            % plot(obj.time, obj.posECEF(3, :), 'DisplayName', 'Z');
            % hold off;
            % legend show;
            % title('ECEF Position Components vs. Time');
            % xlabel('Time (s)');
            % ylabel('Position (m)');
        end
    end
end

