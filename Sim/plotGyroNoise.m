function plotGyroNoise(data)

    numNoise = length(data.gyroNoiseDensity);

    %=============================================
    %% Plot Heading Angle Cmd and Response
    %=============================================
    figure
    tiledlayout(2,1)
    nexttile
    hold on
    for i = 1:numNoise
        plot(data.time(2:end-1), data.flightPathCmd(:,:,i))
        plot(data.time, data.gamma(:,i))
    end



