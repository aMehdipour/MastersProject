function [density, temperature, pressure] = atmosphereModel(altitude)
    if altitude >= 25000
        temperature = -131.21 + 0.00299 * altitude;
        pressure = 2.488 * ((temperature + 273.1) / 216.6)^(-11.388);
    elseif 11000 <= altitude && altitude < 25000
        temperature = -56.46;
        pressure = 22.65 * exp(1.73 - 0.000157 * altitude);
    elseif altitude < 11000
        temperature = 15 - 0.00649 * altitude;
        pressure = 101.29 * ((temperature + 273.1) / 288.08)^5.256;
    end

    temperature = temperature + 273.1;
    density = pressure / (0.2869 * temperature);

    % Allow complex input for use of complex-step derivative in finding beta_r
    if ~isreal(altitude)
        [rho, pressure, temperature] = atmosphereModel(real(altitude));
        rho = complex(rho, 0);
        pressure = complex(pressure, 0);
        temperature = complex(temperature, 0);
    end
end
