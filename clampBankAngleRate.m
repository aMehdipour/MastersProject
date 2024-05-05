function bankAngleRateClamped = clampBankAngleRate(bankAngle, bankAngleCommand, constants, dt)

    % Clamp the bank angle rate to respect the limits
    bankAngleRateLimit   = constants.BANK_ANGLE_RATE_LIMIT * dt;
    bankAngleRateRaw     = bankAngleCommand - bankAngle;
    bankAngleRateClamped = max([min([bankAngleRateRaw, bankAngleRateLimit]), -bankAngleRateLimit]);
