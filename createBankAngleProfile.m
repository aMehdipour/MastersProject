function bankAngleProfile = createBankAngleProfile(initialBankAngle, finalBankAngle, energy)

    % Calculate the bank angle based on the bank angle magnitude profile
    % (Eq. 24 in Low Lifting Entry Guidance)
    bankAngleProfile = initialBankAngle + (finalBankAngle - initialBankAngle) / (finalEnergy - initialEnergy) .* (energy - initialEnergy);
end
