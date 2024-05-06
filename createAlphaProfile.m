function angleOfAttackProfile = createAlphaProfile(initialAngleOfAttack, finalAngleOfAttack, energy)

    % Calculate the bank angle based on the bank angle magnitude profile
    % (Eq. 24 in Low Lifting Entry Guidance)
    angleOfAttackProfile = initialAngleOfAttack + (finalAngleOfAttack - initialAngleOfAttack) / (energy(end) - energy(1)) .* (energy - energy(1));
end
