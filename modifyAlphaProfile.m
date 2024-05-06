function angleOfAttack = modifyAlphaProfile(rangeErrorHistory, angleOfAttackHistory, cntIterations)
    if cntIterations == 0
        % Newton-Raphson method for the first iteration
        lambda = 0.5^cntIterations;
        dz_dsigma = (rangeErrorHistory(cntIterations + 2) - rangeErrorHistory(cntIterations + 1)) / (angleOfAttackHistory(2) - angleOfAttackHistory(1));
        angleOfAttack = angleOfAttackHistory(cntIterations + 1) - lambda * (rangeErrorHistory(cntIterations + 1) / dz_dsigma);
    else
        % Secant method for subsequent iterations
        lambda = 0.5^(cntIterations - 1);
        angleOfAttack = angleOfAttackHistory(cntIterations + 1) - lambda * (rangeErrorHistory(cntIterations + 1) / ((rangeErrorHistory(cntIterations + 1) - rangeErrorHistory(cntIterations)) / (angleOfAttackHistory(cntIterations + 1) - angleOfAttackHistory(cntIterations))));
    end
end
