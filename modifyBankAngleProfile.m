function bankAngle = modifyBankAngleProfile(rangeError, bankAngleHistory, cntIterations)
    % Since we are starting the iteration count at 1 instead of zero (damn you MATLAB)
    % we must subtract 1 from the iteration count in our lambda calculation
    lambda = 0.5^cntIterations;

    % If the range error is greater than our allowable tolerance, we need
    % to update the bank angle and try again.
    bankAngle = bankAngleHistory(cntIterations) - lambda * (rangeError(cntIterations) / (rangeError(cntIterations) - rangeError(cntIterations - 1))) * (bankAngleHistory(cntIterations) - bankAngleHistory(cntIterations - 1));


end
