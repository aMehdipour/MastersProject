function state = createStateObject(parameters, constants)

    maxIters              = parameters.tMax / parameters.dt;
    state.flightPathAngle = NaN(maxIters,1);
    state.rangeToGo       = NaN(maxIters,1);
    state.altitude        = NaN(maxIters,1);
    state.velocity        = NaN(maxIters,1);
    state.r               = NaN(maxIters,1);
    state.v               = NaN(maxIters,1);
    state.lat             = NaN(maxIters,1);
    state.lon             = NaN(maxIters,1);
    state.crossRange      = NaN(maxIters,1);
    state.downRange       = NaN(maxIters,1);
    state.energy          = NaN(maxIters,1);
    state.posECEF         = NaN(3,maxIters);
    state.posNED          = NaN(3,maxIters);
    state.rotNedFromEcef  = NaN(3,3,maxIters);
    state.rotBodyFromEcef = NaN(3,3,maxIters);

