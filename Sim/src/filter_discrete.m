function Y = filter_discrete(In, Mem,a,b)
% Implements a discretized second-order filter

Y = (1/a(1)) * (b(1) * In(1) + b(2) * In(2) + b(3) * In(3)...
    - a(2) * Mem(1) - a(3) * Mem(2));