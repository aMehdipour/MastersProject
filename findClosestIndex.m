function idx = findClosestIndex(arr, value)
    [~, idx] = min(abs(arr - value));
end

