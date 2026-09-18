function R = rowCorrelation(A, B)

% Each row = one spectrum

A = A - mean(A, 2);
B = B - mean(B, 2);

numerator = sum(A .* B, 2);

denominator = sqrt( ...
    sum(A.^2, 2) .* sum(B.^2, 2));

R = numerator ./ (denominator + eps);

end