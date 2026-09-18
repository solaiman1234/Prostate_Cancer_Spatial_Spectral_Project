function value = stableLogDet(C)

C = (C + C')/2;

[R,p] = chol(C);

if p ~= 0

    C = C + 1e-6*eye(size(C));

    R = chol(C);

end

value = ...
    2*sum(log(diag(R)));

end