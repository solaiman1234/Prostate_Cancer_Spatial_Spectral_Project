function Z = gaussianize(X)

[n,p] = size(X);

Z = zeros(n,p);


for j = 1:p

    r = tiedrank(X(:,j));

    u = (r - 0.5) / n;

    %
    % Prevent exactly 0 or 1
    %

    u = max(u,1e-6);

    u = min(u,1-1e-6);


    Z(:,j) = norminv(u);

end

end