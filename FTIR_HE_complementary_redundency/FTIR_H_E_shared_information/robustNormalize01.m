function x = robustNormalize01(x)

x = double(x(:));

%% =========================================================
% Robust limits
%% =========================================================

lowValue  = prctile(x,1);
highValue = prctile(x,99);


%% =========================================================
% Protect against constant images
%% =========================================================

if highValue <= lowValue

    x = zeros(size(x));
    return

end


%% =========================================================
% Clip outliers
%% =========================================================

x(x < lowValue) = lowValue;

x(x > highValue) = highValue;


%% =========================================================
% Normalize to 0-1
%% =========================================================

x = (x-lowValue) / ...
    (highValue-lowValue);

end