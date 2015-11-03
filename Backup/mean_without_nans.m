function [m] = mean_without_nans(x, w)
if nargin == 1
    w = ones(size(x));
end
assert(isvector(x));
x = x .* w;
xx = x(~isnan(x));
if isempty(xx)
    m = nan;
else
    m = sum(xx);
end