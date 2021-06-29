function r = calc_corr_val(x,y)
% this is from Emery & Thompson p 243
% assume x & y are vectors of the same length
N = length(x);

mean_x = nanmean(x);
mean_y = nanmean(y);

std_x = nanstd(x);
std_y = nanstd(y);

covar_xy  = (1 ./ (N - 1)) .* nansum((x - mean_x).*(y - mean_y));

r = covar_xy ./ (std_x * std_y);

end