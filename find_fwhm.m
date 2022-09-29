function fwhm = find_fwhm(vec, I_max)
% Function to find the FWHM of a signal using linear interpolation
% Input:
%   vec = numeric signal
%   I_max = (optional) index of the maximum in case its a local maxima and
%           not a global maxima
% Output:
%   fwhm = signal FWHM in vector length units
if nargin < 2
    [v_max,I_max] = max(vec);
else
    v_max = vec(I_max);
end
[~, I_half] = min(abs(vec-0.5*v_max));
if length(I_half) > 1
    d_from_max = min(abs(I_max-I_half));
    fwhm = 2*d_from_max;
else
    fwhm = 2*abs(I_max-I_half);
end
