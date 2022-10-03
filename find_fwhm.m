function fwhm = find_fwhm(vec, focus)
% Function to find the FWHM of a signal based on the power in the focal
% point and the minimal power beforehand (the minimal power is taken before
% the focus to overcome the atteuation in calculation of axial FWHM)
% Input:
%   vec = numeric signal
%   focus = the index of the focus to calculate FWHM around
% Output:
%   fwhm = signal FWHM in vector length units
if nargin < 2
    focus = round(length(vec)/2);
end
v_half = (vec(focus)-min(vec(1:focus)))/2 + min(vec(1:focus));
[~,I_min] = min(vec(1:focus));
id1 = find(vec(I_min:end)>=v_half,1) + I_min;
id2 = find(vec(I_min:end)>=v_half,1,'last') + I_min;
fwhm = id2-id1;
end
