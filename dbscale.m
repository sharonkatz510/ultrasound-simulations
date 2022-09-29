function scaled_im = dbscale(img, lim)
% Function to scale image in the range of 0-(-lim) db
num_im = im2double(img);
a = (1-10^(lim/(-20)))/(max(max(num_im))-min(min(num_im)));
b = 1-a*max(max(num_im));
scaled_im = 20.*log10(a.*num_im + b);
end