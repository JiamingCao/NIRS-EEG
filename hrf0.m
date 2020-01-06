function hrf = hrf0(length, Fs, t0, beta)
alpha = 6;
t = 0:1/Fs:length;
hrf = ((t-t0).^(alpha-1).*beta.^alpha.*exp(-beta*(t-t0)))/gamma(alpha);
