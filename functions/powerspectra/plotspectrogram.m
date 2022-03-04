function [ftF, ftT, ftP] = plotspectrogram(signal,fs,winlen,overlap,foi)
%PLOTSPECTROGRAM Plot the spectrogram of a signal.
%   Use as:
%       [ftF, ftT, ftP] = plotspectrogram(data,fs,winlen,overlap,foi);
%   Input:
%       - signal, the signal vector, each column is one channel
%       - fs, sampling rate of the siganl
%       - winlen, the window length of FFT
%       - overlap, overlap*winlen will be used as the length of overlap
%       window
%       - foi, the frequency band of interest
%   Output:
%       - ftF, the frequency vector of the spectrogram
%       - ftT, the time vector of the spectrogram
%       - ftP, the power(dB) matrix
%
%   Author   : NIE Yingnan
%   Created  : Aug 9, 2020
%   Modified : Aug 9, 2020

noverlap = ceil(winlen*overlap);

[~, ftF, ftT, ftP] = spectrogram(signal, winlen, noverlap, [], fs);

surf(ftT, ftF, log10(ftP), 'edgecolor', 'none');
view(0, 90);
axis tight;
ylim(foi);
grid off; 
set(gca, 'LineWidth', 2, 'FontSize', 10, 'xtick', [], 'xcolor', 'w');
colormap('jet');
colorbar('location', 'EastOutside', 'FontSize', 10);
set(gca, 'Clim', [-2, 0.5]);

end

