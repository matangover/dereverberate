%% Load the input file and impulse response.
inputFilename = 'audio/singing.mp3';
% inputFilename = 'mozart.mp3';
irFilename = 'audio/stalbans_a_mono.wav';

[signal, fs] = audioread(inputFilename);
[impulseResponse, irFs] = audioread(irFilename);
assert(fs == irFs);

% Truncate for faster processing.
impulseResponse = impulseResponse(1:50000);

figure;
plot(impulseResponse);

%% Apply reverb using convolution.
reverberated = conv(signal, impulseResponse);

%% Plot and play the original and reverberated sound.

paddedSignal = [signal; zeros(length(impulseResponse) - 1, 1)];
%subplot(2,1,1);
figure;
hold on;
plot(reverberated);
plot(paddedSignal);
hold off;
legend('Reverberated', 'Original (dry)');

soundsc([signal; reverberated], fs);

%% Inverse filtering - frequency domain
nfft = length(reverberated);
hf = fft(impulseResponse, nfft);
spectrum = fft(reverberated);

spectrum = spectrum ./ hf;
% To reduce noise amplification:
% spectrum = spectrum .* conj(hf) ./ (abs(hf).^2 + 1e-2);
% See: https://blogs.mathworks.com/steve/2007/08/13/image-deblurring-introduction/
inverseFiltered = real(ifft(spectrum));
soundsc(inverseFiltered, fs);

%% Plot the signals frequency responses
figure;
freqz(reverberated);
title('Reverberated');

figure;
freqz(inverseFiltered);
title('Inverse filtered');

figure;
freqz(impulseResponse);
title('Impulse response');

%% Inverse filtering - time domain
% Doesn't work.
[dry, r] = deconv(reverberated, impulseResponse);

figure;
plot(dry);
soundsc(dry, fs);
% Why does this not work?
% In the beginning, it resembles an periodic signal,
% but explodes very fast.

figure;
freqz(impulseResponse);
% figure;
% zplane(impulseResponse);

inverseFiltered = filter(1, impulseResponse, reverberated);
freqz(1, impulseResponse);
% The inverse filter is not stable because the original filter
% is not minimal phase. The original filter has zeros on or outside of the
% unit circle.
isstable(1, impulseResponse)
figure;
plot(inverseFiltered);