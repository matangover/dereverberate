inputFilename = 'singing.mp3';
irFilename = 'stalbans_a_mono.wav';

[signal, fs] = audioread(inputFilename);
[impulseResponse, irFs] = audioread(irFilename);

assert(fs == irFs);
impulseResponse = impulseResponse(1:50000);
figure;
plot(impulseResponse);
%zeroPaddedInput = [signal; zeros(length(impulseResponse), 1)];
reverberated = conv(signal, impulseResponse);

%%

paddedSignal = [signal; zeros(length(impulseResponse) - 1, 1)];
%subplot(2,1,1);
figure;
hold on;
plot(reverberated);
plot(paddedSignal);
hold off;
legend('Reverberated', 'Original (dry)');

soundsc([signal; reverberated], fs);

%% inverse filtering - frequency domain
nfft = length(reverberated);
hf = fft(impulseResponse, nfft);
spectrum = fft(reverberated);

spectrum = spectrum ./ hf;
% To reduce noise amplification:
% spectrum = spectrum .* conj(hf) ./ (abs(hf).^2 + 1e-2);
% See: https://blogs.mathworks.com/steve/2007/08/13/image-deblurring-introduction/
inverseFiltered = real(ifft(spectrum));
soundsc(inverseFiltered, fs);

%%
figure;
freqz(reverberated);
title('Reverberated');

figure;
freqz(inverseFiltered);
title('Inverse filtered');

figure;
freqz(impulseResponse);
title('Impulse response');

%% inverse filtering - time domain
% doesn't work
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
isstable(1, impulseResponse)
figure;
plot(inverseFiltered);