inputFilename = 'singing.mp3';
irFilename = 'stalbans_a_mono.wav';

[signal, fs] = audioread(inputFilename);
[impulseResponse, irFs] = audioread(irFilename);

assert(fs == irFs);
impulseResponse = impulseResponse(1:50000);
figure;
plot(impulseResponse);
zeroPaddedInput = [signal; zeros(length(impulseResponse), 1)];
reverberated = conv(zeroPaddedInput, impulseResponse);

%%
figure;
plot(signal);
figure;
plot(reverberated);
soundsc([signal;reverberated], fs);

%%
[dry, r] = deconv(reverberated, impulseResponse);

%%
figure;
plot(dry);
soundsc(dry, fs);
% Why does this not work?
% In the beginning, it resembles an periodic signal,
% but explodes very fast.

%%
figure;
freqz(impulseResponse);
% figure;
% zplane(impulseResponse);

%%
inverseFiltered = filter(1, impulseResponse, reverberated);
%%
freqz(1, impulseResponse);
isstable(1, impulseResponse)
figure;
plot(inverseFiltered);