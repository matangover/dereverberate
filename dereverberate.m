%% Dereverberation

inputFilename = 'EchoSample.mp3';
% inputFilename = 'stalbans_omni_sing.mp3';
% inputFilename = 'mozart_reverb_short.mp3';
[signal, fs] = audioread(inputFilename);
% Take only the left channel from the stereo recording
signal = signal(:, 1);
inputLength = length(signal);
%%
% Compute the signal's STFT with nfft=1024
% and a 50% overlapping Hann window.
stftWindowSize = 1024;
overlap = 0.75; % In the original paper they used 0.5, but this makes artifacts in reconstruction. Why?
overlapSamples = overlap*stftWindowSize;
% TODO: zero padding?
window = hann(stftWindowSize);
[M,w,t] = spectrogram(signal, window, overlapSamples); %, fs)
spectrogram(signal, hann(stftWindowSize), overlap*stftWindowSize, 'yaxis');
title("Input - before processing");
% M ('microphone') - input signal (reverberated) frequency-domain vectors.
% Rows - frequencies, Columns - time frames (middle time-point)
[frequencyCount, frameCount] = size(M);
%%
% Constants

% Impulse response block length - has to be 'sufficiently small'
% TODO: experiment with value. Is this really the same as the window size?
D = stftWindowSize;
B = 50; % Number of impulse response blocks - 1. TODO: experiment.
% Minimum gain per frequency. TODO: Should be frequency specific? Experiment.
minGain = zeros(frequencyCount, 1);
bias = zeros(frequencyCount, B);
bias(:) = 1.05; % TODO: Experiment
maxHEstimate = zeros(frequencyCount, B);
% TODO: Experiment - should 'reflect real world systems' -
% e.g. exponential decay / given impulse response. (Page 11 - MaxValue.)
maxHEstimate(:) = 0.9;
%%
S = zeros(size(M)); % Dry signal frequency-domain vectors.
R = zeros(size(M)); % Reverberated components frequency-domain vectors.

H_pow = zeros(frequencyCount, B); % Reverberant system frequency response estimate blocks - power.
pow = @(x) abs(x).^2;

previousSFramesPower = zeros(frequencyCount, B); % Most recent previous frame is first.
previousMFramesPower = ones(frequencyCount, B); % Most recent previous frame is first.

C = zeros(size(H_pow)); % Reverberant system frequency response power estimate blocks - temporary.
H_pow_all = zeros([frameCount size(H_pow)]);

for inputFrameIndex = 1:frameCount
    inputFrame = M(:, inputFrameIndex);
    inputFramePower = pow(inputFrame);
    
    % Update system frequency response estimates.
    for blockIndex=1:B
        newEstimates = inputFramePower ./ previousMFramesPower(:, blockIndex);
        unchanged = newEstimates >= H_pow(:, blockIndex);
        newEstimates(unchanged) = H_pow(unchanged, blockIndex) .* bias(unchanged, blockIndex) + eps;
        C(:, blockIndex) = min(newEstimates, maxHEstimate(:, blockIndex));
        % TODO: Temporal smoothing (Page 11).
        H_pow(:, blockIndex) = C(:, blockIndex);
    end
    H_pow_all(inputFrameIndex, :, :) = H_pow;
    
    G_S = 1 - (sum(previousSFramesPower .* H_pow)) ./ inputFramePower;
    %G_S = sqrt(G_S);
    % Enforce a minimum gain for each frequency.
    G_S = max([G_S minGain], [], 2);
    G_R = 1 - G_S;
    % TODO: Apply gain smooting over time to G_S and G_R.
    S(:, inputFrameIndex) = G_S .* inputFrame;
    R(:, inputFrameIndex) = G_R .* inputFrame;
    
    % Shift all previous frames one column to the right, and insert the
    % current frame at the beginning.
    previousSFramesPower(:, 2:end) = previousSFramesPower(:, 1:end-1);
    previousSFramesPower(:, 1) = pow(S(:, inputFrameIndex));
    
    previousMFramesPower(:, 2:end) = previousMFramesPower(:, 1:end-1);
    previousMFramesPower(:, 1) = inputFramePower;
end
%%
% Display the spectrum power before and after:
spectrogramPlot(M, t, w); % TODO: M.^2
title("Input");
spectrogramPlot(S, t, w);
title("Dry");
spectrogramPlot(R, t, w);
title("Reverberant components");

spectrogramPlot(H_pow, t, w);
title("Final frequency response estimate");
%%
dry = reconstruct(S, window, overlapSamples, inputLength);
soundsc(dry, fs);
spectrogram(dry, window, overlapSamples, 'yaxis');
title("Dry - reconstructed")
%%
reverberant = reconstruct(R, window, overlapSamples, inputLength);
soundsc(reverberant, fs);
spectrogram(reverberant, window, overlapSamples, 'yaxis');
title("Reverberant components - reconstructed");
%%
reverbConstant = 1;
%reverbConstant = 2;
reverberated = reverbConstant*reverberant + dry;
soundsc(reverberated, fs);
spectrogram(reverberated, window, overlapSamples, 'yaxis');
title("Reverberated signal - reconstructed")
%%
for i=[1:200:frameCount frameCount]
    spectrogramPlot(squeeze(H_pow_all(i, :, :)), t, w);
    title("Frequency response estimate - at input frame " + i);
end
%%
% De-reverberate using
% Without block-wise magic - simple frequency domain processing with inverse filter.
% (Divide by filter i
%%
% De-reverberate again, but this time use the final frequency response estimate from the beginning.

S = zeros(size(M)); % Dry signal frequency-domain vectors.
R = zeros(size(M)); % Reverberated components frequency-domain vectors.
% Take the last one frequency response estimate.
H_pow = squeeze(H_pow_all(end, :, :));

previousSFramesPower = ones(frequencyCount, B); % Most recent previous frame is first.
for inputFrameIndex = 1:frameCount
    inputFrame = M(:, inputFrameIndex);
    inputFramePower = pow(inputFrame);
        
    G_S = 1 - (sum(previousSFramesPower .* H_pow)) ./ inputFramePower;
    %G_S = sqrt(G_S);
    % Enforce a minimum gain for each frequency.
    G_S = max([G_S minGain], [], 2);
    G_R = 1 - G_S;
    % TODO: Apply gain smooting over time to G_S and G_R.
    S(:, inputFrameIndex) = G_S .* inputFrame;
    R(:, inputFrameIndex) = G_R .* inputFrame;
    
    % Shift all previous frames one column to the right, and insert the
    % current frame at the beginning.
    previousSFramesPower(:, 2:end) = previousSFramesPower(:, 1:end-1);
    previousSFramesPower(:, 1) = pow(S(:, inputFrameIndex));
end
%%
% Dereverberate using inverse filtering
finalEstimate = squeeze(H_pow_all(end, :, :)); % TODO: sqrt?
firstBlock = ones(frequencyCount, 1);
finalEstimate = [firstBlock finalEstimate];
stftHop = stftWindowSize - overlapSamples;
blockCount = B+1;
irLength = B * stftHop + overlapSamples;
impulseResponse = reconstruct(finalEstimate, rectwin(stftWindowSize), overlapSamples, irLength);
%%
plot(impulseResponse);
spectrogramPlot(finalEstimate, t, w);
figure;
freqz(impulseResponse);

%freqz(1, impulseResponse);
%impz(1, impulseResponse);

%%
dryInverseFiltered = filter(1, impulseResponse, signal);
soundsc(dryInverseFiltered, fs);
figure;
% The output is exponentially increasing - filter is unstable.
spectrogram(dryInverseFiltered, window, overlapSamples, 'yaxis');
title("Dry - inverse filtered");

%%
% Scale the filter so that all zeros are strictly inside the unit circle.
[z,p,k] = tf2zpk(impulseResponse);
len = length(z);
km = k;
z_minP = zeros(len,1);
for i = 1:length(z)
    % if magnitude of zero greater than 1, replace the zero with its
    % inverse conjugate and scale the gain parameter.
    if(zMag(i) > 1)
        z_minP(i) = 1/conj(z(i));
        km = k*conj(z(i));
    else
        z_minP(i) = z(i);
    end
end

scaledImpulseResponse = zpk2tf(z_minP, p, newK);
% Now the inverse filter's poles are guaranteed to be in the unit circle, so it's stable.
%%
% Draw movie frames (very slow - run in Matlab rather than Live Editor).
F(frameCount) = struct('cdata',[],'colormap',[]);
image = spectrogramPlot(squeeze(H_pow_all(1, :, :)), t, w);
for i = 1:frameCount
    image.CData = squeeze(H_pow_all(i, :, :));
    drawnow();
    F(i) = getframe(gcf);
end
%%
% Display movie interactively (slow).
figure;
movie(F);
%%
% Write movie to file.
video = VideoWriter('reverb.avi');
open(video);
writeVideo(video, F);
close(video)

%%
function output = reconstruct(spectrum, window, overlapSamples, outputLength)
    stftWindowSize = length(window);
    [frequencyCount, frameCount] = size(spectrum);
    spectrum(frequencyCount+1:stftWindowSize, :) = conj(flipud(spectrum(1:frequencyCount-2, :)));
    output = zeros(outputLength, 1);
    stftHop = stftWindowSize - overlapSamples;
    for frameIndex = 1:frameCount
      frameStart = (frameIndex-1)*stftHop;
      if frameStart + stftWindowSize > outputLength
          break
      end
      sampleRange = frameStart+1 : frameStart+stftWindowSize;
      reconstructedFrame = window .* ifft(spectrum(:, frameIndex), 'symmetric');
      output(sampleRange) = output(sampleRange) + reconstructedFrame;
    end
end

function image = spectrogramPlot(spectrum, t, w)
    figure;
    image = imagesc(t, w, 10*log10(abs(spectrum)+eps));
    image.Parent.YDir = 'normal';
    colorbar();
end