%% Dereverberation

inputFilename = 'audio/EchoSample.mp3';
% inputFilename = 'audio/stalbans_omni_sing.mp3';
% inputFilename = 'audio/mozart_reverb_short.mp3';
[signal, fs] = audioread(inputFilename);
% Take only the left channel from the stereo recording
signal = signal(:, 1);
inputLength = length(signal);

%% Compute the signal's STFT
stftWindowSize = 1024;
overlap = 0.75; % In the paper he used 0.5 (50% overlap).
overlapSamples = overlap*stftWindowSize;
% TODO: zero padding?
window = hann(stftWindowSize, 'periodic');
% Specify empty Fs so that t is returned in samples and w in rad/sample.
[M,w,t] = spectrogram(signal, window, overlapSamples, [], []); 
spectrogram(signal, window, overlapSamples, 'yaxis');
title("Input - before processing");
% M ('microphone') - input signal (reverberated) frequency-domain vectors.
% Rows - frequencies, Columns - time frames (middle time-point)
[frequencyCount, frameCount] = size(M);

%% Constants
% Impulse response block length - has to be 'sufficiently small'
% TODO: experiment with value. Is this really the same as the window size?
D = stftWindowSize - overlapSamples;
B = 400; % Number of impulse response blocks. TODO: experiment.

% Minimum gain per frequency. TODO: Experiment.
minGain = zeros(frequencyCount, 1);
maxHEstimate = zeros(frequencyCount, B);

% Maximum magnitude estimate per impulse response block and frequency.
% Should 'reflect real world systems'. (Page 11 - MaxValue.)
maxHEstimate(:) = readImpulseResponse('audio/stalbans_a_mono.wav', B, window, overlapSamples);
%maxHEstimate(:) = 0.9;

% Bias used to keep the magnitude estimate from getting stuck on a wrong
% minimum.
bias = zeros(frequencyCount, B);
bias(:) = 1.01; % TODO: Experiment

% From the paper - gamma (Page 10). Lower means more smoothing between
% gains vectors of consecutive frames. Value 0-1 (1 = no smoothing).
gainSmoothingFactor = zeros(frequencyCount, 1);
gainSmoothingFactor(:) = 0.3;

% From the paper - alpha (Page 11). Lower means less smoothing between
% magnitude response block estimates on consecutive frames.
% Value 0-1 (0 = no smoothing).
magnitudeSmoothingFactor = zeros(frequencyCount, B);
magnitudeSmoothingFactor(:) = 0.2;

%% Algorithm implementation: block-wise signal dereverberation.
S = zeros(size(M)); % Dry signal frequency-domain vectors.
R = zeros(size(M)); % Reverberated components frequency-domain vectors.

% Reverberant system frequency response estimate blocks - power.
% H_pow = zeros(frequencyCount, B); 
% Set initial estimates to half of the maximum.
H_pow = maxHEstimate / 2;
pow = @(x) abs(x).^2;

previousSFramesPower = zeros(frequencyCount, B); % Most recent previous frame is first.
previousMFramesPower = ones(frequencyCount, B); % Most recent previous frame is first.

% Reverberant system frequency response power estimate blocks - temporary.
C = zeros(size(H_pow)); 
% Matrix to keep all estimates over time (not needed by algorithm - just
% for visualization purposes).
H_pow_all = zeros([frameCount size(H_pow)]);
% Initially, all components are estimated to belong to the dry signal.
G_S = ones(frequencyCount, 1);
G_R = zeros(frequencyCount, 1);

for inputFrameIndex = 1:frameCount
    inputFrame = M(:, inputFrameIndex);
    inputFramePower = pow(inputFrame);
    
    % Update system frequency response estimates.
    for blockIndex=1:B
        newEstimates = inputFramePower ./ previousMFramesPower(:, blockIndex);
        unchanged = newEstimates >= H_pow(:, blockIndex);
        newEstimates(unchanged) = H_pow(unchanged, blockIndex) .* bias(unchanged, blockIndex) + eps;
        C(:, blockIndex) = min(newEstimates, maxHEstimate(:, blockIndex));
        alpha = magnitudeSmoothingFactor(:, blockIndex);
        H_pow(:, blockIndex) = alpha .* H_pow(:, blockIndex) + (1-alpha) .* C(:, blockIndex);
    end
    H_pow_all(inputFrameIndex, :, :) = H_pow;
    
    newG_S = 1 - (sum(previousSFramesPower .* H_pow, 2)) ./ inputFramePower;
    % TODO: ?? G_S = sqrt(G_S);
    % Enforce a minimum gain for each frequency.
    newG_S = max([newG_S minGain], [], 2);
    G_S = (1-gainSmoothingFactor).*G_S + gainSmoothingFactor.*newG_S;
    newG_R = 1 - G_S;
    G_R = (1-gainSmoothingFactor).*G_R + gainSmoothingFactor.*newG_R;
    S(:, inputFrameIndex) = G_S .* inputFrame;
    R(:, inputFrameIndex) = G_R .* inputFrame;
    
    % Shift all previous frames one column to the right, and insert the
    % current frame at the beginning.
    previousSFramesPower(:, 2:end) = previousSFramesPower(:, 1:end-1);
    previousSFramesPower(:, 1) = pow(S(:, inputFrameIndex));
    
    previousMFramesPower(:, 2:end) = previousMFramesPower(:, 1:end-1);
    previousMFramesPower(:, 1) = inputFramePower;
end

%% Visualize algorithm results
% Display the spectrum power before and after:
spectrogramPlot(M, t, w);
title("Input");
spectrogramPlot(S, t, w);
title("Dry");
spectrogramPlot(R, t, w);
title("Reverberant components");

spectrogramPlot(H_pow, t, w);
title("Final frequency response estimate");

%% Reconstruct and play dry signal
dry = reconstruct(S, window, overlapSamples, inputLength);
soundsc([signal; dry], fs);
figure;
spectrogram(dry, window, overlapSamples, 'yaxis');
title("Dry - reconstructed")

%% Reconstruct and play reverberant components
reverberant = reconstruct(R, window, overlapSamples, inputLength);
soundsc(reverberant, fs);
spectrogram(reverberant, window, overlapSamples, 'yaxis');
title("Reverberant components - reconstructed");

%% Reconstruct and play original reverberant signal using both components
reverbConstant = 1;
%reverbConstant = 5;
reverberated = reverbConstant*reverberant + dry;
soundsc(reverberated, fs);
spectrogram(reverberated, window, overlapSamples, 'yaxis');
title("Reverberated signal - reconstructed")

%% Plot some interim magnitude response blocks.
for i=[1:200:frameCount frameCount]
    spectrogramPlot(squeeze(H_pow_all(i, :, :)), t, w);
    title("Frequency response estimate - at input frame " + i);
end

%% Dereverberate using final frequency response estimate
% De-reverberate again, but this time use the final frequency response
% estimate right from the beginning.

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

%% Dereverberate using inverse filtering (unstable).
finalEstimate = squeeze(H_pow_all(end, :, :)); % TODO: sqrt?
firstBlock = ones(frequencyCount, 1);
finalEstimate = [firstBlock finalEstimate];
stftHop = stftWindowSize - overlapSamples;
blockCount = B+1;
irLength = B * stftHop + overlapSamples;
impulseResponse = reconstruct(finalEstimate, rectwin(stftWindowSize), overlapSamples, irLength);

plot(impulseResponse);
spectrogramPlot(finalEstimate, t, w);
figure;
freqz(impulseResponse);

%freqz(1, impulseResponse);
%impz(1, impulseResponse);

% This doesn't really work. The output is exponentially increasing - the
% inverse filter is unstable.
dryInverseFiltered = filter(1, impulseResponse, signal);
soundsc(dryInverseFiltered, fs);
figure;
spectrogram(dryInverseFiltered, window, overlapSamples, 'yaxis');
title("Dry - inverse filtered");

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
% Now the inverse filter's poles are guaranteed to be in the unit circle,
% so it's stable.

%% Draw movie frames
clear F;
frameHop = 5;
movieFrameCount = floor(frameCount / frameHop);
F(movieFrameCount) = struct('cdata',[],'colormap',[]);
stftHop = stftWindowSize - overlapSamples;
t2 = stftHop*(1:B);
image = spectrogramPlot(squeeze(H_pow_all(1, :, :)), t2, w);
for i = 1:frameHop:frameCount
    image.CData = squeeze(pow2db(H_pow_all(i, :, :)));
    drawnow();
    F((i - 1)/ frameHop + 1) = getframe(gcf);
end

%% Display movie interactively.
figure;
movie(F);

%% Write movie to file.
video = VideoWriter('estimate_evolution.avi');
open(video);
writeVideo(video, F);
close(video)

%% Compute length of estimated impulse response
stftHop = stftWindowSize - overlapSamples;
impulseResponseLengthSamples = stftHop * (B+1) + overlapSamples;
impulseResponseLengthSeconds = impulseResponseLengthSamples / fs;
impulseResponseLengthSeconds

%% Write results to audio files
audiowrite('results/dry.wav', dry ./ (max(abs(dry))), fs);
audiowrite('results/reverberant.wav', reverberant ./ (max(abs(reverberant))), fs);
audiowrite('results/reverberated.wav', reverberated ./ (max(abs(reverberated))), fs);

%% Testing spectrogramPlot
spectrogramPlot(M, t, w);
figure;
spectrogram(signal, window, overlapSamples, 'yaxis', 'power');

%% Debug impulse response read
irDebug = readImpulseResponse('audio/stalbans_a_mono.wav', B, window, overlapSamples);
spectrogramPlot(irDebug, t, w);
spectrogramPlot(H_pow, t, w);

%% Low-pass filter reverberant components
% Doesn't help.
reverberantFiltered = filter(lowpassfilter, reverberant);
spectrogram(reverberantFiltered, window, overlapSamples, 'yaxis');
soundsc(reverberantFiltered, fs);

%% Helper functions

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

function blocks = readImpulseResponse(irFilename, blockCount, window, overlapSamples)

    impulseResponse = audioread(irFilename);
    %assert(irFs == fs);
    allBlocks = spectrogram(impulseResponse, window, overlapSamples);
    pow = @(x) abs(x).^2;
    blocks = pow(allBlocks(:, 1:blockCount));

end