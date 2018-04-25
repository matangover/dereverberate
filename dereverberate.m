%% Dereverberation

inputFilename = 'audio/EchoSample.mp3';
% inputFilename = 'stalbans_omni_sing.mp3';
% inputFilename = 'mozart_reverb_short.mp3';
[signal, fs] = audioread(inputFilename);
% Take only the left channel from the stereo recording
signal = signal(:, 1);
inputLength = length(signal);
%%
% Compute the signal's STFT with nfft=1024
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
%%
% Constants

% Impulse response block length - has to be 'sufficiently small'
% TODO: experiment with value. Is this really the same as the window size?
D = stftWindowSize - overlapSamples;
B = 400; % Number of impulse response blocks - 1. TODO: experiment.

% Minimum gain per frequency. TODO: Should be frequency specific? Experiment.
minGain = zeros(frequencyCount, 1);
bias = zeros(frequencyCount, B);
bias(:) = 1.05; % TODO: Experiment
maxHEstimate = zeros(frequencyCount, B);
% TODO: Experiment - should 'reflect real world systems' -
% e.g. exponential decay / given impulse response. (Page 11 - MaxValue.)
maxHEstimate(:) = readImpulseResponse('audio/stalbans_a_mono.wav', B, window, overlapSamples);
%maxHEstimate(:) = 0.9;
% From the paper - gamma (Page 10). Lower means more smoothing between
% gains vectors of consecutive frames. Value 0-1 (1 = no smoothing).
gainSmoothingFactor = zeros(frequencyCount, 1);
gainSmoothingFactor(:) = 0.3;

% From the paper - alpha (Page 11). Lower means less smoothing between
% magnitude response block estimates on consecutive frames.
% Value 0-1 (0 = no smoothing).
magnitudeSmoothingFactor = zeros(frequencyCount, B);
magnitudeSmoothingFactor(:) = 0.2;

%%
S = zeros(size(M)); % Dry signal frequency-domain vectors.
R = zeros(size(M)); % Reverberated components frequency-domain vectors.

H_pow = zeros(frequencyCount, B); % Reverberant system frequency response estimate blocks - power.
H_pow = maxHEstimate;
pow = @(x) abs(x).^2;

previousSFramesPower = zeros(frequencyCount, B); % Most recent previous frame is first.
previousMFramesPower = ones(frequencyCount, B); % Most recent previous frame is first.

C = zeros(size(H_pow)); % Reverberant system frequency response power estimate blocks - temporary.
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
        % TODO: Temporal smoothing (Page 11).
        alpha = magnitudeSmoothingFactor(:, blockIndex);
        H_pow(:, blockIndex) = alpha .* H_pow(:, blockIndex) + (1-alpha) .* C(:, blockIndex);
    end
    H_pow_all(inputFrameIndex, :, :) = H_pow;
    
    newG_S = 1 - (sum(previousSFramesPower .* H_pow)) ./ inputFramePower;
    %G_S = sqrt(G_S);
    % Enforce a minimum gain for each frequency.
    newG_S = max([newG_S minGain], [], 2);
    G_S = (1-gainSmoothingFactor).*G_S + gainSmoothingFactor.*newG_S;
    newG_R = 1 - G_S;
    G_R = (1-gainSmoothingFactor).*G_R + gainSmoothingFactor.*newG_R;
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
figure;
spectrogram(dry, window, overlapSamples, 'yaxis');
title("Dry - reconstructed")
%%
reverberant = reconstruct(R, window, overlapSamples, inputLength);
soundsc(reverberant, fs);
spectrogram(reverberant, window, overlapSamples, 'yaxis');
title("Reverberant components - reconstructed");
%% Low-pass filter reverberant components
% Doesn't help.
reverberantFiltered = filter(lowpassfilter, reverberant);
spectrogram(reverberantFiltered, window, overlapSamples, 'yaxis');
soundsc(reverberantFiltered, fs);
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
%% Draw movie frames
clear F;
frameHop = 5;
movieFrameCount = floor(frameCount / frameHop);
F(movieFrameCount) = struct('cdata',[],'colormap',[]);
image = spectrogramPlot(squeeze(H_pow_all(1, :, :)), t, w);
for i = 1:frameHop:frameCount
    image.CData = squeeze(pow2db(H_pow_all(i, :, :)));
    drawnow();
    F((i - 1)/ frameHop + 1) = getframe(gcf);
end
%%
% Display movie interactively.
figure;
movie(F);
%%
% Write movie to file.
video = VideoWriter('estimate_evolution.avi');
open(video);
writeVideo(video, F);
close(video)

%%
% Length of estimated impulse response: 
stftHop = stftWindowSize - overlapSamples;
impulseResponseLengthSamples = stftHop * (B+1) + overlapSamples;
impulseResponseLengthSeconds = impulseResponseLengthSamples / fs;
impulseResponseLengthSeconds

%%
audiowrite('results/dry_B400_g03_a02_m09.wav', dry ./ (max(abs(dry))), fs);
audiowrite('results/reverberant_B400_g03_a02_m09.wav', reverberant ./ (max(abs(reverberant))), fs);

%% Testing spectrogramPlot

spectrogramPlot(M, t, w);
figure;
spectrogram(signal, window, overlapSamples, 'yaxis', 'power');
%%
[M2,w2,t2] = spectrogram(signal, window, overlapSamples, [], []);

%% Debug impulse response read
irDebug = readImpulseResponse('audio/stalbans_a_mono.wav', B, window, overlapSamples);
spectrogramPlot(irDebug, w, t);
spectrogramPlot(H_pow, w, t);

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

%%
function blocks = readImpulseResponse(irFilename, blockCount, window, overlapSamples)

    impulseResponse = audioread(irFilename);
    %assert(irFs == fs);
    allBlocks = spectrogram(impulseResponse, window, overlapSamples);
    pow = @(x) abs(x).^2;
    blocks = pow(allBlocks(:, 1:blockCount));

end