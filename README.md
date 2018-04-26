# Dereverberate
This repository implements a sound dereverberation algorithm proposed by Gilbert Soulodre:

> Soulodre, Gilbert A. 2010. “About this dereverberation business: A method for extracting reverberation from audio signals.” In Audio Engineering Society Convention 129. Audio Engineering Society.

Note that this algorithm is patented (US Patent [8,036,767 B2](https://patents.google.com/patent/US8036767B2/en)).

See `dereverberate.m` for the algorithm implementation in MATLAB code. See `audio` for input files used to evaluate the algorithm and the implementation.

Here is a spectrogram representation of an example input (`audio/EchoSample.mp3`) and the algorithm's results on that input:

![Algorithm results](results/speech/results.png?raw=true)
To hear the results, see the `results` folder.
