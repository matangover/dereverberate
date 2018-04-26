# Audio examples
This folder contains sound files that are used by the algorithm or can be used to test the algorithm.

## Impulse response
`stalbans_a_mono.wav` is an impulse response of The Lady Chapel, St Albans Cathedral retrieved from the [OpenAIR project](http://www.openairlib.net/auralizationdb/content/lady-chapel-st-albans-cathedral). The chapel is a very reverberant space. This impulse response is used by the dereverberation algorithm in `reverberate.m` as the maximum limit on the impulse response estimation. It is also used by `convolutionReverb.m` to add artificial reverb to sound recordings.
## Reverberated/dry audio samples
 - `EchoSample.mp3`: an 8-second speech sample containing a large amount of reverb. Taken from a personal website of an Audacity forum member ([original location](http://kozco.com/tech/audacity/clips/EchoSample.mp3)) who [says](http://forum.audacityteam.org/viewtopic.php?p=316232#p316232) "If you have a recording made like this... then you're stuck".
 - `singing.mp3`: a 5-second female operatic voice anechoic recording, from the [OpenAIR project](http://www.openairlib.net/anechoicdb/content/operatic-voice).
 - `stalbans_omni_sing.mp3`: `singing.mp3` with added artificial reverb using convolution with the impulse response from the St Albans Cathedral.
 - `mozart.mp3`: a soprano aria with orchestral accompaniment from Mozart's _Don Giovanni_. Taken from a collection of [anechoic multi-track recordings of symphonic music](http://research.cs.aalto.fi/acoustics/virtual-acoustics/research/acoustic-measurement-and-analysis/85-anechoic-recordings.html) made at Aalto University, and mixed into a single track.
 - `mozart_reverb.mp3`: `mozart.mp3` with added artificial reverb using convolution with the impulse response from the St Albans Cathedral.
 - `mozart_reverb_short.mp3`: the first 18 seconds of `mozart_reverb.mp3`.
