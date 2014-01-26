function [bpm, bpmsr] = dynamictempo(d,sr)
% [bpm, bpmsr] = dynamictempo(d,sr)
%   Make a time-varying estimate of the tempo of track <d,sr>.
%   Return BPM values in bpm, each corresponding to times in tt
% 2014-01-26 Dan Ellis dpwe@ee.columbia.edu HAMR2014-01

% Calculate the spectrogram
[onsetenv, oesr, D] = onsettfsurf(d,sr);

% Collapse to onset function
%mm = mean(max(0,diff(D,1,2)));

% dc-removed mm
%onsetenv = filter([1 -1], [1 -.99],mm);

% Spectrogram of onset fn
W = 8192;
H = 128;
O = abs(specgram(onsetenv,W,oesr,W,W-H));
[nr,nc] = size(O);
bpmsr = oesr/H;

imgsc([1:nc]/bpmsr, [0:nr-1]/(nr-1)*oesr/2, 20*log10(O));
caxis(max(caxis)+[-50 0])

% Find viterbi path
% Limit bins to 8 Hz
eighthzbin = round( (8.0/oesr * W)/2)*2 + 1;
O = O(1:eighthzbin, :);
nOb = size(O,1);

% Create transition matrix
txwidth = 1.0;
txr = max(0.002, exp(-abs([-(nOb-1)/2:(nOb-1)/2])/txwidth));
tpz = toeplitz(rot(txr,(nOb-1)/2));

% Run Viterbi path
vp = viterbi_path(ones(nOb,1), tpz, max(0,20*log10(O)));

% Convert to BPM
bpm = 60*(vp-1)/(nr-1)*oesr/2;

% Overplot on spectrogram
hold on; plot([1:length(vp)]/bpmsr, t/60, '-r'); hold off
ax = axis; axis([ax(1) ax(2) 0 8]);
