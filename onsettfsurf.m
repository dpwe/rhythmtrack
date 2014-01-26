function [onsetenv, oesr, D] = onsettfsurf(d,sr, DISPLAY)
% [onsetenv, oesr, D] = onsettfsurf(d,sr, DISPLAY)
%   Onset envelope and spectrogram surface of D conditioned for onset detection
% 2014-01-26 Dan Ellis

if nargin < 3; DISPLAY = 0; end

sro = 22050;
swin = 256; % 256/22050 = 11.6 ms
shop = 64;  %  32/22050 =  1.45 ms
% so the 50% overlapped window is...
sstep = swin/shop
% mel channels
nmel = 10;
% sample rate for specgram frames (granularity for rest of processing)
oesr = sro/shop;

% resample to sro
if (sr ~= sro)
  gg = gcd(sro,sr);
  d = resample(d,sro/gg,sr/gg);
  sr = sro;
end

% prepad with silence to make onset envelope line up with time
% domain correctly.  

% We prepad so the very first frame is actually a full hop
% *prior* to the first sample (just the extra bit inside), 
% (i.e. centered on -thop/2)
% so when we take the first difference, it will correspond to the
% increment in energy precisely at time 0

% Make an swin point window with shop-sized half-cosines at each
% end, surrounding a flat core
mywin = hanning(2*shop)';
mywin = [mywin(1:shop), ones(1,swin-2*shop), mywin(shop+1:end)];
D = stft([zeros(swin-shop/2,1);d;zeros((swin-shop)/2,1)],swin,mywin,shop);
% This is the correct way to view it:
% imgsc(([1:size(D,2)]-(swin-shop/2+swin/2)/shop)/oesr, [0:size(D,1)-1]*sr/(2*(size(D,1)-1)), D)

% Construct db-magnitude-mel-spectrogram
mlmx = fft2melmx(swin,sr,nmel);
D = 20*log10(max(1e-10,mlmx(:,1:(swin/2+1))*abs(D)));

% Only look at the top 60 dB
D = max(D, max(D(:))-60);

% The raw onset decision waveform
%mm = mean(max(0, diff(D,1,2))); 
mm = mean(max(0, D(:, sstep:end) - D(:, 1:end-(sstep-1))));
eelen = length(mm);

% dc-removed mm
%onsetenv = filter([1 -1], [1 -.99],mm);
% remove mean and normalize amplitude in local 3 sec window
onsetenv = localmvnorm(mm, round(3.0*oesr));

if DISPLAY
  subplot(211)
  imgsc(([1:size(D,2)]-(swin-shop/2+swin/2)/shop)/oesr, ...
        [0:size(D,1)-1]*sr/(2*(size(D,1)-1)), D);
  subplot(212)
  plot([0:length(d)-1]/sr, d, ...
       [0:length(onsetenv)-1]/oesr, onsetenv);
  linkaxes([subplot(211),subplot(212)], 'x')
end
