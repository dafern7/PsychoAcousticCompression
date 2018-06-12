%% Compression 

%target bitrate for this compression is around 170kbps

%{
[s,fs] = audioread('25Or6To4.wav'); 
info = audioinfo('entertainment.wav');
output = compressor(s,fs);
soundsc(output,fs)
%}
function [output] = compressor(s,fs)
[M,N] = size(s);
if N == 2
    s = s(:,end);
end
fn = fs/2;
ford = 511; %filter taps based on the standard coefficients(512)
Y = [];
UD = [];
for n = 1:32
    W = [(((n-1)*fn/32)+1e-7)/fn ((n*(fn/32))-1e-7)/fn]; %bandpass filters
    B = fir1(ford,W);
    y = filter(B,1,s); %filtering
    yd = downsample(y,32); %downsample to get the correct rate per band
    ud = upsample(yd,32); %testing
    Ud = filter(B,1,ud);  %testing
    Y = [Y yd];
    UD = [UD Ud]; %testing
end


%S = sum(Ud,2); %testing
%sound(S,fs)%testing

%maximal points on the power spectrum within each band show the level of
%each band. This will  be compared with the threshold of hearing for this
%psychoacoustic model
MP = [];
for j = 1:32
L = length(Y(:,j));
nfft = 512;
fsb = fft(Y(:,j),nfft);
Px = fsb.*conj(fsb)/(nfft*L);
mp = max(10*log10(Px/1e-12));
MP = [MP mp]; 
fVals = (0:nfft/2-1)/nfft;
end
T = [];
for n = 1:32
f = (n-1)*(fn/32)+1:n*(fn/32);
%threshold of hearing 
t = 3.64*((f/1000).^-0.8)-6.5*exp(-0.6*((f/1000)-3.3).^2)+(1e-3*((f/1000).^4));
t1 = mean(t); %taking the mean of the threshold of hearing within the band
T = [T t1]; %threshold of hearing per band (32 points for comparison)
end
%SMR calculations
for i = 1:32
    if T(i) < 0
        SMR(i) = MP(i);
    else
        SMR(i) = MP(i)/T(i);
    end
end
SMR;
%SMR is equal to SNR, and lower bits mean more quantization, lower SNR
for i = 1:32
    if SMR(i) < 1
        bpsamp(i) = 0;
    elseif SMR(i) < 8
        bpsamp(i) = 2;
    elseif SMR(i) < 16
        bpsamp(i) = 3;
    elseif SMR(i) < 24
        bpsamp(i) = 4;
    elseif SMR(i) < 30
        bpsamp(i) = 5;
    elseif SMR(i) < 36
        bpsamp(i) = 6;
    elseif SMR(i) < 42
        bpsamp(i) = 7;
    elseif SMR(i) < 48
        bpsamp(i) = 8;
    elseif SMR(i) < 54
        bpsamp(i) = 9;
    elseif SMR(i) < 60
        bpsamp(i) = 10;
    elseif SMR(i) < 66
        bpsamp(i) = 11;
    elseif SMR(i) < 72
        bpsamp(i) = 12;
    elseif SMR(i) < 78
        bpsamp(i) = 13;
    elseif SMR(i) < 84
        bpsamp(i) = 14;
    elseif SMR(i) < 90
        bpsamp(i) = 15;
    else
        bpsamp(i) = 16;
    end
end
bpsamp; %allocated the bits in each sample in each band
%sum(bpsamp)*(info.SampleRate/32) %125410 bits/sec 
%info.BitsPerSample*info.SampleRate %705600 bits/sec original bitRate

%quantization of each band for each of the determined bits per band
F = [];
for i = 1:32
    if bpsamp(i) == 16 %bit rate is originally 16, no compression needed
        f = Y(:,i);
    elseif bpsamp(i) == 0
        f = zeros(size(Y(:,i))); %most high frequency was below threshold of hearing
    else
        f = quantizer(Y(:,i),bpsamp(i));  
    end
    F = [F f];
end

%recombining each band, upsampling, then filtering for aliasing.
%basically the opposite of what was done in the beginning when bandpass
%filtering the original signal
%summing each component from each band yields reconstructed signal. (tested
%earlier)
FD = [];
for i = 1:32
    W = [(((i-1)*fn/32)+1e-7)/fn ((i*(fn/32))-1e-7)/fn];
    B = fir1(ford,W);
    fud = upsample(F(:,i),32);
    Fd = filter(B,1,fud);
    FD = [FD Fd];
end

output = sum(FD,2); 
%reconstructed signal
%soundsc(s,fs) %original signal
%Can hear the difference between the original signal and the compressed
%signal pretty significantly mostly because this program estimates the
%critical bands with the 32 bandpass filters

%threshold of hearing graphed alongside the sound input signal to show the
%masking 
nfft = 512;
fsb = fft(s,nfft);
Px = fsb.*conj(fsb)/(nfft*L);
fVals = fs*(0:nfft/2-1)/nfft;
f = fVals;
figure
plot(fVals,10*log10(Px(1:nfft/2)/1e-14));
axis([0 20000 -10 100])
hold on
T = 3.64*((f/1000).^-0.8)-6.5*exp(-0.6*((f/1000)-3.3).^2)+(1e-3*((f/1000).^4));
plot(fVals,T);
end


