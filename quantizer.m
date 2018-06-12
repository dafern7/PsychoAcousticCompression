%% Quantization


%{
s = audioread('25Or6To4.wav');
info = audioinfo('25Or6To4.wav')
fs = 44100;
bit = 8;
%}
%output = quantiser(s,bit);
%soundsc(output,fs)
%quantization method reference: Fabien Petitcolas in his MPEG for MATLAB
%program on his website www.petitcolas.net/fabien/software/mpeg/index.html
function output = quantiser(s,bit)
%this is to ensure that the size of the matrix is 1x(size) not (size)x1
[M,N] = size(s);
if M > 1
    s = s';
end

%denote range, then separate the total range into step size based on bits
%allowed. The higher the number of bits, the more "steps" are allowed, and
%a greater approximation of the original signal is allowed. Inversely, a
%bigger step size yields a lesser approximation of the original signal.
r = max(s)-min(s); %range
d = roundfun(r/2^bit); %this is the "step size" rounded to 1e-4.
b = roundfun(min(s)); %minimum of the signal rounded.
q = transpose(b:d:roundfun(max(s) - d)); %these are the quantization values that are allowed
%based on the number of bits
%this is to make input and quantized input the same size 
input = repmat(s,length(q),1);
quantinput = repmat(q,1,length(s));
%the set of vectors that represents difference between the quantized input and the actual input
qdiff = abs(input-quantinput);
%minimum difference gives quantization 
[Y,index] = min(qdiff); 
quantResults = quantinput(index);
for i = 1:length(quantResults) 
    levels(i) = find(quantResults(i) == q); %finding the indices of the levels to be used in scaling
end
maxlevel = max(levels);
%scale the levels in a similar way to quantization above, except scale up
%to the max level

scaledLevels = transpose(b:d:roundfun(b+d*(maxlevel))); 
output = scaledLevels(levels); 
end




