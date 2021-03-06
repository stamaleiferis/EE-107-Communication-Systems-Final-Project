% run the system 
clear all, close all
alpha = 0.5;    % roll-off factor
K = 5;  %truncation constant
fs = 32;    %samples per bit duration
T = 1;      %bit duration 
ch_coeff = [1 1/2 3/4 -2/7]; %channel taps
noisePower = 0.1;   % AWGN noise power is sigma^2
%% Image pre-processing
% filename = 'onion.png';
filename = 'autumn.tif';
qbits = 8;
[Ztres,r,c,m,n,minval,maxval]=ImagePreProcess_gray(filename,qbits); % returns array of 8x8 blocks
% [Ztresr,Ztresg,Ztresb,r,c,m,n,minval,maxval]
%% Convert to bit stream
sZ = size(Ztres);
possibleGroupNums = factor(sZ(3)); % in ascending order
N = possibleGroupNums(end);
[stream,chunkLength] = convertToBitStream(Ztres,N);
numChunks = length(stream)/chunkLength;
recoverStreamHS = zeros(size(stream));
recoverStreamSRRC = zeros(size(stream));
recoverStreamHS_ZF = recoverStreamHS;
recoverStreamSRRC_ZF = recoverStreamSRRC;
numBits = chunkLength;

for i=1:numChunks
    % Take the desired amount of the stream
    b = double(stream((i-1)*chunkLength+1:i*chunkLength));
    %% Modulate bitstream
    [HS_mod, t1] = bitStreamModulationHS(b,T,fs); % Modulate signal with Half Sine pulses
    [SRRC_mod,t2] = bitStreamModulationSRRC(b,T,fs,K,alpha); % SRRC modulated signal
    %% Pass modulated bitstream through channel
    out_PS_HS = getChannel(ch_coeff,fs, HS_mod);
    out_PS_SRRC = getChannel(ch_coeff,fs, SRRC_mod);
    %% Add noise
    out_PS_HS = addNoise(out_PS_HS, noisePower);
    out_PS_SRRC = addNoise(out_PS_SRRC, noisePower);
    %% Pass signal through matched filter
    HS_MF = HalfSineMatchedFilter(T, fs); %this is the matched filter
    HS_MF_out = conv(out_PS_HS,HS_MF);
    SRRC_MF = SRRCMatchedFilter(alpha, T, K, fs); %this is the matched filter
    SRRC_MF_out = conv(out_PS_SRRC,SRRC_MF);
    %% Pass signal through equalizer
    % Zero-forcing
    ZF_HS = ZF_equalizer(ch_coeff, fs, HS_MF_out);  % half sine
    ZF_SRRC = ZF_equalizer(ch_coeff, fs, SRRC_MF_out);  % SRRC
    % or MMSE
    MMSE_HS = MMSE_equalizer(ch_coeff, fs,noisePower,HS_MF_out);  % half sine
    MMSE_SRRC = MMSE_equalizer(ch_coeff, fs,noisePower,SRRC_MF_out);  % SRRC
    %MMSE_SRRC(2*(K-1)*fs:end); % remove leading zeros
    %% Sample 
    bits_ZF_HS = sampleAndDetect(ZF_HS,T,fs,1,numBits,K);
    bits_ZF_SRRC = sampleAndDetect(ZF_SRRC,T,fs,0,numBits,K);
    bits_MMSE_HS = sampleAndDetect(MMSE_HS,T,fs,1,numBits,K);
    bits_MMSE_SRRC = sampleAndDetect(MMSE_SRRC,T,fs,0,numBits,K);
    %% Recover bitstream
    recoverStreamHS_ZF((i-1)*chunkLength+1:i*chunkLength) = bits_ZF_HS;
    recoverStreamSRRC_ZF((i-1)*chunkLength+1:i*chunkLength) = bits_ZF_SRRC;
    recoverStreamHS((i-1)*chunkLength+1:i*chunkLength) = bits_MMSE_HS;
    recoverStreamSRRC((i-1)*chunkLength+1:i*chunkLength) = bits_MMSE_SRRC;
end
%% Convert back to image
ZtresHS_MMSE = convertFromBitStream(recoverStreamHS, sZ(1), sZ(2), sZ(3), qbits);
ZtresSRRC_MMSE = convertFromBitStream(recoverStreamSRRC, sZ(1), sZ(2), sZ(3), qbits);
ZtresHS_ZF = convertFromBitStream(recoverStreamHS_ZF, sZ(1), sZ(2), sZ(3), qbits);
ZtresSRRC_ZF = convertFromBitStream(recoverStreamSRRC_ZF, sZ(1), sZ(2), sZ(3), qbits);
%% Image post processing
ImagePostProcess_gray(ZtresHS_MMSE,r,c,m,n,minval,maxval)
title('HS Recovered Image, MMSE Equalizer')
ImagePostProcess_gray(ZtresSRRC_MMSE,r,c,m,n,minval,maxval)
title('SRRC Recovered Image, MMSE Equalizer')
ImagePostProcess_gray(ZtresHS_ZF,r,c,m,n,minval,maxval)
title('HS Recovered Image, ZF Equalizer')
ImagePostProcess_gray(ZtresSRRC_ZF,r,c,m,n,minval,maxval)
title('SRRC Recovered Image, ZF Equalizer')