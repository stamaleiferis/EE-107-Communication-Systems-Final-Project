clear all, close all
alpha = 0.5;    % roll-off factor
b = round(rand(100,1));  %random bitStream
% b = [0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
K = 4;  %truncation constant
fs = 32;    %samples per bit duration
T = 1;      %bit duration 
%% Q1: Plot pulse shaping functions and their frequency responses
[pulse_HS, t1] = halfSineWave(T,fs);  %get half sine pulse and time vector
[pulse_SRRC, t2] = SRRC(alpha, T, K, fs);   %get SRRC pulse and time vector

figure
plot(t1, pulse_HS)      % plot Half Sine Pulse
title('Half Sine Wave vs Time')
xlabel('time(s)')
ylabel('Half Sine Wave')          
figure
plot(t2,pulse_SRRC);    % plot SRRC pulse
ttl = sprintf('SRRC vs Time K=%d', K);
title(ttl)
xlabel('time(s)')
ylabel('SRRC')

figure
freqz(pulse_HS);    % plot Half Sine Pulse Frequency Response
title('Half Sine Wave Frequency Response')
figure
freqz(pulse_SRRC);  % plot SRRC Pulse Frequency Response
title('SRRC Frequency Response')

%% Q3: Plot the modulated signals
figure
[HS_mod, t1] = bitStreamModulationHS(b,T,fs); % Modulate signal with Half Sine pulses
subplot(2,1,1)
plot(t1, HS_mod)    % Plot Half Sine modulated signal
title('Half Sine Modulated Signal')
[SRRC_mod,t2] = bitStreamModulationSRRC(b,T,fs,K,alpha); % SRRC modulated signal
subplot(2,1,2)
plot(t2, SRRC_mod)
title('SRRC Modulated Signal')

%% Q2 Plot spectrum of modulated signals 
figure
freqz(HS_mod)
title('Frequency Response of Half Sine Modulated Signal')
figure
freqz(SRRC_mod)
title('Frequency Response of SRRC Modulated Signal')

%% Q4: Plot eye diagram for modulated signals

%eye diagram for half sine modulated signal
eyediagram(HS_mod, fs, T, fs/2)
title('Eye diagram for half sine modulated signal')

%eye diagram for SRRC modulated signal
eyediagram(SRRC_mod(2*(K-1)*T*fs+1:end-2*(K-1)*T*fs), fs, T, 0)
%eyediagram(SRRC_mod, fs, T, 0)
title('Eye diagram for SRRC modulated signal')
%% Q5 Channel impulse and frequency response
ch_coeff = [1 1/2 3/4 -2/7];    % channel coefficients, first sample is the ideal channel, rest are echoes
channel = getChannel(ch_coeff,fs);    % giving it 2 arguments returns the channel zero padded according to the samples per bit duration
t_ch = linspace(0,length(ch_coeff)*T,length(channel));  % time vector to plot the channel response
figure
stem(t_ch,channel) 
title('Channel Impulse Response')
xlabel('Time')
ylabel('Channel')
figure
freqz(channel,1,1024);
[H, w] = freqz(channel,1,1024); 
title('Channel Frequency Response')

%% Channel output for 2 different modualted signals
out_PS_HS = getChannel(ch_coeff,fs, HS_mod);    %giving it 3 arguments, returns the signal convolved with channel response
t_out1 = linspace(0,(length(b)+length(ch_coeff))*T,length(out_PS_HS));  %time vector for plotting the modulated signal convolved with channel
figure
subplot(2,1,1)
plot(t_out1, out_PS_HS)
xlabel('Time')
ylabel('Channel output')
title('Channel output for modulated signal with Half Sine Wave')

out_PS_SRRC = getChannel(ch_coeff,fs, SRRC_mod);  %SRRC modulated signal after passing through channel
t_out2 = linspace(-K*T,(length(b)+length(ch_coeff)*T+K*T),length(out_PS_SRRC));
subplot(2,1,2)
plot(t_out2, out_PS_SRRC);
xlabel('Time')
ylabel('Channel output')
title('Channel output for modulated signal with SRRC')

eyediagram(out_PS_HS,fs,T,fs/2)
title('Eye diagram for HS channel output')
%eyediagram(out_PS_SRRC((2*K-1)*T*fs:end-(2*K-1)*T*fs), fs, T, 0)
eyediagram(out_PS_SRRC, fs, T, 0)
title('Eye diagram for SRRC channel output')
%% Q7: Noise

noisePower = 0.2;   % noise power is sigma^2
% Half Sine: Channel output added with noise
out_PS_HS_withNoise = addNoise(out_PS_HS, noisePower);    %add random noise to signal
figure
subplot(2,1,1)
plot(t_out1, out_PS_HS_withNoise)
title('HS Channel Output with Noise')
% SRRC: Channel output added with noise
out_PS_SRRC_withNoise = addNoise(out_PS_SRRC, noisePower); %add random noise to signal
subplot(2,1,2)
plot(t_out2, out_PS_SRRC_withNoise)
title('SRRC Channel Output with Noise')

%plot eye diagrams
eyediagram(out_PS_HS_withNoise,fs,T,fs/2)
title('Eye diagram for HS channel output with noise')
%eyediagram(out_PS_SRRC_withNoise((K-0.5)*T*fs:end-(K-0.5)*T*fs), fs, T, fs/2+2)
eyediagram(out_PS_SRRC_withNoise((2*K-1)*T*fs:end-(2*K-1)*T*fs), fs, T, 0)
%eyediagram(out_PS_SRRC_withNoise,fs, T, 0)
title('Eye diagram for SRRC channel output with noise')


%% Matched Filter Q8 & Q9

%Half sine matched filter
%It is the same pulse as the Half sine pulse
[HS_MF, t1] = HalfSineMatchedFilter(T, fs);
%SRRC matched filter
%same pulse shape as SRRC, shifted in time by KT to be causal
[SRRC_MF, t2] = SRRCMatchedFilter(alpha, T, K, fs);
%plot the impulse responses
figure
subplot(2,1,1)
plot(t1,HS_MF)  %HS
title('Half Sine Matched Filter')
subplot(2,1,2)
plot(t2, SRRC_MF)   %SRRC
title('SRRC Matched Filter')
%plot the frequency responses
figure
freqz(HS_MF)
title('Half Sine Matched Filter Frequency Response')
figure
freqz(SRRC_MF)
title('SRRC Matched Filter Frequency Response')

%convolve modulated signals with their respective matched filter
HS_MF_out = conv(out_PS_HS_withNoise,HS_MF);  % use signal at channel output instead of just modulated signal
%
%HS_MF_out = ifft(fft(out_PS_HS).*fft(HS_MF,length(out_PS_HS)));
%HS_MF_out = HS_MF_out(fs*T:end-fs*T);
eyediagram(HS_MF_out,fs,T,0)
title('Eye diagram for HS matched filter output')

SRRC_MF_out = conv(out_PS_SRRC,SRRC_MF); 
%eyediagram(SRRC_MF_out((2*K-1)*fs*T:end-T*(2*K-1)*fs),fs,T,0)
eyediagram(SRRC_MF_out,fs, T, 0)
title('Eye diagram for SRRC matched filter output')

%% Q10 & Q11 Zero-Forcing Equalizer
% out_PS_HS is channel output for half sine modulated signal
% out_PS_HS_withNoise
% HS_MF_out  is after channel and matched filter

% out_PS_SRRC is channel output for
% out_PS_SRRC_withNoise
% SRRC_MF_out is after channel and matched filter

% plot zero-forcing equalizer responses
ZF_equalizer(ch_coeff, fs); % plot frequency and impulse response

% Modulated signals after channel and matched filter, pass through equalizer
ZF_HS = ZF_equalizer(ch_coeff, fs, HS_MF_out);  % half sine
ZF_SRRC = ZF_equalizer(ch_coeff, fs, SRRC_MF_out);  % SRRC
figure
subplot(2,1,1)  %plot half sine
plot(ZF_HS)
title('Half sine modulated signal after channel and zero-forcing equalizer')
subplot(2,1,2)
plot(ZF_SRRC)   %plot SRRC
title('SRRC modulated signal after channel and zero-forcing equalizer')
% Modulated signals after channel with noise, pass through equalizer
ZF_HS_noise = ZF_equalizer(ch_coeff, fs, out_PS_HS_withNoise);  % half sine
ZF_SRRC_noise = ZF_equalizer(ch_coeff, fs, out_PS_SRRC_withNoise);  % SRRC
figure
subplot(2,1,1)  %plot half sine
plot(ZF_HS_noise)
title('Half sine modulated signal after channel w/ noise and zero-forcing equalizer')
subplot(2,1,2)
plot(ZF_SRRC_noise)   %plot SRRC
title('SRRC modulated signal after channel w/ noise and zero-forcing equalizer')

% plot eye diagrams at the output of equalizer
eyediagram(ZF_HS(length(ZF_HS)/4:length(ZF_HS)/2),fs,T,1)     % for half sine
title('Eye diagram of the zero-forcing equalizer output for Half Sine modulated signal ')
eyediagram(ZF_SRRC,fs,T,0)      % for SRRC
title('Eye diagram of the zero-forcing equalizer output for SRRC modulated signal ')
% plot eye diagrams for signals with noise at the output of the equalizer
eyediagram(ZF_HS_noise(length(ZF_HS_noise)/4:length(ZF_HS_noise)/2),fs,T,25)     % for half sine
title('Eye diagram of the zero-forcing equalizer output for Half Sine modulated signal w/ noise')
eyediagram(ZF_SRRC_noise,fs,T,0)      % for SRRC
title('Eye diagram of the zero-forcing equalizer output for SRRC modulated signal w/ noise ')
%% Q12 & Q13 MMSE Equalizer
MMSE_equalizer(ch_coeff, fs,noisePower);    % plot MMSE frequency and impulse responses
MMSE_HS = MMSE_equalizer(ch_coeff, fs,noisePower,HS_MF_out);  % pass HS signal with noise through the equalizer
MMSE_SRRC = MMSE_equalizer(ch_coeff, fs,noisePower,SRRC_MF_out);  % pass the SRRC signal with noise through the equalizer
MMSE_SRRC = MMSE_SRRC(K*T*fs:end);  %remove leading zeros
%plot eye diagrams at the output of the equalizer
eyediagram(MMSE_HS(length(MMSE_HS)/4:length(MMSE_HS)/2),fs,T,5)     % for half sine
title('Eye diagram of the MMSE equalizer output for Half Sine ')
eyediagram(MMSE_SRRC(length(MMSE_SRRC)/4:length(MMSE_SRRC)/2),fs,T,21)        % for SRRC
title('Eye diagram of the MMSE equalizer output for SRRC ')
