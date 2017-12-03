clear all, close all
alpha = 2;    % roll-off factor
b = round(rand(20,1));  %random bitStream
% b = [0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
K = 2;  %truncation constant
fs = 32;    %samples per bit duration
T = 1;      %bit duratio 
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

%% Q2: Plot the modulated signals
figure
[HS_mod, t1] = bitStreamModulationHS(b,T,fs); % Modulate signal with Half Sine pulses
subplot(2,1,1)
plot(t1, HS_mod)    % Plot Half Sine modulated signal
title('Half Sine Modulated Signal')
[SRRC_mod,t2] = bitStreamModulationSRRC(b,T,fs,K,alpha); % SRRC modulated signal
subplot(2,1,2)
plot(t2, SRRC_mod)
title('SRRC Modulated Signal')

%% Q3 Plot spectrum of modulated signals 
figure
freqz(HS_mod)
title('Frequency Response of Half Sine Modulated Signal')
figure
freqz(SRRC_mod)
title('Frequency Response of SRRC Modulated Signal')

%% Q4: Plot eye diagram for modulated signals

%eye diagram for half sine modulated signal
eyediagram(HS_mod, fs-1, T, fs/2)
title('Eye diagram for half sine modulated signal')

%eye diagram for SRRC modulated signal
eyediagram(SRRC_mod(2*(K-1)*fs:end-2*(K-1)*fs),(K-1)*fs, T, 0)
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
freqz(channel) 
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

eyediagram(out_PS_HS,fs-1,T,fs/2+1)
title('Eye diagram for HS channel output')

eyediagram(out_PS_SRRC,(K-1)*fs,T,0)
title('Eye diagram for SRRC channel output')
%% Q7: Noise

noisePower = 0.01;   % noise power is sigma^2
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
eyediagram(out_PS_HS_withNoise,fs-1,T,fs/2+1)
title('Eye diagram for HS channel output with noise')
eyediagram(out_PS_SRRC_withNoise(2*(K-1)*fs:end-2*(K-1)*fs),(K-1)*fs,T,0)
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
HS_MF_out = conv(out_PS_HS,HS_MF); 
eyediagram(HS_MF_out,fs,T,0)

SRRC_MF_out = conv(out_PS_SRRC,SRRC_MF); 
eyediagram(SRRC_MF_out(2*(K-1)*fs:end-2*(K-1)*fs),(K-1)*fs,T,0)
