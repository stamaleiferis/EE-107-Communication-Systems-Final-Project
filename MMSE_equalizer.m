function out = MMSE_equalizer(ch_coeff, fs,noise_power,signal)
    
    channel = getChannel(ch_coeff,fs);
      
    if nargin == 3 % plot frequency response 
       [H,f] = freqz(channel,1,1024,fs);
        Q = conj(H)./ (abs(H).^2 + 2*noise_power);
        figure
        subplot(2,1,1)
        plot(f,mag2db(abs(Q)))
        title('MMSE magnitude response')
        subplot(2,1,2)
        plot(f,angle(Q))
        title('MMSE phase')
        figure
        stem(ifft(Q))
        title('MMSE impulse response')
        
    elseif nargin == 4
        H = fft(channel, length(signal));
        Q = conj(H)./ (abs(H).^2 + 2*noise_power);
        out = ifft(Q .* fft(signal));
        out = out(1:end-fs/2+1);
    end

end