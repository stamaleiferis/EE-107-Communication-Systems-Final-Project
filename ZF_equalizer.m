function [signal_out] = ZF_equalizer(ch_coeff, fs,signal)

    channel = getChannel(ch_coeff,fs);
    
    if nargin == 2   % plot equalizer frequency response
        [H,f] = freqz(channel,1,1024,fs);
        ZF = 1./H;   % zero forcing equalizer frequency response
        figure
        subplot(2,1,1)
        plot(f, mag2db(abs(ZF).^2))
        title('Zero-forcing equalizer frequency response')
        subplot(2,1,2)
        plot(f,angle(ZF))
        title('Zero-forcing equalizer phase')
        figure
        plot(ifft(ZF))
        title('Zero-forcing equalizer impulse response')
        
    elseif nargin == 3   
        H = fft(channel,length(signal));    %get channel frequency response
        %[H,f] = freqz(channel,1,length(signal),fs);
        ZF = 1./H;   % zero forcing equalizer frequency response
        signal_out = ifft(ZF.*fft(signal));    % return time domain signal at the output of the equalizer
        %signal_out = signal_out(1:length(signal)); % need to truncate the zeros in the end so that the eye diagram is good
   
    end
  
end