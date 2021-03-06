function bits = sampleAndDetect(signal,T,fs,HS,numBits,K)
    % signal - Half Sine or SRRC to be samples
    % T      - Bit duration
    % fs     - number of samples per bit duration
    % HS     - true is signal to be sampled was modulated with HS
    % numBits- number of bits the signal is carrying
    
    % zero thresholding will be applied
    
    bits = zeros(numBits,1); %recovered bits will be stored here
    if HS
        for i=1:numBits
            sample = signal(i*T*fs);
            if sample > 0
                bits(i) = 1;
            elseif sample < 0
                bits(i) = 0;
            else
                bits(i) = round(rand);
            end
        end
        
        
    else %its srrc
        offset = 0;
        signal = signal(2*(K)*fs-fs:end);
        for i=1:numBits
            sample = signal(i*T*fs);
            %p = plot(signal(1+i*T*fs),'o','MarkerFaceColor','red');
            %drawnow
            if sample > 0
                bits(i) = 1;
            elseif sample < 0
                bits(i) = 0;
            else
                bits(i) = round(rand);
            end
        end
        
    end

end