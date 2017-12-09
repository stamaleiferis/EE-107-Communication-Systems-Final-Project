function bits = sampleAndDetect(signal,T,fs)
    
    samples = [];
    bits = [];
    for i = 1:fs*T:length(signal)-fs*T
       samples(i) = signal((i+fs*T)/2);
    end
    bits = samples;
    bits(bits>0) = 1;
    bits(bits<0) = -1;


end