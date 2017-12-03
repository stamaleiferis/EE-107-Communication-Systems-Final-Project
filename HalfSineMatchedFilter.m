function [HS_MF, t] = HalfSineMatchedFilter(T, fs)
    
    if nargout == 2
        HS_MF = halfSineWave(T,fs); %It is the same pulse as the Half sine pulse
        t = linspace(0,T-T/fs,fs);
    else
        HS_MF = halfSineWave(T,fs); %It is the same pulse as the Half sine pulse
    end
end