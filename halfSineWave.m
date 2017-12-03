function [g1, t1] = halfSineWave(T, fs)
    
    t = linspace(0,T-T/fs,fs);
    if nargout == 1 
        g1 = sin(pi*t/T);
    elseif nargout ==2
        g1 = sin(pi*t/T);
        t1 = t;
    end
end