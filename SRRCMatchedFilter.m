function [SRRC_MF, t2] = SRRCMatchedFilter(alpha, T, K, fs)
    if nargout == 2 
        SRRC_MF = SRRC(alpha, T, K, fs);
        t2 = linspace(-K*T+T,K*T+T,2*fs*K); %shift in time by KT
    else
        SRRC_MF = SRRC(alpha, T, K, fs);
    end
end