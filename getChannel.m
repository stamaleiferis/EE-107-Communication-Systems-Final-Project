function [out] = channel(ch_coeff,fs,signal)

zero = zeros(1, fs-1);
numTaps = length(ch_coeff);
channel = [];
output = 0;
for i = 1:numTaps
    channel = [channel ch_coeff(i) zero];  
end

if nargin == 2 && nargout == 1  % return channel with zero padding
    out = channel;
    return;
    
elseif nargin == 3 && nargout == 1
    out = conv(signal, channel); % return signal with channel effect; 
    return; 
    
end