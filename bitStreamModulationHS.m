function [y, t] = bitStreamModulationHS(b,T,fs)


	b = 2*b-1;					%	convert 0's to -1's
	pulse = halfSineWave(T,fs);		%	get pulse shape

space = zeros(1,fs-1);
convDum = [];

for i = 1:length(b)-1
    convDum = [convDum b(i) space];
end
convDum = [convDum b(end)];

y = conv(pulse,convDum);
    
    if nargout == 2
        t = linspace(0,length(b)*T-T/fs,length(y));    % time vector for HS signal

end