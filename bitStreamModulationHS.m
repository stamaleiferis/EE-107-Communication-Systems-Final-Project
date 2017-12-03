function [y, t] = bitStreamModulationHS(b,T,fs)


	b(b==0) = -1;					%	convert 0's to -1's
	pulse = halfSineWave(T,fs);		%	get pulse shape
% 	y = zeros(1,T*fs*length(b));	% 	initialize output vector
% 	pulse_length = length(pulse);	%		
% 	y = pulse*b(1);
% 
% 	for i=2:length(b)
% 		x = b(i) * pulse;			% 	shape pulse according to bits
% 		y = [y(1:end-1) x];
%     end
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