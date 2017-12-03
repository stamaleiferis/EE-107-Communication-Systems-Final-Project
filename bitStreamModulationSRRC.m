function [y, t] = bitStreamModulationSRRC(b,T,fs,K,alpha)

	b(b==0) = -1;					%	convert 0's to -1's
	pulse = SRRC(alpha, T, K, fs);		%	get pulse shape
	pulse = pulse(1:end);
	spacing = zeros(1, (T*fs)-1);

	dataToConvolve = [];
	for i=1:length(b)-1
		dataToConvolve = [dataToConvolve b(i) spacing];
    end
    dataToConvolve = [dataToConvolve b(end)];
	y = conv(dataToConvolve,pulse);

    if nargout == 2
        t = linspace(0,length(b)*(2*K)*T,length(y));
    end

end