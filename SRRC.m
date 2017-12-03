function [g2, t2] = SRRC(alpha, T, K, fs)

	t = linspace(-K*T,K*T,2*K*fs);

	for i=1:length(t)

		if t(i)==0
			x(i) = 1 - alpha + 4*alpha/pi;

		elseif t(i) == abs(T/(4*alpha))
			x(i) = (alpha/sqrt(2)) * ((1+2/pi)*sin(pi/(4*alpha)) + (1-2/pi)*cos(pi/(4*alpha)));

		else
			x(i) = (sin(pi*(t(i)/T)*(1-alpha))+4*alpha*(t(i)/T)*cos(pi*(t(i)/T)*(1+alpha)))/((pi*t(i)/T)*(1-(4*alpha*t(i)/T)^2));
		
		end

	end

	HS = halfSineWave(T,fs);
	A = sqrt(sum(HS.^2)/sum(x.^2)); % Amplitude correction factor
    
    if nargout == 2
        g2 = A*x;
        t2 = t;
    elseif nargout == 1
        g2 = A*x;
    end
end