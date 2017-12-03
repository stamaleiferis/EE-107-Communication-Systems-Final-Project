function signalWithNoise = addNoise(signal, noisePower)

sigma = sqrt(noisePower);
noise = sigma * randn(length(signal),1)';

signalWithNoise = signal + noise;


end