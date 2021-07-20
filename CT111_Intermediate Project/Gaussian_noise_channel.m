function gec=Gaussian_noise_channel(c,k,K,rate,SNRinwatt)
    gec=2*c-1;                % to convert 0 to -1 and keeping 1 as 1 only as per BPSK
    sigma2=1/(2*rate*SNRinwatt);
    sigma=sqrt(sigma2);
    noise=sigma*randn(1,(1/rate)*(k+K-1));
    gec=gec+noise;
end