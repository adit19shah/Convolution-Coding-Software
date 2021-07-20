function bec=bec_channel(c,k,K,rate,SNRinwatt)
    pError=qfunc(sqrt(2*rate*SNRinwatt));            % probability of the bit getting erased
    bec=c;
    for i=1:(1/rate)*(k+K-1)
      if rand<pError
        bec(i)=25;           % replaced with some constant to show that this bit was erased
      end
    end
end