function bsc=bsc_channel(c,k,K,rate,SNRinwatt)
    pFlip=qfunc(sqrt(2*rate*SNRinwatt));            % probability of the bit getting flipped
    bsc=c;
    for i=1:(1/rate)*(k+K-1)
      if rand<pFlip
        bsc(i)=1-bsc(i);
      end
    end
end