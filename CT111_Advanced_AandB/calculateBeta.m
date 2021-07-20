function beta=calculateBeta(i,gamma,beta,k)
    if(i==k+3)
        beta(1:3,i)=0;
        beta(4,i)=1;
    else
        beta=calculateBeta(i+1,gamma,beta,k);
        beta(4,i)=gamma(8,i)*beta(4,i+1)+gamma(7,i)*beta(2,i+1);
        beta(3,i)=gamma(6,i)*beta(4,i+1)+gamma(5,i)*beta(2,i+1);
        beta(2,i)=gamma(4,i)*beta(3,i+1)+gamma(3,i)*beta(1,i+1);
        beta(1,i)=gamma(2,i)*beta(3,i+1)+gamma(1,i)*beta(1,i+1);
    end
end