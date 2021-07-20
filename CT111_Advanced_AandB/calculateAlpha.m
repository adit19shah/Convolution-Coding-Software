function alpha=calculateAlpha(i,gamma,alpha)
    if(i==1)
        alpha(1:3,1)=0;
        alpha(4,1)=1;
    else
        alpha=calculateAlpha(i-1,gamma,alpha);
        alpha(4,i)=gamma(8,i-1)*alpha(4,i-1)+gamma(6,i-1)*alpha(3,i-1);
        alpha(3,i)=gamma(4,i-1)*alpha(2,i-1)+gamma(2,i-1)*alpha(1,i-1);
        alpha(2,i)=gamma(7,i-1)*alpha(4,i-1)+gamma(5,i-1)*alpha(3,i-1);
        alpha(1,i)=gamma(3,i-1)*alpha(2,i-1)+gamma(1,i-1)*alpha(1,i-1);
    end
end