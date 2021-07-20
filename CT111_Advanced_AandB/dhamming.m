function y=dhamming(A,B)
    y=0;
    for i=1:length(A)
        y=y+mod(A(i)+B(i),2);
    end
end


