function z=dbec(A,B)
    y=0;
    for i=1:length(A)
        if(A(1,i)~=B(1,i))
            y=y+1;
        end
    end
    z=y;
end