function c=Encoder(m,k,K,r)
    new_m=[zeros(1,K-1),m,zeros(1,K-1)];    % Appending and Prepending zeros
    if(K==3)
        G=[1 1 1;1 0 1];        % Given in project Manual
    elseif(K==4 && r==0.5)
        G=[1 1 1 1;1 1 0 1];    % Source: Proakis book
    elseif(K==4 && (r==1/3||r==0.33||r==0.34))
        G=[1 1 1 1;1 1 0 1;1 0 1 1];    %Source:Proakis book
    end
    %Memory allocation for encoded message
     c=zeros(1,(1/r)*(k+K-1));
    % Core Encoding Procedure in a way similar to Linear Shift Register
     for i=K:length(new_m)
         temp=new_m(1,i-(K-1):i);
         c(1,(1/r)*(i-K)+1:(1/r)*(i-K)+(1/r))=mod(temp*flip(G'),2);
     end
end