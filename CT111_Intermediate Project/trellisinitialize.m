function [path_m,next,branch_m,s]=trellisinitialize(k,K) 
    if(K==3)
        % 4 states of the path_m as the constraint lenght is 3 in this case(2^(3-1))
        s(1)=1;          % corresponds to state [0 0]
        s(2)=2;          % corresponds to state [0 1]
        s(3)=3;          % corresponds to state [1 0]
        s(4)=4;          % corresponds to state [1 1]
    
    elseif(K==4)
        % 8 states of the path_m as the constraint lenght is 4 in this case (2^(4-1))
        s(1)=1;          % corresponds to state [0 0 0]
        s(2)=2;          % corresponds to state [0 0 1]
        s(3)=3;          % corresponds to state [0 1 0]
        s(4)=4;          % corresponds to state [0 1 1]
        s(5)=5;          % corresponds to state [1 0 0]
        s(6)=6;          % corresponds to state [1 0 1]
        s(7)=7;          % corresponds to state [1 1 0]
        s(8)=8;          % corresponds to state [1 1 1]   
    end
    path_m=-ones(2^(K-1),k+K);         %initializing the memory for reducing complexity
    branch_m=zeros(2*(2^(K-1)),k+K-1)+9;           % just initializing it with some high constant 
    % Some places of the path_m initialized with appropriate values
    path_m(2^(K-1),1)=0;             
    path_m(1:2^(K-1)-1,1)=inf;
    next=zeros(2^(K-1),2)+Inf;
    % Giving some intial conditions before starting the loop
    next(1,1)=0;
