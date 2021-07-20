% CT111 Final Project   Developed by: Adit Shah   Topic: Convolution coding

% To clear all the variables,clear command window and close all previous
% graphs
clearvars; close all; clc;

% To take message size as input from the user
k=input("Enter the length of the message to be transmitted      : "); 
K=input("Enter the constraint length of the message             : ");
rate=input("Please enter the rate with which you want to encode(1/2 or 1/3): "); 

% Generating message m of length k
fprintf("The message string to be transmitted is: \n");
m=randi(2,1,k)-1;       %randi(2,1,k) generates random matrix of dimensions 1 x k and containing integer values between [1,2]
disp(m);
c=Encoder(m,k,K,rate);  % Encodes the randomly generated message m

% Fixing the SNR ratio
SNR=5;                    % SNR in db
SNRinwatt=10^(SNR/10);    % converting SNR from db to watt

% when the encoded message passes through the Gaussian Noise channel
gec=Gaussian_noise_channel(c,k,K,rate,SNRinwatt);

% When the encoded message passes through the BSC channel
bsc=bsc_channel(c,k,K,rate,SNRinwatt);

% when the encoded message passes through the BEC channel
bec=bec_channel(c,k,K,rate,SNRinwatt);

% Decoder for the BSC Channel [Using Viterbi Algorithm]

% Initializing some necessary values for Trellis
[path_m,next,branch_m,s]=trellisinitialize(k,K);

if(K==3)
    % Iterating over the bits of received message in pairs of (1/rate=2)
    for t=1:2:(2*k+4)
      temp=[bsc(t) bsc(t+1)];             % pair formation of 2 bits 
      % S_1
      curr=min(next(1,1),next(2,1));
      % Finding path metrics of that state 
      path_m(4,(t+1)/2)=curr;       
      branch_m(7,(t+1)/2)=dhamming([1,1],temp);
      branch_m(8,(t+1)/2)=dhamming([0,0],temp);
      %Finding branch metric of branch taken if bit 1 comes next 
      %and adding to path metric of current state which is stored in "curr"
      temp_next_12=dhamming([1 1],temp)+curr; 
      %Finding branch metric of branch taken if bit 0 comes next 
      % and then adding it to path metric of current state
      next(1,1)=dhamming([0 0],temp)+curr;        

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(3,(t+1)/2)=curr;
      branch_m(5,(t+1)/2)=dhamming([0,0],temp);
      branch_m(6,(t+1)/2)=dhamming([1,1],temp);
      temp_next_22=dhamming([0 0],temp)+curr;
      next(2,1)=dhamming([1 1],temp)+curr;

      %s_3
      curr=min(next(2,2),next(1,2));
      path_m(2,(t+1)/2)=curr;
      branch_m(3,(t+1)/2)=dhamming([0,1],temp);
      branch_m(4,(t+1)/2)=dhamming([1,0],temp);
      temp_next_32=dhamming([0 1],temp)+curr;
      next(3,1)=dhamming([1 0],temp)+curr;

      %s_4
      curr=min(next(3,2),next(4,2));
      path_m(1,(t+1)/2)=curr;
      branch_m(1,(t+1)/2)=dhamming([1,0],temp);
      branch_m(2,(t+1)/2)=dhamming([0,1],temp);
      next(4,2)=dhamming([1 0],temp)+curr;
      next(4,1)=dhamming([0 1],temp)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
    end
    path_m(4,k+3)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state
    
elseif(K==4 && rate==1/2)
    for t=1:2:(2*k+6)
      temp=[bsc(t) bsc(t+1)];                     % forming pairs of 2 bits from the sequence received from the channel 
      % S_1
      curr=min(next(1,1),next(2,1));        % next(1,1):from state-1 carrying 0
      path_m(8,(t+1)/2)=curr;                    % Finding path metrics of that state
      branch_m(15,(t+1)/2)=dhamming([1,1],temp);
      branch_m(16,(t+1)/2)=dhamming([0,0],temp);
      temp_next_12=branch_m(15,(t+1)/2)+curr;     %Finding branch metrics and adding to path metric of current state
      next(1,1)=branch_m(16,(t+1)/2)+curr;        %Finding branch metrics and adding to path metric of curent state

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(7,(t+1)/2)=curr;
      branch_m(13,(t+1)/2)=dhamming([0,0],temp);
      branch_m(14,(t+1)/2)=dhamming([1,1],temp);
      temp_next_22=branch_m(13,(t+1)/2)+curr;
      next(2,1)=branch_m(14,(t+1)/2)+curr;
      
      %s_3
      curr=min(next(5,1),next(6,1));
      path_m(6,(t+1)/2)=curr;
      branch_m(11,(t+1)/2)=dhamming([0,1],temp);
      branch_m(12,(t+1)/2)=dhamming([1,0],temp);
      temp_next_32=branch_m(11,(t+1)/2)+curr;
      next(3,1)=branch_m(12,(t+1)/2)+curr;
      
      %s_4
      curr=min(next(7,1),next(8,1));
      path_m(5,(t+1)/2)=curr;
      branch_m(9,(t+1)/2)=dhamming([1,0],temp);
      branch_m(10,(t+1)/2)=dhamming([0,1],temp);
      temp_next_42=branch_m(9,(t+1)/2)+curr;
      next(4,1)=branch_m(10,(t+1)/2)+curr;
      
       %s_5
      curr=min(next(2,2),next(1,2));
      path_m(4,(t+1)/2)=curr;
      branch_m(7,(t+1)/2)=dhamming([0,0],temp);
      branch_m(8,(t+1)/2)=dhamming([1,1],temp);
      temp_next_52=branch_m(7,(t+1)/2)+curr;
      next(5,1)=branch_m(8,(t+1)/2)+curr;
      
      %s_6
      curr=min(next(3,2),next(4,2));
      path_m(3,(t+1)/2)=curr;
      branch_m(5,(t+1)/2)=dhamming([1,1],temp);
      branch_m(6,(t+1)/2)=dhamming([0,0],temp);
      temp_next_62=branch_m(5,(t+1)/2)+curr;
      next(6,1)=branch_m(6,(t+1)/2)+curr;

      %s_7
      curr=min(next(5,2),next(6,2));
      path_m(2,(t+1)/2)=curr;
      branch_m(3,(t+1)/2)=dhamming([1,0],temp);
      branch_m(4,(t+1)/2)=dhamming([0,1],temp);
      temp_next_72=branch_m(3,(t+1)/2)+curr;
      next(7,1)=branch_m(4,(t+1)/2)+curr;

      %s_8
      curr=min(next(7,2),next(8,2));
      path_m(1,(t+1)/2)=curr;
      branch_m(1,(t+1)/2)=dhamming([0,1],temp);
      branch_m(2,(t+1)/2)=dhamming([1,0],temp);
      next(8,2)=branch_m(1,(t+1)/2)+curr;
      next(8,1)=branch_m(2,(t+1)/2)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
      next(4,2)=temp_next_42;
      next(5,2)=temp_next_52;
      next(6,2)=temp_next_62;
      next(7,2)=temp_next_72;
    end
    path_m(8,k+4)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state
    
elseif(K==4&&(rate==1/3||rate==0.33||rate==0.34))
    for t=1:3:(3*(k+K-1))
      temp=[bsc(t) bsc(t+1) bsc(t+2)];                     % forming pairs of 2 bits from the sequence received from the channel 
      % S_1
      curr=min(next(1,1),next(2,1));        % next(1,1):from state-1 carrying 0
      path_m(8,(t+2)/3)=curr;                    % Finding path metrics of that state
      branch_m(15,(t+2)/3)=dhamming([1 1 1],temp);
      branch_m(16,(t+2)/3)=dhamming([0 0 0],temp);
      temp_next_12=branch_m(15,(t+2)/3)+curr;     %Finding branch metrics and adding to path metric of current state
      next(1,1)=branch_m(16,(t+2)/3)+curr;        %Finding branch metrics and adding to path metric of curent state

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(7,(t+2)/3)=curr;
      branch_m(13,(t+2)/3)=dhamming([0 0 0],temp);
      branch_m(14,(t+2)/3)=dhamming([1 1 1],temp);
      temp_next_22=branch_m(13,(t+2)/3)+curr;
      next(2,1)=branch_m(14,(t+2)/3)+curr;
      
      %s_3
      curr=min(next(5,1),next(6,1));
      path_m(6,(t+2)/3)=curr;
      branch_m(11,(t+2)/3)=dhamming([0 1 0],temp);
      branch_m(12,(t+2)/3)=dhamming([1 0 1],temp);
      temp_next_32=branch_m(11,(t+2)/3)+curr;
      next(3,1)=branch_m(12,(t+2)/3)+curr;
      
      %s_4
      curr=min(next(7,1),next(8,1));
      path_m(5,(t+2)/3)=curr;
      branch_m(9,(t+2)/3)=dhamming([1 0 1],temp);
      branch_m(10,(t+2)/3)=dhamming([0 1 0],temp);
      temp_next_42=branch_m(9,(t+2)/3)+curr;
      next(4,1)=branch_m(10,(t+2)/3)+curr;
      
       %s_5
      curr=min(next(2,2),next(1,2));
      path_m(4,(t+2)/3)=curr;
      branch_m(7,(t+2)/3)=dhamming([0 0 1],temp);
      branch_m(8,(t+2)/3)=dhamming([1 1 0],temp);
      temp_next_52=branch_m(7,(t+2)/3)+curr;
      next(5,1)=branch_m(8,(t+2)/3)+curr;
      
      %s_6
      curr=min(next(3,2),next(4,2));
      path_m(3,(t+2)/3)=curr;
      branch_m(5,(t+2)/3)=dhamming([1 1 0],temp);
      branch_m(6,(t+2)/3)=dhamming([0 0 1],temp);
      temp_next_62=branch_m(5,(t+2)/3)+curr;
      next(6,1)=branch_m(6,(t+2)/3)+curr;

      %s_7
      curr=min(next(5,2),next(6,2));
      path_m(2,(t+2)/3)=curr;
      branch_m(3,(t+2)/3)=dhamming([1 0 0],temp);
      branch_m(4,(t+2)/3)=dhamming([0 1 1],temp);
      temp_next_72=branch_m(3,(t+2)/3)+curr;
      next(7,1)=branch_m(4,(t+2)/3)+curr;

      %s_8
      curr=min(next(7,2),next(8,2));
      path_m(1,(t+2)/3)=curr;
      branch_m(1,(t+2)/3)=dhamming([0 1 1],temp);
      branch_m(2,(t+2)/3)=dhamming([1 0 0],temp);
      next(8,2)=branch_m(1,(t+2)/3)+curr;
      next(8,1)=branch_m(2,(t+2)/3)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
      next(4,2)=temp_next_42;
      next(5,2)=temp_next_52;
      next(6,2)=temp_next_62;
      next(7,2)=temp_next_72;
    end
    path_m(8,k+4)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state
end

% Traceback starting from s_1 i.e. [0 0] present in the last column
received=traceback(s,K,k,path_m,branch_m);
disp(received);      % Display the decoded sequence of BSC decoder

%The working of below two decoders is same as described above with minor changes.
% i.e. in case of decoder working for sequence received from Gaussian noise channel,
% soft-decision decoding(Eucledian distance) is considered instead of Hard-decision
%(Hamming Distance). BEC decoder is also quite similar to BSC decoder with few minor
% changes. Hence the explanatory comments are not repeated in further decoders.

% Decoder for the BEC Channel
[path_m,next,branch_m]=trellisinitialize(k,K);
if(K==3)
    % Loop that fills the path_m matrix made above with appropriate values of path metrics
    for t=1:2:(2*k+4)
      temp=[bec(t) bec(t+1)];             % forming pairs of 2 bits from the sequence received from the channel 
      % S_1
      curr=min(next(1,1),next(2,1));
      path_m(4,(t+1)/2)=curr;       % Finding path metrics of that state
      branch_m(7,(t+1)/2)=dbec([1,1],temp);
      branch_m(8,(t+1)/2)=dbec([0,0],temp);
      temp_next_12=dbec([1 1],temp)+curr;   %Finding branch metrics and adding to path metric of current state
      next(1,1)=dbec([0 0],temp)+curr;        %Finding branch metrics and adding to path metric of curent state

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(3,(t+1)/2)=curr;
      branch_m(5,(t+1)/2)=dbec([0,0],temp);
      branch_m(6,(t+1)/2)=dbec([1,1],temp);
      temp_next_22=dbec([0 0],temp)+curr;
      next(2,1)=dbec([1 1],temp)+curr;

      %s_3
      curr=min(next(2,2),next(1,2));
      path_m(2,(t+1)/2)=curr;
      branch_m(3,(t+1)/2)=dbec([0,1],temp);
      branch_m(4,(t+1)/2)=dbec([1,0],temp);
      temp_next_32=dbec([0 1],temp)+curr;
      next(3,1)=dbec([1 0],temp)+curr;

      %s_4
      curr=min(next(3,2),next(4,2));
      path_m(1,(t+1)/2)=curr;
      branch_m(1,(t+1)/2)=dbec([1,0],temp);
      branch_m(2,(t+1)/2)=dbec([0,1],temp);
      next(4,2)=dbec([1 0],temp)+curr;
      next(4,1)=dbec([0 1],temp)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
    end
    path_m(4,k+3)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state
       
elseif(K==4 && rate==0.5)
    for t=1:2:(2*k+6)
      temp=[bec(t) bec(t+1)];                     % forming pairs of 2 bits from the sequence received from the channel 
      % S_1
      curr=min(next(1,1),next(2,1));        % next(1,1):from state-1 carrying 0
      path_m(8,(t+1)/2)=curr;                    % Finding path metrics of that state
      branch_m(15,(t+1)/2)=dbec([1,1],temp);
      branch_m(16,(t+1)/2)=dbec([0,0],temp);
      temp_next_12=branch_m(15,(t+1)/2)+curr;     %Finding branch metrics and adding to path metric of current state
      next(1,1)=branch_m(16,(t+1)/2)+curr;        %Finding branch metrics and adding to path metric of curent state

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(7,(t+1)/2)=curr;
      branch_m(13,(t+1)/2)=dbec([0,0],temp);
      branch_m(14,(t+1)/2)=dbec([1,1],temp);
      temp_next_22=branch_m(13,(t+1)/2)+curr;
      next(2,1)=branch_m(14,(t+1)/2)+curr;
      
      %s_3
      curr=min(next(5,1),next(6,1));
      path_m(6,(t+1)/2)=curr;
      branch_m(11,(t+1)/2)=dbec([0,1],temp);
      branch_m(12,(t+1)/2)=dbec([1,0],temp);
      temp_next_32=branch_m(11,(t+1)/2)+curr;
      next(3,1)=branch_m(12,(t+1)/2)+curr;
      
      %s_4
      curr=min(next(7,1),next(8,1));
      path_m(5,(t+1)/2)=curr;
      branch_m(9,(t+1)/2)=dbec([1,0],temp);
      branch_m(10,(t+1)/2)=dbec([0,1],temp);
      temp_next_42=branch_m(9,(t+1)/2)+curr;
      next(4,1)=branch_m(10,(t+1)/2)+curr;
      
       %s_5
      curr=min(next(2,2),next(1,2));
      path_m(4,(t+1)/2)=curr;
      branch_m(7,(t+1)/2)=dbec([0,0],temp);
      branch_m(8,(t+1)/2)=dbec([1,1],temp);
      temp_next_52=branch_m(7,(t+1)/2)+curr;
      next(5,1)=branch_m(8,(t+1)/2)+curr;
      
      %s_6
      curr=min(next(3,2),next(4,2));
      path_m(3,(t+1)/2)=curr;
      branch_m(5,(t+1)/2)=dbec([1,1],temp);
      branch_m(6,(t+1)/2)=dbec([0,0],temp);
      temp_next_62=branch_m(5,(t+1)/2)+curr;
      next(6,1)=branch_m(6,(t+1)/2)+curr;

      %s_7
      curr=min(next(5,2),next(6,2));
      path_m(2,(t+1)/2)=curr;
      branch_m(3,(t+1)/2)=dbec([1,0],temp);
      branch_m(4,(t+1)/2)=dbec([0,1],temp);
      temp_next_72=branch_m(3,(t+1)/2)+curr;
      next(7,1)=branch_m(4,(t+1)/2)+curr;

      %s_8
      curr=min(next(7,2),next(8,2));
      path_m(1,(t+1)/2)=curr;
      branch_m(1,(t+1)/2)=dbec([0,1],temp);
      branch_m(2,(t+1)/2)=dbec([1,0],temp);
      next(8,2)=branch_m(1,(t+1)/2)+curr;
      next(8,1)=branch_m(2,(t+1)/2)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
      next(4,2)=temp_next_42;
      next(5,2)=temp_next_52;
      next(6,2)=temp_next_62;
      next(7,2)=temp_next_72;
    end
    path_m(8,k+4)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state
    
 elseif(K==4&&(rate==1/3||rate==0.33||rate==0.34))
    for t=1:3:(3*(k+K-1))
      temp=[bec(t) bec(t+1) bec(t+2)];                     % forming pairs of 2 bits from the sequence received from the channel 
      % S_1
      curr=min(next(1,1),next(2,1));        % next(1,1):from state-1 carrying 0
      path_m(8,(t+2)/3)=curr;                    % Finding path metrics of that state
      branch_m(15,(t+2)/3)=dbec2([1,1,1],temp);
      branch_m(16,(t+2)/3)=dbec2([0,0,0],temp);
      temp_next_12=branch_m(15,(t+2)/3)+curr;     %Finding branch metrics and adding to path metric of current state
      next(1,1)=branch_m(16,(t+2)/3)+curr;        %Finding branch metrics and adding to path metric of curent state

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(7,(t+2)/3)=curr;
      branch_m(13,(t+2)/3)=dbec2([0,0,0],temp);
      branch_m(14,(t+2)/3)=dbec2([1,1,1],temp);
      temp_next_22=branch_m(13,(t+2)/3)+curr;
      next(2,1)=branch_m(14,(t+2)/3)+curr;
      
      %s_3
      curr=min(next(5,1),next(6,1));
      path_m(6,(t+2)/3)=curr;
      branch_m(11,(t+2)/3)=dbec2([0,1,0],temp);
      branch_m(12,(t+2)/3)=dbec2([1,0,1],temp);
      temp_next_32=branch_m(11,(t+2)/3)+curr;
      next(3,1)=branch_m(12,(t+2)/3)+curr;
      
      %s_4
      curr=min(next(7,1),next(8,1));
      path_m(5,(t+2)/3)=curr;
      branch_m(9,(t+2)/3)=dbec2([1,0,1],temp);
      branch_m(10,(t+2)/3)=dbec2([0,1,0],temp);
      temp_next_42=branch_m(9,(t+2)/3)+curr;
      next(4,1)=branch_m(10,(t+2)/3)+curr;
      
       %s_5
      curr=min(next(2,2),next(1,2));
      path_m(4,(t+2)/3)=curr;
      branch_m(7,(t+2)/3)=dbec2([0,0,1],temp);
      branch_m(8,(t+2)/3)=dbec2([1,1,0],temp);
      temp_next_52=branch_m(7,(t+2)/3)+curr;
      next(5,1)=branch_m(8,(t+2)/3)+curr;
      
      %s_6
      curr=min(next(3,2),next(4,2));
      path_m(3,(t+2)/3)=curr;
      branch_m(5,(t+2)/3)=dbec2([1,1,0],temp);
      branch_m(6,(t+2)/3)=dbec2([0,0,1],temp);
      temp_next_62=branch_m(5,(t+2)/3)+curr;
      next(6,1)=branch_m(6,(t+2)/3)+curr;

      %s_7
      curr=min(next(5,2),next(6,2));
      path_m(2,(t+2)/3)=curr;
      branch_m(3,(t+2)/3)=dbec2([1,0,0],temp);
      branch_m(4,(t+2)/3)=dbec2([0,1,1],temp);
      temp_next_72=branch_m(3,(t+2)/3)+curr;
      next(7,1)=branch_m(4,(t+2)/3)+curr;

      %s_8
      curr=min(next(7,2),next(8,2));
      path_m(1,(t+2)/3)=curr;
      branch_m(1,(t+2)/3)=dbec2([0,1,1],temp);
      branch_m(2,(t+2)/3)=dbec2([1,0,0],temp);
      next(8,2)=branch_m(1,(t+2)/3)+curr;
      next(8,1)=branch_m(2,(t+2)/3)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
      next(4,2)=temp_next_42;
      next(5,2)=temp_next_52;
      next(6,2)=temp_next_62;
      next(7,2)=temp_next_72;
    end
    path_m(8,k+4)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state
end

% Traceback
received=traceback(s,K,k,path_m,branch_m);
disp(received)          % Display the decoded sequence of BEC decoder

% Decoder for the Gaussian Noise Channel
[path_m,next,branch_m]=trellisinitialize(k,K);

if(K==3)
    % Loop that fills the path_m matrix made above with appropriate values of path metrics
    for t=1:2:(2*k+4)
      temp=[gec(t) gec(t+1)];             % forming pairs of 2 bits from the sequence received from the channel 
      % S_1
      curr=min(next(1,1),next(2,1));
      path_m(4,(t+1)/2)=curr;       % Finding path metrics of that state
      branch_m(7,(t+1)/2)=deuclid([1,1],temp);
      branch_m(8,(t+1)/2)=deuclid([-1,-1],temp);
      temp_next_12=deuclid([1 1],temp)+curr;   %Finding branch metrics and adding to path metric of current state
      next(1,1)=deuclid([-1 -1],temp)+curr;        %Finding branch metrics and adding to path metric of curent state

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(3,(t+1)/2)=curr;
      branch_m(5,(t+1)/2)=deuclid([-1,-1],temp);
      branch_m(6,(t+1)/2)=deuclid([1,1],temp);
      temp_next_22=deuclid([-1 -1],temp)+curr;
      next(2,1)=deuclid([1 1],temp)+curr;

      %s_3
      curr=min(next(2,2),next(1,2));
      path_m(2,(t+1)/2)=curr;
      branch_m(3,(t+1)/2)=deuclid([-1,1],temp);
      branch_m(4,(t+1)/2)=deuclid([1,-1],temp);
      temp_next_32=deuclid([-1 1],temp)+curr;
      next(3,1)=deuclid([1 -1],temp)+curr;

      %s_4
      curr=min(next(3,2),next(4,2));
      path_m(1,(t+1)/2)=curr;
      branch_m(1,(t+1)/2)=deuclid([1,-1],temp);
      branch_m(2,(t+1)/2)=deuclid([-1,1],temp);
      next(4,2)=deuclid([1 -1],temp)+curr;
      next(4,1)=deuclid([-1 1],temp)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
    end
    path_m(4,k+3)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state
        
elseif(K==4 && rate==1/2)
    for t=1:2:(2*k+6)
      temp=[gec(t) gec(t+1)];                     % forming pairs of 2 bits from the sequence received from the channel 
      % S_1
      curr=min(next(1,1),next(2,1));        % next(1,1):from state-1 carrying 0
      path_m(8,(t+1)/2)=curr;                    % Finding path metrics of that state
      branch_m(15,(t+1)/2)=deuclid([1,1],temp);
      branch_m(16,(t+1)/2)=deuclid([-1,-1],temp);
      temp_next_12=branch_m(15,(t+1)/2)+curr;     %Finding branch metrics and adding to path metric of current state
      next(1,1)=branch_m(16,(t+1)/2)+curr;        %Finding branch metrics and adding to path metric of curent state

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(7,(t+1)/2)=curr;
      branch_m(13,(t+1)/2)=deuclid([-1,-1],temp);
      branch_m(14,(t+1)/2)=deuclid([1,1],temp);
      temp_next_22=branch_m(13,(t+1)/2)+curr;
      next(2,1)=branch_m(14,(t+1)/2)+curr;
      
      %s_3
      curr=min(next(5,1),next(6,1));
      path_m(6,(t+1)/2)=curr;
      branch_m(11,(t+1)/2)=deuclid([-1,1],temp);
      branch_m(12,(t+1)/2)=deuclid([1,-1],temp);
      temp_next_32=branch_m(11,(t+1)/2)+curr;
      next(3,1)=branch_m(12,(t+1)/2)+curr;
      
      %s_4
      curr=min(next(7,1),next(8,1));
      path_m(5,(t+1)/2)=curr;
      branch_m(9,(t+1)/2)=deuclid([1,-1],temp);
      branch_m(10,(t+1)/2)=deuclid([-1,1],temp);
      temp_next_42=branch_m(9,(t+1)/2)+curr;
      next(4,1)=branch_m(10,(t+1)/2)+curr;
      
       %s_5
      curr=min(next(2,2),next(1,2));
      path_m(4,(t+1)/2)=curr;
      branch_m(7,(t+1)/2)=deuclid([-1,-1],temp);
      branch_m(8,(t+1)/2)=deuclid([1,1],temp);
      temp_next_52=branch_m(7,(t+1)/2)+curr;
      next(5,1)=branch_m(8,(t+1)/2)+curr;
      
      %s_6
      curr=min(next(3,2),next(4,2));
      path_m(3,(t+1)/2)=curr;
      branch_m(5,(t+1)/2)=deuclid([1,1],temp);
      branch_m(6,(t+1)/2)=deuclid([-1,-1],temp);
      temp_next_62=branch_m(5,(t+1)/2)+curr;
      next(6,1)=branch_m(6,(t+1)/2)+curr;

      %s_7
      curr=min(next(5,2),next(6,2));
      path_m(2,(t+1)/2)=curr;
      branch_m(3,(t+1)/2)=deuclid([1,-1],temp);
      branch_m(4,(t+1)/2)=deuclid([-1,1],temp);
      temp_next_72=branch_m(3,(t+1)/2)+curr;
      next(7,1)=branch_m(4,(t+1)/2)+curr;

      %s_8
      curr=min(next(7,2),next(8,2));
      path_m(1,(t+1)/2)=curr;
      branch_m(1,(t+1)/2)=deuclid([-1,1],temp);
      branch_m(2,(t+1)/2)=deuclid([1,-1],temp);
      next(8,2)=branch_m(1,(t+1)/2)+curr;
      next(8,1)=branch_m(2,(t+1)/2)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
      next(4,2)=temp_next_42;
      next(5,2)=temp_next_52;
      next(6,2)=temp_next_62;
      next(7,2)=temp_next_72;
    end
    path_m(8,k+4)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state
    
    elseif(K==4 && (rate==1/3||rate==0.33||rate==0.34))
    for t=1:3:3*(k+K-1)
      temp=[gec(t) gec(t+1) gec(t+2)];                     % forming pairs of 2 bits from the sequence received from the channel 
      % S_1
      curr=min(next(1,1),next(2,1));        % next(1,1):from state-1 carrying 0
      path_m(8,(t+2)/3)=curr;                    % Finding path metrics of that state
      branch_m(15,(t+2)/3)=deuclid([1,1,1],temp);
      branch_m(16,(t+2)/3)=deuclid([-1,-1,-1],temp);
      temp_next_12=branch_m(15,(t+2)/3)+curr;     %Finding branch metrics and adding to path metric of current state
      next(1,1)=branch_m(16,(t+2)/3)+curr;        %Finding branch metrics and adding to path metric of curent state

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(7,(t+2)/3)=curr;
      branch_m(13,(t+2)/3)=deuclid([-1,-1,-1],temp);
      branch_m(14,(t+2)/3)=deuclid([1,1,1],temp);
      temp_next_22=branch_m(13,(t+2)/3)+curr;
      next(2,1)=branch_m(14,(t+2)/3)+curr;
      
      %s_3
      curr=min(next(5,1),next(6,1));
      path_m(6,(t+2)/3)=curr;
      branch_m(11,(t+2)/3)=deuclid([-1,1,-1],temp);
      branch_m(12,(t+2)/3)=deuclid([1,-1,1],temp);
      temp_next_32=branch_m(11,(t+2)/3)+curr;
      next(3,1)=branch_m(12,(t+2)/3)+curr;
      
      %s_4
      curr=min(next(7,1),next(8,1));
      path_m(5,(t+2)/3)=curr;
      branch_m(9,(t+2)/3)=deuclid([1,-1,1],temp);
      branch_m(10,(t+2)/3)=deuclid([-1,1,-1],temp);
      temp_next_42=branch_m(9,(t+2)/3)+curr;
      next(4,1)=branch_m(10,(t+2)/3)+curr;
      
       %s_5
      curr=min(next(2,2),next(1,2));
      path_m(4,(t+2)/3)=curr;
      branch_m(7,(t+2)/3)=deuclid([-1,-1,1],temp);
      branch_m(8,(t+2)/3)=deuclid([1,1,-1],temp);
      temp_next_52=branch_m(7,(t+2)/3)+curr;
      next(5,1)=branch_m(8,(t+2)/3)+curr;
      
      %s_6
      curr=min(next(3,2),next(4,2));
      path_m(3,(t+2)/3)=curr;
      branch_m(5,(t+2)/3)=deuclid([1,1,-1],temp);
      branch_m(6,(t+2)/3)=deuclid([-1,-1,1],temp);
      temp_next_62=branch_m(5,(t+2)/3)+curr;
      next(6,1)=branch_m(6,(t+2)/3)+curr;

      %s_7
      curr=min(next(5,2),next(6,2));
      path_m(2,(t+2)/3)=curr;
      branch_m(3,(t+2)/3)=deuclid([1,-1,-1],temp);
      branch_m(4,(t+2)/3)=deuclid([-1,1,1],temp);
      temp_next_72=branch_m(3,(t+2)/3)+curr;
      next(7,1)=branch_m(4,(t+2)/3)+curr;

      %s_8
      curr=min(next(7,2),next(8,2));
      path_m(1,(t+2)/3)=curr;
      branch_m(1,(t+2)/3)=deuclid([-1,1,1],temp);
      branch_m(2,(t+2)/3)=deuclid([1,-1,-1],temp);
      next(8,2)=branch_m(1,(t+2)/3)+curr;
      next(8,1)=branch_m(2,(t+2)/3)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
      next(4,2)=temp_next_42;
      next(5,2)=temp_next_52;
      next(6,2)=temp_next_62;
      next(7,2)=temp_next_72;
    end
    path_m(8,k+4)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state

end

% Traceback

received=traceback(s,K,k,path_m,branch_m);
disp(received);        % display the decoded sequence of Gaussian decoder
 
%% Monte-Carlo simulations:
k=500;
Nsim=20000;
ber_bsc=zeros(1,17);
ber_bec=zeros(1,17);
ber_gec=zeros(1,17);
SNR=zeros(1,17);
K=3;
rate=0.5;
for u=0:0.5:8       % loop that runs 20000 simulations for different values of SNR DB
    for z=1:Nsim
        or=randi(2,1,k)-1;
        c=Encoder(or,k,K,rate);

        % Fixing the SNR value for that iteration
        SNR(2*u+1)=u;                      
        SNRinwatt=10^(SNR(2*u+1)/10); 

        % when the encoded message passes through the Gaussian Noise channel
        gec=Gaussian_noise_channel(c,k,K,rate,SNRinwatt);

        % When the encoded message passes through the BSC channel
        bsc=bsc_channel(c,k,K,rate,SNRinwatt);

        % when the encoded message passes through the BEC channel
        bec=bec_channel(c,k,K,rate,SNRinwatt);

        % Decoder for the BSC Channel
        % Initializing some required values of the path_m
        [path_m,next,branch_m,s]=trellisinitialize(k,K);

        % Loop that fills the path_m matrix made above with appropriate values of path metrics
         for t=1:2:(2*k+4)
              temp=[bsc(t) bsc(t+1)];            
              % S_1
              curr=min(next(1,1),next(2,1));
              path_m(4,(t+1)/2)=curr;       
              branch_m(7,(t+1)/2)=dhamming([1,1],temp);
              branch_m(8,(t+1)/2)=dhamming([0,0],temp);
              temp_next_12=branch_m(7,(t+1)/2)+curr;   
              next(1,1)=branch_m(8,(t+1)/2)+curr;        

              %s_2
              curr=min(next(3,1),next(4,1));
              path_m(3,(t+1)/2)=curr;
              branch_m(5,(t+1)/2)=dhamming([0,0],temp);
              branch_m(6,(t+1)/2)=dhamming([1,1],temp);
              temp_next_22=branch_m(5,(t+1)/2)+curr;
              next(2,1)=branch_m(6,(t+1)/2)+curr;

              %s_3
              curr=min(next(2,2),next(1,2));
              path_m(2,(t+1)/2)=curr;
              branch_m(3,(t+1)/2)=dhamming([0,1],temp);
              branch_m(4,(t+1)/2)=dhamming([1,0],temp);
              temp_next_32=branch_m(3,(t+1)/2)+curr;
              next(3,1)=branch_m(4,(t+1)/2)+curr;

              %s_4
              curr=min(next(3,2),next(4,2));
              path_m(1,(t+1)/2)=curr;
              branch_m(1,(t+1)/2)=dhamming([1,0],temp);
              branch_m(2,(t+1)/2)=dhamming([0,1],temp);
              next(4,2)=branch_m(1,(t+1)/2)+curr;
              next(4,1)=branch_m(2,(t+1)/2)+curr;

              next(1,2)=temp_next_12;
              next(2,2)=temp_next_22;
              next(3,2)=temp_next_32;
        end
        path_m(4,k+3)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state

        % Traceback starting from s_1 i.e. [0 0] present in the last column
        received=traceback(s,K,k,path_m,branch_m);
        ber_bsc(2*u+1)=ber_bsc(2*u+1)+(dhamming([or zeros(1,K-1)],received)/500);      

        % Decoder for the BEC Channel
        [path_m,next,branch_m]=trellisinitialize(k,K);
        % Loop that fills the path_m matrix made above with appropriate values of path metrics
        for t=1:2:(2*k+4)
          temp=[bec(t) bec(t+1)];             
          % S_1
          curr=min(next(1,1),next(2,1));
          path_m(4,(t+1)/2)=curr;       
          branch_m(7,(t+1)/2)=dbec([1,1],temp);
          branch_m(8,(t+1)/2)=dbec([0,0],temp);
          temp_next_12=branch_m(7,(t+1)/2)+curr;   
          next(1,1)=branch_m(8,(t+1)/2)+curr;        

          %s_2
          curr=min(next(3,1),next(4,1));
          path_m(3,(t+1)/2)=curr;
          branch_m(5,(t+1)/2)=dbec([0,0],temp);
          branch_m(6,(t+1)/2)=dbec([1,1],temp);
          temp_next_22=branch_m(5,(t+1)/2)+curr;
          next(2,1)=branch_m(6,(t+1)/2)+curr;

          %s_3
          curr=min(next(2,2),next(1,2));
          path_m(2,(t+1)/2)=curr;
          branch_m(3,(t+1)/2)=dbec([0,1],temp);
          branch_m(4,(t+1)/2)=dbec([1,0],temp);
          temp_next_32=branch_m(3,(t+1)/2)+curr;
          next(3,1)=branch_m(4,(t+1)/2)+curr;

          %s_4
          curr=min(next(3,2),next(4,2));
          path_m(1,(t+1)/2)=curr;
          branch_m(1,(t+1)/2)=dbec([1,0],temp);
          branch_m(2,(t+1)/2)=dbec([0,1],temp);
          next(4,2)=branch_m(1,(t+1)/2)+curr;
          next(4,1)=branch_m(2,(t+1)/2)+curr;

          next(1,2)=temp_next_12;
          next(2,2)=temp_next_22;
          next(3,2)=temp_next_32;
        end
        path_m(4,k+3)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state

        % Traceback
        received=traceback(s,K,k,path_m,branch_m);
        ber_bec(2*u+1)=ber_bec(2*u+1)+(dhamming([or zeros(1,K-1)],received)/500);         

        % Decoder for the Gaussian Noise Channel
        [path_m,next,branch_m]=trellisinitialize(k,K);
        
        for t=1:2:(2*k+4)
          temp=[gec(t) gec(t+1)];             
          % S_1
          curr=min(next(1,1),next(2,1));
          path_m(4,(t+1)/2)=curr;       
          branch_m(7,(t+1)/2)=deuclid([1,1],temp);
          branch_m(8,(t+1)/2)=deuclid([-1,-1],temp);
          temp_next_12=branch_m(7,(t+1)/2)+curr;   
          next(1,1)=branch_m(8,(t+1)/2)+curr;       

          %s_2
          curr=min(next(3,1),next(4,1));
          path_m(3,(t+1)/2)=curr;
          branch_m(5,(t+1)/2)=deuclid([-1,-1],temp);
          branch_m(6,(t+1)/2)=deuclid([1,1],temp);
          temp_next_22=branch_m(5,(t+1)/2)+curr;
          next(2,1)=branch_m(6,(t+1)/2)+curr;

          %s_3
          curr=min(next(2,2),next(1,2));
          path_m(2,(t+1)/2)=curr;
          branch_m(3,(t+1)/2)=deuclid([-1,1],temp);
          branch_m(4,(t+1)/2)=deuclid([1,-1],temp);
          temp_next_32=branch_m(3,(t+1)/2)+curr;
          next(3,1)=branch_m(4,(t+1)/2)+curr;

          %s_4
          curr=min(next(3,2),next(4,2));
          path_m(1,(t+1)/2)=curr;
          branch_m(1,(t+1)/2)=deuclid([1,-1],temp);
          branch_m(2,(t+1)/2)=deuclid([-1,1],temp);
          next(4,2)=branch_m(1,(t+1)/2)+curr;
          next(4,1)=branch_m(2,(t+1)/2)+curr;

          next(1,2)=temp_next_12;
          next(2,2)=temp_next_22;
          next(3,2)=temp_next_32;
        end
        path_m(4,k+3)=min(next(1,1),next(2,1));         % last row last column of path_m i.e. 00 state

        % Traceback
        received=traceback(s,K,k,path_m,branch_m);
        ber_gec(2*u+1)=ber_gec(2*u+1)+(dhamming([or zeros(1,K-1)],received)/500);
    end
end
% Snippet to plot the obtained results
figure(1)
semilogy(SNR,ber_bsc/Nsim,'o-','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on; 
semilogy(SNR,ber_gec/Nsim,'^-','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
hold on; 
semilogy(SNR,ber_bec/Nsim,'d-','linewidth',2,'color',[0 0.4 0.9],'markerfacecolor',[0 0.4 0.9],'markeredgecolor',[0 0.4 0.9]);

xlabel('SNR per Bit in dB');
ylabel('Probability of Bit Error'); 
grid on;

legend('BSC','Gaussian Noise','BEC');
axis([0 8 1e-7 1]);
set(gca,'xtick',0:0.5:8);

pbErrBSC_Kequalto3=ber_bsc/Nsim;
pbErrAWGN_Kequalto3=ber_gec/Nsim;
pbErrBEC_Kequalto3=ber_bec/Nsim;



%% Intermediate-1 Results
ber_bsc=zeros(1,17);
ber_bec=zeros(1,17);
ber_gec=zeros(1,17);
SNR=zeros(1,17);
K=4;
rate=0.5;
for u=0:0.5:8       % loop that runs 20000 simulations for different values of SNR DB
    for z=1:Nsim
     disp(z);
    disp(u);% for debugging,to check time complexity
    or=randi(2,1,k)-1;
    c=Encoder(or,k,K,rate);
    
    % Fixing the SNR value for that iteration
    SNR(2*u+1)=u;                      
    SNRinwatt=10^(SNR(2*u+1)/10); 
    
    % when the encoded message passes through the Gaussian Noise channel
    gec=Gaussian_noise_channel(c,k,K,rate,SNRinwatt);

    % When the encoded message passes through the BSC channel
    bsc=bsc_channel(c,k,K,rate,SNRinwatt);

    % when the encoded message passes through the BEC channel
    bec=bec_channel(c,k,K,rate,SNRinwatt);
    
    % Decoder for the BSC Channel
    % Initializing some required values of the path_m
    [path_m,next,branch_m,s]=trellisinitialize(k,K);
    
    % Loop that fills the path_m matrix made above with appropriate values of path metrics
     for t=1:2:(2*k+6)
      temp=[bsc(t) bsc(t+1)];                     % forming pairs of 2 bits from the sequence received from the channel 
      % S_1
      curr=min(next(1,1),next(2,1));        % next(1,1):from state-1 carrying 0
      path_m(8,(t+1)/2)=curr;                    % Finding path metrics of that state
      branch_m(15,(t+1)/2)=dhamming([1,1],temp);
      branch_m(16,(t+1)/2)=dhamming([0,0],temp);
      temp_next_12=branch_m(15,(t+1)/2)+curr;     %Finding branch metrics and adding to path metric of current state
      next(1,1)=branch_m(16,(t+1)/2)+curr;        %Finding branch metrics and adding to path metric of curent state

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(7,(t+1)/2)=curr;
      branch_m(13,(t+1)/2)=dhamming([0,0],temp);
      branch_m(14,(t+1)/2)=dhamming([1,1],temp);
      temp_next_22=branch_m(13,(t+1)/2)+curr;
      next(2,1)=branch_m(14,(t+1)/2)+curr;
      
      %s_3
      curr=min(next(5,1),next(6,1));
      path_m(6,(t+1)/2)=curr;
      branch_m(11,(t+1)/2)=dhamming([0,1],temp);
      branch_m(12,(t+1)/2)=dhamming([1,0],temp);
      temp_next_32=branch_m(11,(t+1)/2)+curr;
      next(3,1)=branch_m(12,(t+1)/2)+curr;
      
      %s_4
      curr=min(next(7,1),next(8,1));
      path_m(5,(t+1)/2)=curr;
      branch_m(9,(t+1)/2)=dhamming([1,0],temp);
      branch_m(10,(t+1)/2)=dhamming([0,1],temp);
      temp_next_42=branch_m(9,(t+1)/2)+curr;
      next(4,1)=branch_m(10,(t+1)/2)+curr;
      
       %s_5
      curr=min(next(2,2),next(1,2));
      path_m(4,(t+1)/2)=curr;
      branch_m(7,(t+1)/2)=dhamming([0,0],temp);
      branch_m(8,(t+1)/2)=dhamming([1,1],temp);
      temp_next_52=branch_m(7,(t+1)/2)+curr;
      next(5,1)=branch_m(8,(t+1)/2)+curr;
      
      %s_6
      curr=min(next(3,2),next(4,2));
      path_m(3,(t+1)/2)=curr;
      branch_m(5,(t+1)/2)=dhamming([1,1],temp);
      branch_m(6,(t+1)/2)=dhamming([0,0],temp);
      temp_next_62=branch_m(5,(t+1)/2)+curr;
      next(6,1)=branch_m(6,(t+1)/2)+curr;

      %s_7
      curr=min(next(5,2),next(6,2));
      path_m(2,(t+1)/2)=curr;
      branch_m(3,(t+1)/2)=dhamming([1,0],temp);
      branch_m(4,(t+1)/2)=dhamming([0,1],temp);
      temp_next_72=branch_m(3,(t+1)/2)+curr;
      next(7,1)=branch_m(4,(t+1)/2)+curr;

      %s_8
      curr=min(next(7,2),next(8,2));
      path_m(1,(t+1)/2)=curr;
      branch_m(1,(t+1)/2)=dhamming([0,1],temp);
      branch_m(2,(t+1)/2)=dhamming([1,0],temp);
      next(8,2)=branch_m(1,(t+1)/2)+curr;
      next(8,1)=branch_m(2,(t+1)/2)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
      next(4,2)=temp_next_42;
      next(5,2)=temp_next_52;
      next(6,2)=temp_next_62;
      next(7,2)=temp_next_72;
    end
    path_m(8,k+4)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state

    % Traceback starting from s_1 i.e. [0 0] present in the last column
    received=traceback(s,K,k,path_m,branch_m);
    ber_bsc(2*u+1)=ber_bsc(2*u+1)+(dhamming([or zeros(1,K-1)],received)/500);      

    % Decoder for the BEC Channel
    [path_m,next,branch_m]=trellisinitialize(k,K);
    % Loop that fills the path_m matrix made above with appropriate values of path metrics
    for t=1:2:(2*k+6)
      temp=[bec(t) bec(t+1)];                     % forming pairs of 2 bits from the sequence received from the channel 
      % S_1
      curr=min(next(1,1),next(2,1));        % next(1,1):from state-1 carrying 0
      path_m(8,(t+1)/2)=curr;                    % Finding path metrics of that state
      branch_m(15,(t+1)/2)=dbec([1,1],temp);
      branch_m(16,(t+1)/2)=dbec([0,0],temp);
      temp_next_12=branch_m(15,(t+1)/2)+curr;     %Finding branch metrics and adding to path metric of current state
      next(1,1)=branch_m(16,(t+1)/2)+curr;        %Finding branch metrics and adding to path metric of curent state

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(7,(t+1)/2)=curr;
      branch_m(13,(t+1)/2)=dbec([0,0],temp);
      branch_m(14,(t+1)/2)=dbec([1,1],temp);
      temp_next_22=branch_m(13,(t+1)/2)+curr;
      next(2,1)=branch_m(14,(t+1)/2)+curr;
      
      %s_3
      curr=min(next(5,1),next(6,1));
      path_m(6,(t+1)/2)=curr;
      branch_m(11,(t+1)/2)=dbec([0,1],temp);
      branch_m(12,(t+1)/2)=dbec([1,0],temp);
      temp_next_32=branch_m(11,(t+1)/2)+curr;
      next(3,1)=branch_m(12,(t+1)/2)+curr;
      
      %s_4
      curr=min(next(7,1),next(8,1));
      path_m(5,(t+1)/2)=curr;
      branch_m(9,(t+1)/2)=dbec([1,0],temp);
      branch_m(10,(t+1)/2)=dbec([0,1],temp);
      temp_next_42=branch_m(9,(t+1)/2)+curr;
      next(4,1)=branch_m(10,(t+1)/2)+curr;
      
       %s_5
      curr=min(next(2,2),next(1,2));
      path_m(4,(t+1)/2)=curr;
      branch_m(7,(t+1)/2)=dbec([0,0],temp);
      branch_m(8,(t+1)/2)=dbec([1,1],temp);
      temp_next_52=branch_m(7,(t+1)/2)+curr;
      next(5,1)=branch_m(8,(t+1)/2)+curr;
      
      %s_6
      curr=min(next(3,2),next(4,2));
      path_m(3,(t+1)/2)=curr;
      branch_m(5,(t+1)/2)=dbec([1,1],temp);
      branch_m(6,(t+1)/2)=dbec([0,0],temp);
      temp_next_62=branch_m(5,(t+1)/2)+curr;
      next(6,1)=branch_m(6,(t+1)/2)+curr;

      %s_7
      curr=min(next(5,2),next(6,2));
      path_m(2,(t+1)/2)=curr;
      branch_m(3,(t+1)/2)=dbec([1,0],temp);
      branch_m(4,(t+1)/2)=dbec([0,1],temp);
      temp_next_72=branch_m(3,(t+1)/2)+curr;
      next(7,1)=branch_m(4,(t+1)/2)+curr;

      %s_8
      curr=min(next(7,2),next(8,2));
      path_m(1,(t+1)/2)=curr;
      branch_m(1,(t+1)/2)=dbec([0,1],temp);
      branch_m(2,(t+1)/2)=dbec([1,0],temp);
      next(8,2)=branch_m(1,(t+1)/2)+curr;
      next(8,1)=branch_m(2,(t+1)/2)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
      next(4,2)=temp_next_42;
      next(5,2)=temp_next_52;
      next(6,2)=temp_next_62;
      next(7,2)=temp_next_72;
    end
    path_m(8,k+4)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state

    % Traceback
    received=traceback(s,K,k,path_m,branch_m);
    ber_bec(2*u+1)=ber_bec(2*u+1)+(dhamming([or zeros(1,K-1)],received)/500);         

    % Decoder for the Gaussian Noise Channel
    [path_m,next,branch_m]=trellisinitialize(k,K);
    % Loop that fills the path_m matrix made above with appropriate values of path metrics
   for t=1:2:(2*k+6)
      temp=[gec(t) gec(t+1)];                     % forming pairs of 2 bits from the sequence received from the channel 
      % S_1
      curr=min(next(1,1),next(2,1));        % next(1,1):from state-1 carrying 0
      path_m(8,(t+1)/2)=curr;                    % Finding path metrics of that state
      branch_m(15,(t+1)/2)=deuclid([1,1],temp);
      branch_m(16,(t+1)/2)=deuclid([-1,-1],temp);
      temp_next_12=branch_m(15,(t+1)/2)+curr;     %Finding branch metrics and adding to path metric of current state
      next(1,1)=branch_m(16,(t+1)/2)+curr;        %Finding branch metrics and adding to path metric of curent state

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(7,(t+1)/2)=curr;
      branch_m(13,(t+1)/2)=deuclid([-1,-1],temp);
      branch_m(14,(t+1)/2)=deuclid([1,1],temp);
      temp_next_22=branch_m(13,(t+1)/2)+curr;
      next(2,1)=branch_m(14,(t+1)/2)+curr;
      
      %s_3
      curr=min(next(5,1),next(6,1));
      path_m(6,(t+1)/2)=curr;
      branch_m(11,(t+1)/2)=deuclid([-1,1],temp);
      branch_m(12,(t+1)/2)=deuclid([1,-1],temp);
      temp_next_32=branch_m(11,(t+1)/2)+curr;
      next(3,1)=branch_m(12,(t+1)/2)+curr;
      
      %s_4
      curr=min(next(7,1),next(8,1));
      path_m(5,(t+1)/2)=curr;
      branch_m(9,(t+1)/2)=deuclid([1,-1],temp);
      branch_m(10,(t+1)/2)=deuclid([-1,1],temp);
      temp_next_42=branch_m(9,(t+1)/2)+curr;
      next(4,1)=branch_m(10,(t+1)/2)+curr;
      
       %s_5
      curr=min(next(2,2),next(1,2));
      path_m(4,(t+1)/2)=curr;
      branch_m(7,(t+1)/2)=deuclid([-1,-1],temp);
      branch_m(8,(t+1)/2)=deuclid([1,1],temp);
      temp_next_52=branch_m(7,(t+1)/2)+curr;
      next(5,1)=branch_m(8,(t+1)/2)+curr;
      
      %s_6
      curr=min(next(3,2),next(4,2));
      path_m(3,(t+1)/2)=curr;
      branch_m(5,(t+1)/2)=deuclid([1,1],temp);
      branch_m(6,(t+1)/2)=deuclid([-1,-1],temp);
      temp_next_62=branch_m(5,(t+1)/2)+curr;
      next(6,1)=branch_m(6,(t+1)/2)+curr;

      %s_7
      curr=min(next(5,2),next(6,2));
      path_m(2,(t+1)/2)=curr;
      branch_m(3,(t+1)/2)=deuclid([1,-1],temp);
      branch_m(4,(t+1)/2)=deuclid([-1,1],temp);
      temp_next_72=branch_m(3,(t+1)/2)+curr;
      next(7,1)=branch_m(4,(t+1)/2)+curr;

      %s_8
      curr=min(next(7,2),next(8,2));
      path_m(1,(t+1)/2)=curr;
      branch_m(1,(t+1)/2)=deuclid([-1,1],temp);
      branch_m(2,(t+1)/2)=deuclid([1,-1],temp);
      next(8,2)=branch_m(1,(t+1)/2)+curr;
      next(8,1)=branch_m(2,(t+1)/2)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
      next(4,2)=temp_next_42;
      next(5,2)=temp_next_52;
      next(6,2)=temp_next_62;
      next(7,2)=temp_next_72;
    end
    path_m(8,k+4)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state
    
    % Traceback
    received=traceback(s,K,k,path_m,branch_m);
    ber_gec(2*u+1)=ber_gec(2*u+1)+(dhamming([or zeros(1,K-1)],received)/500);
    end
end
pbErrBSC_Kequalto4=ber_bsc/Nsim;
pbErrBEC_Kequalto4=ber_bec/Nsim;
pbErrAWGN_Kequalto4=ber_gec/Nsim;
figure(2)
semilogy(SNR,pbErrBSC_Kequalto3,'o-','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on;
semilogy(SNR,pbErrBSC_Kequalto4,'o:','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on; 
semilogy(SNR,pbErrAWGN_Kequalto3,'^-','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
hold on;
semilogy(SNR,pbErrAWGN_Kequalto4,'^:','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
hold on; 
semilogy(SNR,pbErrBEC_Kequalto3,'d-','linewidth',2,'color',[0 0.4 0.9],'markerfacecolor',[0 0.4 0.9],'markeredgecolor',[0 0.4 0.9]);
hold on;
semilogy(SNR,pbErrBEC_Kequalto4,'d:','linewidth',2,'color',[0 0.4 0.9],'markerfacecolor',[0 0.4 0.9],'markeredgecolor',[0 0.4 0.9]);

xlabel('SNR per Bit in dB');
ylabel('Probability of Bit Error'); 
grid on;

legend('BSC K=3','BSC K=4','AWGN K=3','AWGN K=4','BEC K=3','BEC K=4');
axis([0 8 1e-7 1]);
set(gca,'xtick',0:0.5:8);


%% Intermediate 2 Results
ber_bsc=zeros(1,17);
ber_bec=zeros(1,17);
ber_gec=zeros(1,17);
SNR=zeros(1,17);
K=4;
rate=1/3;
for u=0:0.5:8       % loop that runs 20000 simulations for different values of SNR DB
    for z=1:Nsim
     disp(z);
    disp(u);% for debugging,to check time complexity
    or=randi(2,1,k)-1;
    c=Encoder(or,k,K,rate);
    
    % Fixing the SNR value for that iteration
    SNR(2*u+1)=u;                      
    SNRinwatt=10^(SNR(2*u+1)/10); 
    
    % when the encoded message passes through the Gaussian Noise channel
    gec=Gaussian_noise_channel(c,k,K,rate,SNRinwatt);

    % When the encoded message passes through the BSC channel
    bsc=bsc_channel(c,k,K,rate,SNRinwatt);

    % when the encoded message passes through the BEC channel
    bec=bec_channel(c,k,K,rate,SNRinwatt);
    
    % Decoder for the BSC Channel
    % Initializing some required values of the path_m
    [path_m,next,branch_m,s]=trellisinitialize(k,K);
    
    % Loop that fills the path_m matrix made above with appropriate values of path metrics
   for t=1:3:(3*(k+K-1))
      temp=[bsc(t) bsc(t+1) bsc(t+2)];                     % forming pairs of 2 bits from the sequence received from the channel 
      % S_1
      curr=min(next(1,1),next(2,1));        % next(1,1):from state-1 carrying 0
      path_m(8,(t+2)/3)=curr;                    % Finding path metrics of that state
      branch_m(15,(t+2)/3)=dhamming([1 1 1],temp);
      branch_m(16,(t+2)/3)=dhamming([0 0 0],temp);
      temp_next_12=branch_m(15,(t+2)/3)+curr;     %Finding branch metrics and adding to path metric of current state
      next(1,1)=branch_m(16,(t+2)/3)+curr;        %Finding branch metrics and adding to path metric of curent state

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(7,(t+2)/3)=curr;
      branch_m(13,(t+2)/3)=dhamming([0 0 0],temp);
      branch_m(14,(t+2)/3)=dhamming([1 1 1],temp);
      temp_next_22=branch_m(13,(t+2)/3)+curr;
      next(2,1)=branch_m(14,(t+2)/3)+curr;
      
      %s_3
      curr=min(next(5,1),next(6,1));
      path_m(6,(t+2)/3)=curr;
      branch_m(11,(t+2)/3)=dhamming([0 1 0],temp);
      branch_m(12,(t+2)/3)=dhamming([1 0 1],temp);
      temp_next_32=branch_m(11,(t+2)/3)+curr;
      next(3,1)=branch_m(12,(t+2)/3)+curr;
      
      %s_4
      curr=min(next(7,1),next(8,1));
      path_m(5,(t+2)/3)=curr;
      branch_m(9,(t+2)/3)=dhamming([1 0 1],temp);
      branch_m(10,(t+2)/3)=dhamming([0 1 0],temp);
      temp_next_42=branch_m(9,(t+2)/3)+curr;
      next(4,1)=branch_m(10,(t+2)/3)+curr;
      
       %s_5
      curr=min(next(2,2),next(1,2));
      path_m(4,(t+2)/3)=curr;
      branch_m(7,(t+2)/3)=dhamming([0 0 1],temp);
      branch_m(8,(t+2)/3)=dhamming([1 1 0],temp);
      temp_next_52=branch_m(7,(t+2)/3)+curr;
      next(5,1)=branch_m(8,(t+2)/3)+curr;
      
      %s_6
      curr=min(next(3,2),next(4,2));
      path_m(3,(t+2)/3)=curr;
      branch_m(5,(t+2)/3)=dhamming([1 1 0],temp);
      branch_m(6,(t+2)/3)=dhamming([0 0 1],temp);
      temp_next_62=branch_m(5,(t+2)/3)+curr;
      next(6,1)=branch_m(6,(t+2)/3)+curr;

      %s_7
      curr=min(next(5,2),next(6,2));
      path_m(2,(t+2)/3)=curr;
      branch_m(3,(t+2)/3)=dhamming([1 0 0],temp);
      branch_m(4,(t+2)/3)=dhamming([0 1 1],temp);
      temp_next_72=branch_m(3,(t+2)/3)+curr;
      next(7,1)=branch_m(4,(t+2)/3)+curr;

      %s_8
      curr=min(next(7,2),next(8,2));
      path_m(1,(t+2)/3)=curr;
      branch_m(1,(t+2)/3)=dhamming([0 1 1],temp);
      branch_m(2,(t+2)/3)=dhamming([1 0 0],temp);
      next(8,2)=branch_m(1,(t+2)/3)+curr;
      next(8,1)=branch_m(2,(t+2)/3)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
      next(4,2)=temp_next_42;
      next(5,2)=temp_next_52;
      next(6,2)=temp_next_62;
      next(7,2)=temp_next_72;
    end
    path_m(8,k+4)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state
    
    % Traceback starting from s_1 i.e. [0 0] present in the last column
    received=traceback(s,K,k,path_m,branch_m);
    ber_bsc(2*u+1)=ber_bsc(2*u+1)+(dhamming([or zeros(1,K-1)],received)/500);      

    % Decoder for the BEC Channel
    [path_m,next,branch_m]=trellisinitialize(k,K);
    % Loop that fills the path_m matrix made above with appropriate values of path metrics
    
    for t=1:3:(3*(k+K-1))
      temp=[bec(t) bec(t+1) bec(t+2)];                     % forming pairs of 2 bits from the sequence received from the channel 
      % S_1
      curr=min(next(1,1),next(2,1));        % next(1,1):from state-1 carrying 0
      path_m(8,(t+2)/3)=curr;                    % Finding path metrics of that state
      branch_m(15,(t+2)/3)=dbec2([1,1,1],temp);
      branch_m(16,(t+2)/3)=dbec2([0,0,0],temp);
      temp_next_12=branch_m(15,(t+2)/3)+curr;     %Finding branch metrics and adding to path metric of current state
      next(1,1)=branch_m(16,(t+2)/3)+curr;        %Finding branch metrics and adding to path metric of curent state

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(7,(t+2)/3)=curr;
      branch_m(13,(t+2)/3)=dbec2([0,0,0],temp);
      branch_m(14,(t+2)/3)=dbec2([1,1,1],temp);
      temp_next_22=branch_m(13,(t+2)/3)+curr;
      next(2,1)=branch_m(14,(t+2)/3)+curr;
      
      %s_3
      curr=min(next(5,1),next(6,1));
      path_m(6,(t+2)/3)=curr;
      branch_m(11,(t+2)/3)=dbec2([0,1,0],temp);
      branch_m(12,(t+2)/3)=dbec2([1,0,1],temp);
      temp_next_32=branch_m(11,(t+2)/3)+curr;
      next(3,1)=branch_m(12,(t+2)/3)+curr;
      
      %s_4
      curr=min(next(7,1),next(8,1));
      path_m(5,(t+2)/3)=curr;
      branch_m(9,(t+2)/3)=dbec2([1,0,1],temp);
      branch_m(10,(t+2)/3)=dbec2([0,1,0],temp);
      temp_next_42=branch_m(9,(t+2)/3)+curr;
      next(4,1)=branch_m(10,(t+2)/3)+curr;
      
       %s_5
      curr=min(next(2,2),next(1,2));
      path_m(4,(t+2)/3)=curr;
      branch_m(7,(t+2)/3)=dbec2([0,0,1],temp);
      branch_m(8,(t+2)/3)=dbec2([1,1,0],temp);
      temp_next_52=branch_m(7,(t+2)/3)+curr;
      next(5,1)=branch_m(8,(t+2)/3)+curr;
      
      %s_6
      curr=min(next(3,2),next(4,2));
      path_m(3,(t+2)/3)=curr;
      branch_m(5,(t+2)/3)=dbec2([1,1,0],temp);
      branch_m(6,(t+2)/3)=dbec2([0,0,1],temp);
      temp_next_62=branch_m(5,(t+2)/3)+curr;
      next(6,1)=branch_m(6,(t+2)/3)+curr;

      %s_7
      curr=min(next(5,2),next(6,2));
      path_m(2,(t+2)/3)=curr;
      branch_m(3,(t+2)/3)=dbec2([1,0,0],temp);
      branch_m(4,(t+2)/3)=dbec2([0,1,1],temp);
      temp_next_72=branch_m(3,(t+2)/3)+curr;
      next(7,1)=branch_m(4,(t+2)/3)+curr;

      %s_8
      curr=min(next(7,2),next(8,2));
      path_m(1,(t+2)/3)=curr;
      branch_m(1,(t+2)/3)=dbec2([0,1,1],temp);
      branch_m(2,(t+2)/3)=dbec2([1,0,0],temp);
      next(8,2)=branch_m(1,(t+2)/3)+curr;
      next(8,1)=branch_m(2,(t+2)/3)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
      next(4,2)=temp_next_42;
      next(5,2)=temp_next_52;
      next(6,2)=temp_next_62;
      next(7,2)=temp_next_72;
    end
    path_m(8,k+4)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state

    % Traceback
    received=traceback(s,K,k,path_m,branch_m);
    ber_bec(2*u+1)=ber_bec(2*u+1)+(dhamming([or zeros(1,K-1)],received)/500);         

    % Decoder for the Gaussian Noise Channel
       [path_m,next,branch_m]=trellisinitialize(k,K);
    % Loop that fills the path_m matrix made above with appropriate values of path metrics
   for t=1:3:3*(k+K-1)
      temp=[gec(t) gec(t+1) gec(t+2)];                     % forming pairs of 2 bits from the sequence received from the channel 
      % S_1
      curr=min(next(1,1),next(2,1));        % next(1,1):from state-1 carrying 0
      path_m(8,(t+2)/3)=curr;                    % Finding path metrics of that state
      branch_m(15,(t+2)/3)=deuclid([1,1,1],temp);
      branch_m(16,(t+2)/3)=deuclid([-1,-1,-1],temp);
      temp_next_12=branch_m(15,(t+2)/3)+curr;     %Finding branch metrics and adding to path metric of current state
      next(1,1)=branch_m(16,(t+2)/3)+curr;        %Finding branch metrics and adding to path metric of curent state

      %s_2
      curr=min(next(3,1),next(4,1));
      path_m(7,(t+2)/3)=curr;
      branch_m(13,(t+2)/3)=deuclid([-1,-1,-1],temp);
      branch_m(14,(t+2)/3)=deuclid([1,1,1],temp);
      temp_next_22=branch_m(13,(t+2)/3)+curr;
      next(2,1)=branch_m(14,(t+2)/3)+curr;
      
      %s_3
      curr=min(next(5,1),next(6,1));
      path_m(6,(t+2)/3)=curr;
      branch_m(11,(t+2)/3)=deuclid([-1,1,-1],temp);
      branch_m(12,(t+2)/3)=deuclid([1,-1,1],temp);
      temp_next_32=branch_m(11,(t+2)/3)+curr;
      next(3,1)=branch_m(12,(t+2)/3)+curr;
      
      %s_4
      curr=min(next(7,1),next(8,1));
      path_m(5,(t+2)/3)=curr;
      branch_m(9,(t+2)/3)=deuclid([1,-1,1],temp);
      branch_m(10,(t+2)/3)=deuclid([-1,1,-1],temp);
      temp_next_42=branch_m(9,(t+2)/3)+curr;
      next(4,1)=branch_m(10,(t+2)/3)+curr;
      
       %s_5
      curr=min(next(2,2),next(1,2));
      path_m(4,(t+2)/3)=curr;
      branch_m(7,(t+2)/3)=deuclid([-1,-1,1],temp);
      branch_m(8,(t+2)/3)=deuclid([1,1,-1],temp);
      temp_next_52=branch_m(7,(t+2)/3)+curr;
      next(5,1)=branch_m(8,(t+2)/3)+curr;
      
      %s_6
      curr=min(next(3,2),next(4,2));
      path_m(3,(t+2)/3)=curr;
      branch_m(5,(t+2)/3)=deuclid([1,1,-1],temp);
      branch_m(6,(t+2)/3)=deuclid([-1,-1,1],temp);
      temp_next_62=branch_m(5,(t+2)/3)+curr;
      next(6,1)=branch_m(6,(t+2)/3)+curr;

      %s_7
      curr=min(next(5,2),next(6,2));
      path_m(2,(t+2)/3)=curr;
      branch_m(3,(t+2)/3)=deuclid([1,-1,-1],temp);
      branch_m(4,(t+2)/3)=deuclid([-1,1,1],temp);
      temp_next_72=branch_m(3,(t+2)/3)+curr;
      next(7,1)=branch_m(4,(t+2)/3)+curr;

      %s_8
      curr=min(next(7,2),next(8,2));
      path_m(1,(t+2)/3)=curr;
      branch_m(1,(t+2)/3)=deuclid([-1,1,1],temp);
      branch_m(2,(t+2)/3)=deuclid([1,-1,-1],temp);
      next(8,2)=branch_m(1,(t+2)/3)+curr;
      next(8,1)=branch_m(2,(t+2)/3)+curr;

      next(1,2)=temp_next_12;
      next(2,2)=temp_next_22;
      next(3,2)=temp_next_32;
      next(4,2)=temp_next_42;
      next(5,2)=temp_next_52;
      next(6,2)=temp_next_62;
      next(7,2)=temp_next_72;
    end
    path_m(8,k+4)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state
    
    % Traceback
    received=traceback(s,K,k,path_m,branch_m);
    ber_gec(2*u+1)=ber_gec(2*u+1)+(dhamming([or zeros(1,K-1)],received)/500);
    end
end
pbErrBSC_reqtoinverseof3=ber_bsc/Nsim;
pbErrBEC_reqtoinverseof3=ber_bec/Nsim;
pbErrAWGN_reqtoinverseof3=ber_gec/Nsim;
figure(3)
semilogy(SNR,pbErrBSC_Kequalto3,'o-','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on;
semilogy(SNR,pbErrBSC_reqtoinverseof3,'o:','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on; 
semilogy(SNR,pbErrAWGN_Kequalto3,'^-','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
hold on;
semilogy(SNR,pbErrAWGN_reqtoinverseof3,'^:','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
hold on; 
semilogy(SNR,pbErrBEC_Kequalto3,'d-','linewidth',2,'color',[0 0.4 0.9],'markerfacecolor',[0 0.4 0.9],'markeredgecolor',[0 0.4 0.9]);
hold on;
semilogy(SNR,pbErrBEC_reqtoinverseof3,'d:','linewidth',2,'color',[0 0.4 0.9],'markerfacecolor',[0 0.4 0.9],'markeredgecolor',[0 0.4 0.9]);

xlabel('SNR per Bit in dB');
ylabel('Probability of Bit Error'); 
grid on;

legend('BSC K=3','BSC K=4 r=0.33','AWGN K=3','AWGN K=4 r=0.33','BEC K=3','BEC K=4 r=0.33');
axis([0 8 1e-7 1]);
set(gca,'xtick',0:0.5:8);
