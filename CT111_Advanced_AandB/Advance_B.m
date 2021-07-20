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
m=randi(2,1,k)-1;
disp(m);
c=Encoder(m,k,K,rate);  % Encodes the randomly generated message m


% Fixing the SNR ratio
SNR=5;                    % SNR in db
SNRinwatt=10^(SNR/10);    % converting SNR from db to watt

%probability of bit getting flipped/erased
p=qfunc(sqrt(2*rate*SNRinwatt));            

 % when the encoded message passes through the Gaussian Noise channel
gec=Gaussian_noise_channel(c,k,K,rate,SNRinwatt);


% When the encoded message passes through the BSC channel
bsc=bsc_channel(c,k,K,rate,SNRinwatt);

% when the encoded message passes through the BEC channel
bec=bec_channel(c,k,K,rate,SNRinwatt);

% Decoder for the Gaussian Noise Channel
[path_m,next,branch_m]=trellisinitialize(k,K);

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
    
% Decoding using BCJR Algorithm for Gaussian Noise Channel
alpha=zeros(size(path_m));
beta=zeros(size(path_m));
gamma=ones(size(branch_m));
    for i=1:k+2
        temp=[gec(2*i-1) gec(2*i)];
        gamma(8,i)=0.5*exp(-SNRinwatt*deuclid(temp,[-1,-1]));
        gamma(7,i)=0.5*exp(-SNRinwatt*deuclid(temp,[1,1]));
        gamma(6,i)=0.5*exp(-SNRinwatt*deuclid(temp,[1,1]));
        gamma(5,i)=0.5*exp(-SNRinwatt*deuclid(temp,[-1,-1]));
        gamma(4,i)=0.5*exp(-SNRinwatt*deuclid(temp,[1,-1]));
        gamma(3,i)=0.5*exp(-SNRinwatt*deuclid(temp,[-1,1]));
        gamma(2,i)=0.5*exp(-SNRinwatt*deuclid(temp,[-1,1]));
        gamma(1,i)=0.5*exp(-SNRinwatt*deuclid(temp,[1,-1]));
    end
    alpha=calculateAlpha(k+3,gamma,alpha);
    beta=calculateBeta(1,gamma,beta,k);
app=zeros(1,k+2);
app_num=zeros(1,k+2);
app_den=zeros(1,k+2);
for i=1:k+2
    app_num(i)=(beta(1,i+1)*gamma(1,i)*alpha(1,i))+(beta(1,i+1)*gamma(3,i)*alpha(2,i))+(beta(2,i+1)*gamma(5,i)*alpha(3,i))+(beta(2,i+1)*gamma(7,i)*alpha(4,i));
    app_den(i)=(beta(3,i+1)*gamma(2,i)*alpha(1,i))+(beta(3,i+1)*gamma(4,i)*alpha(2,i))+(beta(4,i+1)*gamma(6,i)*alpha(3,i))+(beta(4,i+1)*gamma(8,i)*alpha(4,i)); 
    app(i)=log(app_num(i)/app_den(i));
end
bcjr_output=-ones(1,k+2);
bcjr_output(logical(app>=0))=1;
bcjr_output(logical(app<0))=0;
%disp(bcjr_output);
dhamming(bcjr_output,[m 0 0])

% Decoding using BCJR Algorithm for BSC Channel
alpha=zeros(size(path_m));
beta=zeros(size(path_m));
gamma=ones(size(branch_m));
    for i=1:k+2
        temp=[bsc(2*i-1) bsc(2*i)];
        gamma(8,i)=p^(dhamming(temp,[0 0]))*(1-p)^((1/rate)-dhamming(temp,[0 0]));
        gamma(7,i)=p^(dhamming(temp,[1 1]))*(1-p)^((1/rate)-dhamming(temp,[1 1]));
        gamma(6,i)=p^(dhamming(temp,[1 1]))*(1-p)^((1/rate)-dhamming(temp,[1 1]));
        gamma(5,i)=p^(dhamming(temp,[0 0]))*(1-p)^((1/rate)-dhamming(temp,[0 0]));
        gamma(4,i)=p^(dhamming(temp,[1 0]))*(1-p)^((1/rate)-dhamming(temp,[1 0]));
        gamma(3,i)=p^(dhamming(temp,[0 1]))*(1-p)^((1/rate)-dhamming(temp,[0 1]));
        gamma(2,i)=p^(dhamming(temp,[0 1]))*(1-p)^((1/rate)-dhamming(temp,[0 1]));
        gamma(1,i)=p^(dhamming(temp,[1 0]))*(1-p)^((1/rate)-dhamming(temp,[1 0]));
    end
    alpha=calculateAlpha(k+3,gamma,alpha);
    beta=calculateBeta(1,gamma,beta,k);
app=zeros(1,k+2);
app_num=zeros(1,k+2);
app_den=zeros(1,k+2);
for i=1:k+2
    app_num(i)=(beta(1,i+1)*gamma(1,i)*alpha(1,i))+(beta(1,i+1)*gamma(3,i)*alpha(2,i))+(beta(2,i+1)*gamma(5,i)*alpha(3,i))+(beta(2,i+1)*gamma(7,i)*alpha(4,i));
    app_den(i)=(beta(3,i+1)*gamma(2,i)*alpha(1,i))+(beta(3,i+1)*gamma(4,i)*alpha(2,i))+(beta(4,i+1)*gamma(6,i)*alpha(3,i))+(beta(4,i+1)*gamma(8,i)*alpha(4,i)); 
    app(i)=log(app_num(i)/app_den(i));
end
bcjr_output=-ones(1,k+2);
bcjr_output(logical(app>=0))=1;
bcjr_output(logical(app<0))=0;
%disp(bcjr_output);
dhamming(bcjr_output,[m 0 0])

% Decoding using BCJR Algorithm for BEC Channel
alpha=zeros(size(path_m));
beta=zeros(size(path_m));
gamma=ones(size(branch_m));
    for i=1:k+2
        temp=[bec(2*i-1) bec(2*i)];
        gamma(8,i)=p^(dhamming(temp,[0 0]))*(1-p)^((1/rate)-dhamming(temp,[0 0]));
        gamma(7,i)=p^(dhamming(temp,[1 1]))*(1-p)^((1/rate)-dhamming(temp,[1 1]));
        gamma(6,i)=p^(dhamming(temp,[1 1]))*(1-p)^((1/rate)-dhamming(temp,[1 1]));
        gamma(5,i)=p^(dhamming(temp,[0 0]))*(1-p)^((1/rate)-dhamming(temp,[0 0]));
        gamma(4,i)=p^(dhamming(temp,[1 0]))*(1-p)^((1/rate)-dhamming(temp,[1 0]));
        gamma(3,i)=p^(dhamming(temp,[0 1]))*(1-p)^((1/rate)-dhamming(temp,[0 1]));
        gamma(2,i)=p^(dhamming(temp,[0 1]))*(1-p)^((1/rate)-dhamming(temp,[0 1]));
        gamma(1,i)=p^(dhamming(temp,[1 0]))*(1-p)^((1/rate)-dhamming(temp,[1 0]));
    end
    alpha=calculateAlpha(k+3,gamma,alpha);
    beta=calculateBeta(1,gamma,beta,k);
app=zeros(1,k+2);
app_num=zeros(1,k+2);
app_den=zeros(1,k+2);
for i=1:k+2
    app_num(i)=(beta(1,i+1)*gamma(1,i)*alpha(1,i))+(beta(1,i+1)*gamma(3,i)*alpha(2,i))+(beta(2,i+1)*gamma(5,i)*alpha(3,i))+(beta(2,i+1)*gamma(7,i)*alpha(4,i));
    app_den(i)=(beta(3,i+1)*gamma(2,i)*alpha(1,i))+(beta(3,i+1)*gamma(4,i)*alpha(2,i))+(beta(4,i+1)*gamma(6,i)*alpha(3,i))+(beta(4,i+1)*gamma(8,i)*alpha(4,i)); 
    app(i)=log(app_num(i)/app_den(i));
end
bcjr_output=-ones(1,k+2);
bcjr_output(logical(app>=0))=1;
bcjr_output(logical(app<0))=0;
%disp(bcjr_output);
dhamming(bcjr_output,[m 0 0])
%% Monte Carlo Simulations
close all;
k=50;
Nsim=20000;
ber_gec_viterbi=zeros(1,17);
ber_gec_bcjr=zeros(1,17);
SNR=zeros(1,17);
SNRinwatt=zeros(1,17);
K=3;
rate=0.5;
for u=0:0.5:8      % loop that runs 20000 simulations for different values of SNR DB
    for z=1:Nsim
        disp(z);
        disp(u);% for debugging,to check time complexity
        or=randi(2,1,k)-1;
        c=Encoder(or,k,K,rate);

        % Fixing the SNR value for that iteration
        SNR(2*u+1)=u;                      
        SNRinwatt(2*u+1)=10^(SNR(2*u+1)/10); 

        % when the encoded message passes through the Gaussian Noise channel
        gec=Gaussian_noise_channel(c,k,K,rate,SNRinwatt(2*u+1));
        
        % Decoder for the Gaussian Noise Channel
        [path_m,next,branch_m,s]=trellisinitialize(k,K);
        % Loop that fills the path_m matrix made above with appropriate values of path metrics
        for t=1:2:(2*k+4)
          temp=[gec(t) gec(t+1)];                   % forming pairs of 2 bits from the sequence received from the channel 
          % S_1
          curr=min(next(1,1),next(2,1));
          path_m(4,(t+1)/2)=curr;                    % Finding path metrics of that state
          branch_m(7,(t+1)/2)=deuclid([1,1],temp);
          branch_m(8,(t+1)/2)=deuclid([-1,-1],temp);
          temp_next_12=branch_m(7,(t+1)/2)+curr;     %Finding branch metrics and adding to path metric of current state
          next(1,1)=branch_m(8,(t+1)/2)+curr;        %Finding branch metrics and adding to path metric of curent state

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
        path_m(4,k+3)=min(next(1,1),next(2,1));  % last row last column of path_m i.e. 00 state

        % Traceback
        received=traceback(s,K,k,path_m,branch_m);
        ber_gec_viterbi(2*u+1)=ber_gec_viterbi(2*u+1)+dhamming([or zeros(1,K-1)],received); 
        
        % BCJR Decoding for AWGN Channel
        alpha=zeros(size(path_m));
        beta=zeros(size(path_m));
        gamma=ones(size(branch_m));
            for i=1:k+2
                temp=[gec(2*i-1) gec(2*i)];
            gamma(8,i)=0.5*(sqrt(SNRinwatt(2*u+1)/pi))*exp(-SNRinwatt(2*u+1)*deuclid(temp,[-1,-1]));
            gamma(7,i)=0.5*(sqrt(SNRinwatt(2*u+1)/pi))*exp(-SNRinwatt(2*u+1)*deuclid(temp,[1,1]));
            gamma(6,i)=0.5*(sqrt(SNRinwatt(2*u+1)/pi))*exp(-SNRinwatt(2*u+1)*deuclid(temp,[1,1]));
            gamma(5,i)=0.5*(sqrt(SNRinwatt(2*u+1)/pi))*exp(-SNRinwatt(2*u+1)*deuclid(temp,[-1,-1]));
            gamma(4,i)=0.5*(sqrt(SNRinwatt(2*u+1)/pi))*exp(-SNRinwatt(2*u+1)*deuclid(temp,[1,-1]));
            gamma(3,i)=0.5*(sqrt(SNRinwatt(2*u+1)/pi))*exp(-SNRinwatt(2*u+1)*deuclid(temp,[-1,1]));
            gamma(2,i)=0.5*(sqrt(SNRinwatt(2*u+1)/pi))*exp(-SNRinwatt(2*u+1)*deuclid(temp,[-1,1]));
            gamma(1,i)=0.5*(sqrt(SNRinwatt(2*u+1)/pi))*exp(-SNRinwatt(2*u+1)*deuclid(temp,[1,-1]));
            end
        alpha=calculateAlpha(k+3,gamma,alpha);
        beta=calculateBeta(1,gamma,beta,k);
        app=zeros(1,k+2);
        app_num=zeros(1,k+2);
        app_den=zeros(1,k+2);
        for i=1:k+2
            app_num(i)=(beta(1,i+1)*gamma(1,i)*alpha(1,i))+(beta(1,i+1)*gamma(3,i)*alpha(2,i))+(beta(2,i+1)*gamma(5,i)*alpha(3,i))+(beta(2,i+1)*gamma(7,i)*alpha(4,i));
            app_den(i)=(beta(3,i+1)*gamma(2,i)*alpha(1,i))+(beta(3,i+1)*gamma(4,i)*alpha(2,i))+(beta(4,i+1)*gamma(6,i)*alpha(3,i))+(beta(4,i+1)*gamma(8,i)*alpha(4,i)); 
            app(i)=log(app_num(i)/app_den(i));
        end
        bcjr_output=-ones(1,k+2);
        bcjr_output(logical(app>=0))=1;
        bcjr_output(logical(app<0))=0; 
        ber_gec_bcjr(2*u+1)=ber_gec_bcjr(2*u+1)+dhamming([or zeros(1,K-1)],bcjr_output); 
    end
end
pbErrAWGN_viterbi=ber_gec_viterbi/(Nsim*k);
pbErrAWGN_bcjr=ber_gec_bcjr/(Nsim*k);
semilogy(SNR,pbErrAWGN_viterbi,'^-','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
hold on;
semilogy(SNR,pbErrAWGN_bcjr,'--o','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
xlabel('SNR per Bit in dB');
ylabel('Probability of Bit Error'); 
grid on;

legend('Viterbi Decoder for AWGN Channel','BCJR Decoder for AWGN Channel');
axis([0 8 1e-7 1]);
set(gca,'xtick',0:0.5:8);

%%
k=500;
Nsim=20000;
ber_bec_viterbi=zeros(1,17);
ber_bec_bcjr=zeros(1,17);
ber_bsc_viterbi=zeros(1,17);
ber_bsc_bcjr=zeros(1,17);
SNR=zeros(1,17);
SNRinwatt=zeros(1,17);
K=3;
rate=0.5;
for u=0:0.5:8      % loop that runs 20000 simulations for different values of SNR DB
    for z=1:Nsim
        disp(z);
        disp(u);% for debugging,to check time complexity
        or=randi(2,1,k)-1;
        c=Encoder(or,k,K,rate);

        % Fixing the SNR value for that iteration
        SNR(2*u+1)=u;                      
        SNRinwatt(2*u+1)=10^(SNR(2*u+1)/10); 
        
        %probability of bit getting flipped/erased
        p=qfunc(sqrt(2*rate*SNRinwatt));

        % When the encoded message passes through the BSC channel
        bsc=bsc_channel(c,k,K,rate,SNRinwatt(2*u+1));

        % when the encoded message passes through the BEC channel
        bec=bec_channel(c,k,K,rate,SNRinwatt(2*u+1));
        
        % Viterbi Decoding for BSC
        % Initializing some required values of the path_m
        [path_m,next,branch_m,s]=trellisinitialize(k,K);

        % Loop that fills the path_m matrix made above with appropriate values of path metrics
         for t=1:2:(2*k+4)
              temp=[bsc(t) bsc(t+1)];             % forming pairs of 2 bits from the sequence received from the channel 
              % S_1
              curr=min(next(1,1),next(2,1));
              path_m(4,(t+1)/2)=curr;       % Finding path metrics of that state
              branch_m(7,(t+1)/2)=dhamming([1,1],temp);
              branch_m(8,(t+1)/2)=dhamming([0,0],temp);
              temp_next_12=branch_m(7,(t+1)/2)+curr;   %Finding branch metrics and adding to path metric of current state
              next(1,1)=branch_m(8,(t+1)/2)+curr;        %Finding branch metrics and adding to path metric of curent state

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
        ber_bsc_viterbi(2*u+1)=ber_bsc_viterbi(2*u+1)+dhamming([or zeros(1,K-1)],received);      
        
% Decoding using BCJR Algorithm for BSC Channel
        alpha=zeros(size(path_m));
        beta=zeros(size(path_m));
        gamma=ones(size(branch_m));
            for i=1:k+2
                temp=[bsc(2*i-1) bsc(2*i)];
                gamma(8,i)=p(2*u+1)^(dhamming(temp,[0 0]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[0 0]));
                gamma(7,i)=p(2*u+1)^(dhamming(temp,[1 1]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[1 1]));
                gamma(6,i)=p(2*u+1)^(dhamming(temp,[1 1]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[1 1]));
                gamma(5,i)=p(2*u+1)^(dhamming(temp,[0 0]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[0 0]));
                gamma(4,i)=p(2*u+1)^(dhamming(temp,[1 0]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[1 0]));
                gamma(3,i)=p(2*u+1)^(dhamming(temp,[0 1]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[0 1]));
                gamma(2,i)=p(2*u+1)^(dhamming(temp,[0 1]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[0 1]));
                gamma(1,i)=p(2*u+1)^(dhamming(temp,[1 0]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[1 0]));
            end
            alpha=calculateAlpha(k+3,gamma,alpha);
            beta=calculateBeta(1,gamma,beta,k);
        app=zeros(1,k+2);
        app_num=zeros(1,k+2);
        app_den=zeros(1,k+2);
        for i=1:k+2
            app_num(i)=(beta(1,i+1)*gamma(1,i)*alpha(1,i))+(beta(1,i+1)*gamma(3,i)*alpha(2,i))+(beta(2,i+1)*gamma(5,i)*alpha(3,i))+(beta(2,i+1)*gamma(7,i)*alpha(4,i));
            app_den(i)=(beta(3,i+1)*gamma(2,i)*alpha(1,i))+(beta(3,i+1)*gamma(4,i)*alpha(2,i))+(beta(4,i+1)*gamma(6,i)*alpha(3,i))+(beta(4,i+1)*gamma(8,i)*alpha(4,i)); 
            app(i)=log(app_num(i)/app_den(i));
        end
        bcjr_output=-ones(1,k+2);
        bcjr_output(logical(app>=0))=1;
        bcjr_output(logical(app<0))=0;
        ber_bsc_bcjr(2*u+1)=ber_bsc_bcjr(2*u+1)+dhamming([or zeros(1,K-1)],bcjr_output); 
%         
        % Viterbi Decoding for BEC Channel
        [path_m,next,branch_m]=trellisinitialize(k,K);
        % Loop that fills the path_m matrix made above with appropriate values of path metrics
        for t=1:2:(2*k+4)
          temp=[bec(t) bec(t+1)];             % forming pairs of 2 bits from the sequence received from the channel 
          % S_1
          curr=min(next(1,1),next(2,1));
          path_m(4,(t+1)/2)=curr;       % Finding path metrics of that state
          branch_m(7,(t+1)/2)=dbec([1,1],temp);
          branch_m(8,(t+1)/2)=dbec([0,0],temp);
          temp_next_12=branch_m(7,(t+1)/2)+curr;   %Finding branch metrics and adding to path metric of current state
          next(1,1)=branch_m(8,(t+1)/2)+curr;        %Finding branch metrics and adding to path metric of curent state

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
        ber_bec_viterbi(2*u+1)=ber_bec_viterbi(2*u+1)+dhamming([or zeros(1,K-1)],received); 
        
        % Decoding using BCJR Algorithm for BSC Channel
        alpha=zeros(size(path_m));
        beta=zeros(size(path_m));
        gamma=ones(size(branch_m));
            for i=1:k+2
                temp=[bec(2*i-1) bec(2*i)];
                gamma(8,i)=p(2*u+1)^(dhamming(temp,[0 0]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[0 0]));
                gamma(7,i)=p(2*u+1)^(dhamming(temp,[1 1]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[1 1]));
                gamma(6,i)=p(2*u+1)^(dhamming(temp,[1 1]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[1 1]));
                gamma(5,i)=p(2*u+1)^(dhamming(temp,[0 0]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[0 0]));
                gamma(4,i)=p(2*u+1)^(dhamming(temp,[1 0]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[1 0]));
                gamma(3,i)=p(2*u+1)^(dhamming(temp,[0 1]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[0 1]));
                gamma(2,i)=p(2*u+1)^(dhamming(temp,[0 1]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[0 1]));
                gamma(1,i)=p(2*u+1)^(dhamming(temp,[1 0]))*(1-p(2*u+1))^((1/rate)-dhamming(temp,[1 0]));
            end
            alpha=calculateAlpha(k+3,gamma,alpha);
            beta=calculateBeta(1,gamma,beta,k);
        app=zeros(1,k+2);
        app_num=zeros(1,k+2);
        app_den=zeros(1,k+2);
        for i=1:k+2
            app_num(i)=(beta(1,i+1)*gamma(1,i)*alpha(1,i))+(beta(1,i+1)*gamma(3,i)*alpha(2,i))+(beta(2,i+1)*gamma(5,i)*alpha(3,i))+(beta(2,i+1)*gamma(7,i)*alpha(4,i));
            app_den(i)=(beta(3,i+1)*gamma(2,i)*alpha(1,i))+(beta(3,i+1)*gamma(4,i)*alpha(2,i))+(beta(4,i+1)*gamma(6,i)*alpha(3,i))+(beta(4,i+1)*gamma(8,i)*alpha(4,i)); 
            app(i)=log(app_num(i)/app_den(i));
        end
        bcjr_output=-ones(1,k+2);
        bcjr_output(logical(app>=0))=1;
        bcjr_output(logical(app<0))=0;
        ber_bec_bcjr(2*u+1)=ber_bec_bcjr(2*u+1)+dhamming([or zeros(1,K-1)],bcjr_output);

    end
end
% Plotting the results for BEC and BSC Channels
pbErrBSC_viterbi=ber_bsc_viterbi/(Nsim*k);
pbErrBSC_bcjr=ber_bsc_bcjr/(Nsim*k);
pbErrBEC_viterbi=ber_bec_viterbi/(Nsim*k);
pbErrBEC_bcjr=ber_bec_bcjr/(Nsim*k);
semilogy(SNR,pbErrBSC_viterbi,'o-','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on;
semilogy(SNR,pbErrBSC_bcjr,'o:','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on; 
semilogy(SNR,pbErrBEC_viterbi,'d-','linewidth',2,'color',[0 0.4 0.9],'markerfacecolor',[0 0.4 0.9],'markeredgecolor',[0 0.4 0.9]);
hold on;
semilogy(SNR,pbErrBEC_bcjr,'d:','linewidth',2,'color',[0 0.4 0.9],'markerfacecolor',[0 0.4 0.9],'markeredgecolor',[0 0.4 0.9]);

xlabel('SNR per Bit in dB');
ylabel('Probability of Bit Error'); 
grid on;
legend('BSC Viterbi Decoding','BSC BCJR Decoding','BEC Viterbi Decoding','BEC BCJR Decoding');
axis([0 8 1e-7 1]);
set(gca,'xtick',0:0.5:8);

       
