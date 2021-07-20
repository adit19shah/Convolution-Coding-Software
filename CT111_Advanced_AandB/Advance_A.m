% CT111 Final Project(Advanced)   Developed by: Adit Shah   Topic: Convolution coding
% To clear all the variables,clear command window and close all previous
% graphs
clearvars;
close all;
clc; 
%% Advanced A Results for K=3 and rate=1/2
close all;
clearvars;
k=500;
Nsim=20000;
ber_bsc=zeros(1,17);
ber_bec=zeros(1,17);
ber_gec=zeros(1,17);
SNR=zeros(1,17);
SNRinwatt=zeros(1,17);
K=3;
rate=0.5;
for u=0:0.5:8       % loop that runs 20000 simulations for different values of SNR DB
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

        % When the encoded message passes through the BSC channel
        bsc=bsc_channel(c,k,K,rate,SNRinwatt(2*u+1));

        % when the encoded message passes through the BEC channel
        bec=bec_channel(c,k,K,rate,SNRinwatt(2*u+1));

        % Decoder for the BSC Channel
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
        ber_bsc(2*u+1)=ber_bsc(2*u+1)+(dhamming([or zeros(1,K-1)],received)/500);             

        % Decoder for the Gaussian Noise Channel
        [path_m,next,branch_m,s]=trellisinitialize(k,K);
        % Loop that fills the path_m matrix made above with appropriate values of path metrics
        for t=1:2:(2*k+4)
          temp=[gec(t) gec(t+1)];             % forming pairs of 2 bits from the sequence received from the channel 
          % S_1
          curr=min(next(1,1),next(2,1));
          path_m(4,(t+1)/2)=curr;       % Finding path metrics of that state
          branch_m(7,(t+1)/2)=deuclid([1,1],temp);
          branch_m(8,(t+1)/2)=deuclid([-1,-1],temp);
          temp_next_12=branch_m(7,(t+1)/2)+curr;   %Finding branch metrics and adding to path metric of current state
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
        ber_gec(2*u+1)=ber_gec(2*u+1)+(dhamming([or zeros(1,K-1)],received)/500);
    end
end
% Snippet to plot the obtained results
pbErrBSC_Kequalto3=ber_bsc/Nsim;
pbErrBEC_Kequalto3=ber_bec/Nsim;
pbErrAWGN_Kequalto3=ber_gec/Nsim;
p=qfunc(sqrt(2*rate*SNRinwatt));
%upper_bound=(D.^5)./((1-2*N*D).^2);
%BSC Channel
D=(sqrt(4*p.*(1-p)));
N=1;
upper_bound_BSC=(D.^5)./((1-2*N*D).^2);

%Gaussian Noise Channel
D=exp(-rate*SNRinwatt);
N=1;
upper_bound_AWGN=(D.^5)./((1-2*N*D).^2);
figure(1)
semilogy(SNR,ber_bsc/Nsim,'o-','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on; 
semilogy(SNR,ber_gec/Nsim,'^-','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
hold on;
semilogy(SNR,upper_bound_BSC,'o:','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on;
semilogy(SNR,upper_bound_AWGN,'^:','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);

xlabel('SNR per Bit in dB');
ylabel('Probability of Bit Error'); 
grid on;

legend('BSC K=3','Gaussian Noise K=3','UpperBoundBSC','UpperBoundAWGN');
axis([0 8 1e-7 1e5]);
set(gca,'xtick',0:0.5:8);


%% Advanced-A Results for K=4 and rate=1/2
ber_bsc=zeros(1,17);
ber_bec=zeros(1,17);
ber_gec=zeros(1,17);
SNR=zeros(1,17);
SNRinwatt=zeros(1,17);
K=4;
rate=0.5;
Nsim=10000;
for u=0:0.5:8       % loop that runs 20000 simulations for different values of SNR DB
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

    % When the encoded message passes through the BSC channel
    bsc=bsc_channel(c,k,K,rate,SNRinwatt(2*u+1));

    % when the encoded message passes through the BEC channel
    bec=bec_channel(c,k,K,rate,SNRinwatt(2*u+1));
    
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

    % Decoder for the Gaussian Noise Channel
    [path_m,next,branch_m,s]=trellisinitialize(k,K);
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
p=qfunc(sqrt(2*rate*SNRinwatt));
%BSC Channel
D=(sqrt(4*p.*(1-p)));
N=1;
upper_bound_BSC=((D.^7)+2*(N*D.^6)-2*(N*D.^8)-(2*(N^2)*D.^7)+((N^2)*D.^11)+((N^2)*D.^9))./((1-2*N*D-(N*D).^3).^2);

%Gaussian Noise Channel
D=exp(-rate*SNRinwatt);
N=1;
upper_bound_AWGN=((D.^7)+2*(N*D.^6)-2*(N*D.^8)-(2*(N^2)*D.^7)+((N^2)*D.^11)+((N^2)*D.^9))./((1-2*N*D-(N*D).^3).^2);

figure(2)
semilogy(SNR,pbErrBSC_Kequalto4,'o-','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on;
semilogy(SNR,upper_bound_BSC,'o:','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on; 
semilogy(SNR,pbErrAWGN_Kequalto4,'^-','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
hold on;
semilogy(SNR,upper_bound_AWGN,'^:','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
hold on; 

xlabel('SNR per Bit in dB');
ylabel('Probability of Bit Error'); 
grid on;

legend('BSC K=4','BSCUpperBounds','Gaussian Noise K=4','AWGNUpperBOunds');
axis([0 8 1e-7 1]);
set(gca,'xtick',0:0.5:8);
