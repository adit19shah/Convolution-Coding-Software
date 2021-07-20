function plotresults(SNR,ber_bsc,ber_gec,ber_bec,Nsim)
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


