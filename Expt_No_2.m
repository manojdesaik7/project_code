clc;clear;close all;
x = 1; % signal to transmit Eb = 1
TRIAL = 10000; %number of simulation runs per EbN0 %50000
for EbN0 = 0:1:20 %dB
    linear_EbN0 = 10^(EbN0/10); 
    nvar = 1/(linear_EbN0); %calculation of N0, Here, Eb = 1
    error1 = 0; %set error counter to 0 for MRC
    error2 = 0; %set error counter to 0 for EGC
    error3 = 0; %set error counter to 0 for SC
    
        for trial = 1:TRIAL % monte carlo trials.. count the errors
            n1 = sqrt(nvar/2)*randn; %noise for the first RX1
            n2 = sqrt(nvar/2)*randn; %noise for the second RX2
            h1 = sqrt(0.5)*abs(randn + j*randn); %rayleigh amplitude RX1
            h2 = sqrt(0.5)*abs(randn + j*randn); %rayleigh amplitude RX2
            
            %Equal Gain combining
            y1 = x*h1+n1; % Signal 1
            y2 = x*h2+n2; % Signal 2
            y_equal = 0.5*(y1+y2); %constant gain of 0.5 
            
            %Maximal Ratio combining
            a1 = (abs(h1))^2;
            a2 = (abs(h2))^2;
            y_maximal = x*(a1*h1+a2*h2)+a1*n1+a2*n2;
            
         %Selection Combining
          if a1>=a2
              y_SC=x*h1*a1+a1*n1;
          else
              y_SC=x*h2*a2+a2*n2;
          end
          
            if y_equal < 0 %define decision region as 0 
                error1 = error1 + 1;
            end
            if y_maximal < 0 
                error2 = error2 + 1;
            end
            if y_SC < 0 
                error3 = error2 + 1;
            end
        end
    BER1(EbN0+1) = error1/(TRIAL);
    BER2(EbN0+1) = error2/(TRIAL);
    BER3(EbN0+1) = error3/(TRIAL);
end

% plot simulations

EbNo=0:1:20; %changed from 10
mu = 10.^(EbNo./10);
ber_theory = (1/2)*(1 - sqrt(mu ./ (mu + 1)));%Theoritical expression 
semilogy(EbNo,BER1,'r-.*');
hold on;
semilogy(EbNo,BER2,'g:s');
semilogy(EbNo,BER3,'b-o');
semilogy(EbNo,ber_theory,'y','LineWidth',2); % plot EG BER vs EbNo 
legend('EGC','MRC','SC','THEORY');
xlabel('SNR(dB)') %Label for x-axis
ylabel('BER') %Label for y-axis
title('SNR V/S BER');