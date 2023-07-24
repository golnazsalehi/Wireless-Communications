nc = 1e5; 
fft_size = 1e5;
OFDM_num = 1e3;
SNR = 10.^((-10:20)/10);
a = 1;
N0 = a^2./SNR;
Pe_vect_MRC = zeros(1,length(SNR));
Tx_num = 10;
for i = 1:length(SNR)
  Pe_per_frame = 0;
  for k = 1:OFDM_num
      L = 200;
      x = randi([0 1],1,nc);
      x_hat_vect = zeros(Tx_num,nc);
      h = sqrt(0.5)*(randn(Tx_num,L)+1i*randn(Tx_num,L));
      for j=1:Tx_num
        y = Tx(x, L, h(j,:), N0(i), fft_size);
        x_hat_vect(j,:) = Rx_MRC( L, y, h(j,:) ,N0(i), fft_size);
      end
      z = sum(x_hat_vect);
      x_hat = real(z)>0;
      Pe_per_frame = Pe_per_frame + sum(x~=x_hat);
  end
  Pe_vect_MRC(i) = Pe_per_frame/(nc*OFDM_num);

end
figure(1)
semilogy(-10:20,Pe_vect_MRC,'LineWidth',2);
hold on
grid on
xlabel('SNR(dB)', 'Interpreter', 'latex');
ylabel('BER', 'Interpreter', 'latex');
title('MRC','Interpreter', 'latex');
%%
function y = Tx(x, L, h, N0, fft_size)
    coding = [-1 1];
    x_encoded = coding(x+1);
    X = ifft(x_encoded,fft_size);
    X = [X(length(X)-L+1+1:length(X)) X];
    noise = sqrt(N0/2).*(randn(1,(length(X) + length(h)-1))+1i*randn(1,(length(X) + length(h)-1)));
    y = conv(X,h) + noise;
end
function z = Rx_MRC( L, Y,h, N0, fft_size)
    y = Y(L:fft_size+L-1);
    H = fft([h zeros(1,length(y)-length(h))]);
    yfft = fft(y,fft_size);
    z = yfft.* abs(H).*exp(-1i*angle(H))./sqrt(N0);
end
