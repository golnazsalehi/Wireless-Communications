%% ZF
nc = 1e5; 
fft_size = 1e5;
OFDM_num = 1e3;
SNR = 10.^((-10:20)/10);
a = 1;
N0 = a^2./SNR;
Pe_vect_ZF = zeros(1,length(SNR));

for i = 1:length(SNR)
  Pe_per_frame = 0;
  for k = 1:OFDM_num
    L = 200;
    h = sqrt(0.5)*(randn(1,L)+1i*randn(1,L));
    x = randi([0 1],1,nc);
    y = Tx(x, L, N0(i), h, fft_size);
    x_hat = Rx_ZF( L, y, h, fft_size);
    Pe_per_frame = Pe_per_frame + sum(x~=x_hat);
  end
  Pe_vect_ZF(i) = Pe_per_frame/(nc*OFDM_num);

end
figure(3)
semilogy(-10:20,Pe_vect_ZF,'LineWidth',2);
grid on
xlabel('SNR(dB)', 'Interpreter', 'latex');
ylabel('BER', 'Interpreter', 'latex');
title('Zero Forcing','Interpreter', 'latex');
%% MMSE
nc = 1e5; 
fft_size = 1e5;
OFDM_num = 1e3;
SNR = 10.^((-10:20)/10);
a = 1;
N0 = a^2./SNR;
Pe_vect_MMSE = zeros(1,length(SNR));

for i = 1:length(SNR)
  Pe_per_frame = 0;
  for k = 1:OFDM_num
    L = 200;
    h = sqrt(0.5)*(randn(1,L)+1i*randn(1,L));
    x = randi([0 1],1,nc);
    y = Tx(x, L, N0(i), h, fft_size);
    x_hat = Rx_MMSE( L, y,h,N0(i), a, fft_size);
    Pe_per_frame = Pe_per_frame + sum(x~=x_hat);
  end
  Pe_vect_MMSE(i) = Pe_per_frame/(nc*OFDM_num);

end
figure(4)
semilogy(-10:20,Pe_vect_MMSE,'LineWidth',2);
grid on
xlabel('SNR(dB)', 'Interpreter', 'latex');
ylabel('BER', 'Interpreter', 'latex');
title('MMSE','Interpreter', 'latex');
%% Compare
semilogy(-10:20,Pe_vect_ZF,'LineWidth',2);
hold on
semilogy(-10:20,Pe_vect_MMSE,'LineWidth',2);
grid on
xlabel('SNR(dB)', 'Interpreter', 'latex');
ylabel('BER', 'Interpreter', 'latex');
title('Comparing MMSE and Zero forcing','Interpreter', 'latex');
legend('ZF','MMSE')
%% Functions
function y = Tx(x, L, N0, h, fft_size)
    coding = [-1 1];
    x_encoded = coding(x+1);
    X = ifft(x_encoded,fft_size);
    X = [X(length(X)-L+1+1:length(X)) X];
    noise = sqrt(N0/2).*(randn(1,(length(X) + length(h)-1))+1i*randn(1,(length(X) + length(h)-1)));
    y = conv(X,h) + noise;

end
function x_hat = Rx_ZF( L, Y,h, fft_size)
    y = Y(L:fft_size+L-1);
    x_ = fft(y)./ fft([h zeros(1,length(y)-length(h))]);
    x_hat = real(x_)>0;
end

function x_hat = Rx_MMSE( L, Y,h,N0, a, fft_size)
    y = Y(L:fft_size+L-1);
    H = fft([h zeros(1,length(y)-length(h))]);
    W = conj(H)./(abs(H).^2 + N0/a^2);
    x_ = fft(y).*W;
    x_hat = real(x_)>0;
end