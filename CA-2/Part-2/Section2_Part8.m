%% MMSE
nc = 1e5; 
fft_size = 1e5;
OFDM_num = 1e3;
SNR = 10.^((-20:10)/10);
a = 1;
N0 = a^2./SNR;
Pe_vect_MMSE_clipped = zeros(1,length(SNR));

for i = 1:length(SNR)
  Pe_per_frame = 0;
  for k = 1:OFDM_num
    L = 200;
    h = sqrt(0.5)*(randn(1,L)+1i*randn(1,L));
    x = randi([0 1],1,nc);
    y = Tx_clipped(x, L, N0(i), h, fft_size);
    x_hat = Rx_MMSE( L, y,h,N0(i), a, fft_size);
    Pe_per_frame = Pe_per_frame + sum(x~=x_hat);
  end
  Pe_vect_MMSE_clipped(i) = Pe_per_frame/(nc*OFDM_num);

end
figure
semilogy(-20:10,Pe_vect_MMSE_clipped,'LineWidth',2);
grid on
hold on 
semilogy(-20:10,Pe_vect_MMSE,'LineWidth',2);
xlabel('SNR(dB)', 'Interpreter', 'latex');
ylabel('BER', 'Interpreter', 'latex');
title('Comparing Clipped MMSE and Unclipped MMSE','Interpreter', 'latex');
legend('Clipped MMSE','Unclipped MMSE')
%%
function y = Tx_clipped(x, L, N0, h, fft_size)
    coding = [-1 1];
    x_encoded = coding(x+1);
    X = ifft(x_encoded,fft_size);
    ABS = abs(X);
    Max = max(ABS);    
    for i=1:fft_size
        if abs(X(i))>0.8*Max
            X(i) = 0.8*Max*exp(1i*angle(X(i)));
        end
    end
    X = [X(length(X)-L+1+1:length(X)) X];
    noise = sqrt(N0/2).*(randn(1,(length(X) + length(h)-1))+1i*randn(1,(length(X) + length(h)-1)));
    y = conv(X,h) + noise;
end

function x_hat = Rx_MMSE( L, Y,h,N0, a, fft_size)
    y = Y(L:fft_size+L-1);
    H = fft([h zeros(1,length(y)-length(h))]);
    W = conj(H)./(abs(H).^2 + N0/a^2);
    x_ = fft(y).*W;
    x_hat = real(x_)>0;
end