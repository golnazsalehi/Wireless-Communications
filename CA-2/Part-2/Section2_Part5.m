
nc = 1e5; 
fft_size = 1e5;
OFDM_num = 1e3;
SNR = 10.^((-10:20)/10);
a = 1;
N0 = a^2./SNR;
Pe_vect = zeros(1,length(SNR));
Pmax = nc;
C = zeros(1,length(SNR));

for i = 1:length(SNR)
  avg_Capacity = 0;
  Pe_per_frame = 0;
  for k = 1:OFDM_num
    L = 200;
    h = sqrt(0.5)*(randn(1,L)+1i*randn(1,L));
    [lambda,capacity] = Waterfilling(h,N0(i),nc,Pmax);
    avg_Capacity = avg_Capacity + capacity;
    x = randi([0 1],1,nc);
    y = Tx(x, L, N0(i), h, fft_size,lambda,nc);
    x_hat = Rx( L, y, fft_size);
    Pe_per_frame = Pe_per_frame + sum(x~=x_hat);
  end
  C(i) = avg_Capacity/OFDM_num;
  Pe_vect(i) = Pe_per_frame/(nc*OFDM_num);

end
figure(2)
subplot(1,2,1)
semilogy(-10:20,Pe_vect,'LineWidth',2);
grid on
xlabel('Pmax/(N0 nc)', 'Interpreter', 'latex');
ylabel('BER', 'Interpreter', 'latex');
subplot(1,2,2)
semilogy(-10:20,C,'r','LineWidth',2);
grid on
xlabel('Pmax/(N0 nc)', 'Interpreter', 'latex');
ylabel('Capacity', 'Interpreter', 'latex');
sgtitle('BER and Capacity using Water-Filling', 'Interpreter', 'latex');

%%
function [lambda,capacity] = Waterfilling(h,N0,nc,Pmax)
    H = fft([h zeros(1,nc-length(h))]);
    f = @(l) -sum( log2(1+ max(1/l - (N0./(abs(H).^2)),0).*abs(H).^2/N0 ) );
    x0 = 0.2;
    n0 = N0;
    save( "work.mat",'H','n0','Pmax')
    lambda = fmincon(f,x0,[],[],[],[],[],[],'subject');
    capacity = -f(lambda);
end
function y = Tx(x, L, N0, h, fft_size,lambda,nc)
    H = fft([h zeros(1,nc-length(h))]);
    P = max( (1/lambda-N0./abs(H).^2) , 0);
    coding = [-1 1];
    x_encoded = coding(x+1);
    x_encoded = x_encoded.*sqrt(P).*exp(-1i*angle(H));
    X = ifft(x_encoded,fft_size);
    X = [X(length(X)-L+1+1:length(X)) X];
    noise = sqrt(N0/2).*(randn(1,(length(X) + length(h)-1))+1i*randn(1,(length(X) + length(h)-1)));
    y = conv(X,h) + noise;
end
function x_hat = Rx( L, Y,fft_size)
    y = Y(L:fft_size+L-1);
    %y = y./sqrt(fft_size);
    x_ = fft(y);
    x_hat = real(x_)>0;
end