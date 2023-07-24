%% Part 1: Narrowband channel
%% Q1
% section 1
SNR = 10.^((-20:20)/10);
P_e = zeros(41,1);
for i =1:41
    f = @(h) qfunc(sqrt(2*SNR(i))*h).*exp(-h.^2);
    P_e(i) = integral(f,-Inf,Inf)/sqrt(pi);
end
figure
semilogy((-20:20),P_e,'LineWidth',2)
grid on
xlabel('SNR(dB)')
ylabel('BER')
% section 2
P_e = qfunc(sqrt(2*SNR));
figure
semilogy((-20:20),P_e,'LineWidth',2)
grid on
xlabel('SNR(dB)')
ylabel('BER')
%% Q2
% section 1
SNR = 10.^((-20:20)/10);
P_e = 1./(2+SNR);
figure(1)
semilogy((-20:20),P_e,'LineWidth',2)
grid on
xlabel('SNR(dB)')
ylabel('BER')
% section 2
a = 1;
N = 2e5;
x = randi([0 1],[1 N]);
coding = [[0 a];[a 0]];
coded_x = coding(x+1,:);
SNR = 10.^((-20:20)/10);
N0_vector = a^2./SNR;
h = (randn(N,2)+1i*randn(N,2))*sqrt(0.5); % channel
P_e_2 = zeros(1,length(SNR));
for i = 1:length(SNR)
    w = sqrt(N0_vector(i)/2)*(randn(N,2)+1i*randn(N,2)); % noise
    y = coded_x.*h + w;
    Y = abs(y).^2;
    x_hat = ~((Y(:,1)-Y(:,2))>0);
    P_e_2(i) = sum(x == x_hat')/N;
end
figure(1)
hold on
semilogy((-20:20),P_e_2,'LineWidth',1)
legend('Theoretic', 'Simulation')
%% Q3
% section 1
SNR = 10.^((-20:20)/10);
P_e = zeros(length(SNR),1);
for i =1:length(SNR)
    f = @(hr,hi) qfunc(sqrt(2*SNR(i)*(hr.^2+hi.^2))).*exp(-(hr.^2+hi.^2));
    P_e(i) = integral2(f,-Inf,Inf,-Inf,Inf)/pi;
end
figure(2)
semilogy((-20:20),P_e,'LineWidth',2)
grid on
xlabel('SNR(dB)')
ylabel('BER')
% section 2
a = 2;
N = 2e5;
x = randi([0 1],[1 N]);
coding = [-1;1]*a;
coded_x = coding(x+1)';
SNR = 10.^((-20:20)/10);
h = (randn(1,N)+1i*randn(1,N))*sqrt(0.5); % channel
P_e_3 = zeros(1,length(SNR));
N0_vector = a^2./SNR;
for i = 1:length(SNR)
    w = sqrt(N0_vector(i)/2)*(randn(1,N)+1i*randn(1,N)); % noise
    y = coded_x.*h + w;
    Y = y.*conj(h)./abs(h).^2;
    x_hat = real(Y)<0;
    P_e_3(i) = sum(x == x_hat)/N;
end
figure(2)
hold on
semilogy((-20:20),P_e_3,'LineWidth',1)
grid on
xlabel('SNR(dB)')
ylabel('BER')
legend('Theoretic', 'Simulation')
% section 3
figure(3)
semilogy((-20:20),P_e_2,'LineWidth',2)
hold on
semilogy((-20:20),P_e_3,'LineWidth',2)
grid on
xlabel('SNR(dB)')
ylabel('BER')
legend('with coding, channel unknown', 'without coding, channel known')
%% Q4
% theoretic
SNR = 10.^((-20:20)/10);
P_e = zeros(41,1);
for i =1:41
    F = @(hr,hi) (1-qfunc(sqrt(SNR(i)*(hr.^2 + hi.^2.)))).^2 .* exp(-hr.^2 - hi.^2.);
    P_e(i) = 1-integral2(F,-Inf,Inf,-Inf,Inf)/pi;
end
figure(4)
semilogy((-20:20),P_e,'LineWidth',2)
grid on
xlabel('SNR(dB)')
ylabel('BER')
% simulation
a = 1;
QPSK = a*[1+1i,-1+1i,1-1i,-1-1i]/sqrt(2);
QPSK_map = [0 0; 1 0;0 1;1 1];
N = 1e6;
x = randi([0 1],[N 2]);
[~,Locb] =ismember(x,QPSK_map,'rows');
x_transmit = QPSK(Locb);
SNR = 10.^((-20:20)/10);
h = (randn(1,N)+1i*randn(1,N))*sqrt(0.5); % channel
P_e_4 = zeros(1,length(SNR));
N0_vector = a^2./SNR;
for i = 1:length(SNR)
    w = sqrt(N0_vector(i)/2)*(randn(1,N)+1i*randn(1,N)); % noise
    y = x_transmit.*h + w;
    Y = y.*conj(h)./abs(h).^2;
    x_hat = zeros(1,N);
    x_hat = x_hat + (real(Y)>0 & imag(Y)>0);
    x_hat = x_hat + (real(Y)<0 & imag(Y)>0)*2;
    x_hat = x_hat + (real(Y)<0 & imag(Y)<0)*4;
    x_hat = x_hat + (real(Y)>0 & imag(Y)<0)*3;
    P_e_4(i) = sum(Locb ~= x_hat')/N;
end
figure(4)
hold on
semilogy((-20:20),P_e_4,'LineWidth',1)
grid on
xlabel('SNR(dB)')
ylabel('Pe')
legend('Theoretic', 'Simulation')
% comapring part 2 and part 4
figure(5)
semilogy((-20:20),P_e_2,'LineWidth',1)
hold on
semilogy((-20:20),P_e_4,'LineWidth',1)
grid on
xlabel('SNR(dB)')
ylabel('Pe')
legend('BPSK with coding and CSI is known', 'QPSK, CSI is known')
%% Q5
% theoretic
SNR = 10.^((-10:10)/10);
P_e = zeros(21,5);
L = (1:5);
for j =1:5
    for i =1:21
        f = @(z) qfunc(sqrt(z*2*SNR(i))).*z.^(L(j)-1).*exp(-z)/gamma(L(j));
        P_e(i,j) = integral(f,0,Inf);
    end
end
figure(6)
for j =1:5
    semilogy((-10:10),P_e(:,j),'LineWidth',2)
    hold on
    grid on
end
xlabel('SNR(dB)')
ylabel('Pe')
% simulation
A = 2;
coding = [-1;1]*A;
N = 6e5;
L = (1:5);
SNR = 10.^((-10:10)/10);
N0_vector = A^2./SNR;
P_e_5 = zeros(21,5);
for j =1:5
    x = randi([0 1],[N 1]);
    x1 = coding(x+1)';
    x_ = repmat(x1,L(j),1);
    x_transmit = x_(:)';
    h = (randn(1,N*L(j))+1i*randn(1,N*L(j)))*sqrt(0.5); % channel
    for i = length(SNR):-1:1
        w = sqrt(N0_vector(i)/2)*(randn(1,N*L(j))+1i*randn(1,N*L(j))); % noise
        y = x_transmit.*h + w;
        a = abs(h)/sqrt(N0_vector(i));
        phase = exp(-1i*angle(h));
        z = y.*(a.*phase);
        if j ~= 1
            z = reshape(z,[L(j),N]);
            z = sum(z);
        end
        x_hat = real(z)>0;
        P_e_5(i,j) = sum(x ~= x_hat')/N;
    end
end
figure(6)
for j =1:5
    semilogy((-10:10),P_e_5(:,j),'LineWidth',2)
    hold on
    grid on
end
legend('L = 1, theoretic', 'L = 2, theoretic', 'L = 3, theoretic', 'L = 4, theoretic', 'L = 5, theoretic', ...
    'L = 1, simulation', 'L = 2, simulation', 'L = 3, simulation', 'L = 4, simulation', 'L = 5, simulation')
%% Q6
% simulation
a = 1;
N = 2e4;
x = randi([0 1],N,1);
coding = [-1;1]*a;
coded_x = coding(x+1)';
coded_x = reshape(coded_x,2,[]);
SNR = 10.^((-10:10)/10);
rep_x = repmat(coded_x,2,2);
rep_x = rep_x(1:2,:);
rep_x1 = rep_x(:,1:N/2);
rep_x2 = rep_x(:,N/2+1:end);
rep_x2 = flip(rep_x2,1);
rep_x2 = conj(rep_x2);
rep_x2 = [rep_x2(1,:)*-1; rep_x2(2,:)];
x_transmit = [rep_x1.',rep_x2.'].';
x_transmit = x_transmit(:);
h = (randn(2,N/2)+1i*randn(2,N/2))*sqrt(0.5); % channel
save_h = h;
h  = repmat(h,2,2);
h  = h(:,1:N/2);
h = h(:);
P_e_6 = zeros(1,length(SNR));
N0_vector = a^2./SNR;
for i = 1:length(SNR)
    received = x_transmit.*h;
    received = reshape(received.',2,[]);
    received = sum(received);
    w = sqrt(N0_vector(i)/2)*(randn(1,N)+1i*randn(1,N)); % noise
    y =  received+ w;
    y_ = reshape(y,2,[]);
    y_(2,:) = conj(y_(2,:));
    x_hat = [];
    for k = 1:N/2
        H_ = [conj(save_h(:,k)),[save_h(2,k);-save_h(1,k)]];
        Y = y_(:,k);
        S_hat = H_*Y;
        x_hat = [x_hat;S_hat>0];
    end
    P_e_6(i) = sum(x ~= x_hat)/N;
end
figure(7)
semilogy((-10:10),P_e_6,'LineWidth',2)
grid on
hold on
P_e = zeros(21,1);
for i =1:21
        f = @(z) qfunc(sqrt(z*2*SNR(i))).*z.^(2-1).*exp(-z)/gamma(2);
        P_e(i) = integral(f,0,Inf);
end
semilogy((-10:10),P_e,'LineWidth',2)
xlabel('SNR(dB)')
ylabel('Pe')
legend('2Tx-1Rx - Spatial Diversity','L = 2, Time diversity')