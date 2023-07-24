%% Part 1
r = randi([10,1000],10e5,1);
Pr = @(d) -60 - 40.*log10(d./10);
received_power = Pr(r);
Plot(received_power,'Simplified Path Loss Model')
%% Part 2
N0 = 10^((-175-30)/10);
N0W = N0*1e6;
Pn = 10*log10(N0W);
d = linspace(10,1000,1000);
P = Pr(d);
SNR = P - Pn;
plot(log10(d),SNR)
xlabel('log(d)')
ylabel('SNR(dB)')
grid on
%% Part 3
sigma = 5;
d = linspace(10,1000,1000);
indx = randi([1,1000],1e5,1);
r = d(indx);
Pr = @(d,sigma) -60 - 40.*log10(d./10) - sigma.*randn(1,length(d));
received_power = Pr(r,sigma);
figure
Plot(received_power, ' Shadowing and Path Loss')

SNR = zeros(length(d),1);
for i = 1:length(d)
    place = r == d(i);
    power = received_power(place) - Pn;
    SNR(i) = mean(power);
end
figure
plot(log10(d), SNR)
xlabel('log(d)')
ylabel('SNR(dB)')
grid on
%% Part 4
%-------------Plot Outage Probability----------------
d = linspace(10,1000,1000);
% SNR < 18 dB => Pr - Pn < 18 dB => Pr < 18 + Pn dB = 18 - 145 = -127 dB
% -60 - 40log10(d/d0) + X < -127 => X < -127 + 60 + 40log10(d/d0) = -67 + 40log10(d/d0) 
% answer = Q([-67 + 40log10(d/d0)]/sigma)
x = -67 + 40*log10(d./10);
Outage_Prob = qfunc(-x./sigma);
plot(log10(d), Outage_Prob)
xlabel('log(d)')
ylabel('Outage Probability')
grid on
%--------------Outage Probability-------------------
% SNR < 18 dB => Pr - Pn < 18 dB => Pr < 18 + Pn dB = 18 - 145 = -127 dB
% we should find the total number of the expected SNRs that we obtained in
% the last part (vector SNR) which are less that 18 dB
outage = SNR < 18;
probability = sum(outage)/ length(d);
%% Part 5
%---------------Simulation----------------
R =  10^(2.65607)*10;
S = pi*R^2;
%---------------Theory--------------------
Pr = @(d) -30 - 40.*log10(d./10);
e = exp(1);
D = 1000;
Pn = 10*log10(10^(-175/10)*10^6);
a = (18 + Pn - Pr(D))/sigma;
b = 40*log10(e)/sigma;
C = qfunc(a) + exp((2-2*a*b)/b^2)*qfunc((2-a*b)/b);
S_th = pi*D^2*C;
%% Functions
function [] = Plot(x,s)
    [h,stats] = cdfplot(x);
    legend('Empirical CDF')
    xlabel('Received power (dB)')
    ylabel('F(Received power)')
    title(s)
end