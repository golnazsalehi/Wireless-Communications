%% Part 1
%----------Defining the Parameters-------------
a = 1e-6;
b = 10e-6;
pathDelays = unifrnd(a,b,1,15); 
Alpha = zeros(1,15);
theta = zeros(1,15);
for i = 1:15
    T = pathDelays(i)/1e-6;
    sigma = sqrt(1e-3*T^(-4)/2);
    Alpha(i) = raylrnd(sigma);
    theta(i) = unifrnd(0,pi/2);
end

fc = 3e9;
wave_length =  3*1e8/fc;
v = 30;
f_D = v/wave_length.*cos(theta);
Phi_D = 2*pi*f_D.*pathDelays;
%----------Impulse Response Calculation---------
f = 0:1e6;
f1 = repmat(f,15,1);
Tau = repmat(pathDelays,length(f),1)';
Phi = repmat(Phi_D,length(f),1)';
A = repmat(Alpha,length(f),1)';
H = sum( A.*exp(-1i*(Phi+2*pi*f1.*Tau)) );
subplot(1,2,1)
plot(f,abs(H),'LineWidth',2)  
xlabel('f(Hz)')
ylabel('|H(f)|')
grid on
subplot(1,2,2)
plot(f,angle(H),'r','LineWidth',2)
xlabel('f(Hz)')
ylabel('Phase(H(f))')
grid on
%% Part 2
a = 1e-6;
b = 10e-6;
pathDelays = unifrnd(a,b,15,1e5);
T = pathDelays./1e-6;
sigma = sqrt(1e-3*T.^(-4)/2);
Alpha = raylrnd(sigma);
Alpha = Alpha.^2;
power = sum(Alpha);
Plot(power,'Channel Power')

x = 0:0.1:10;
y = expcdf(x,2);
figure;
plot(x,y)
xlabel('Observation')
ylabel('Cumulative Probability')
grid on
%% Functions
function [] = Plot(x,s)
    [h,stats] = cdfplot(x);
    legend('Empirical CDF')
    xlabel('Received power (dB)')
    ylabel('F(Received power)')
    title(s)
end