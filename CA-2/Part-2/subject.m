function [c,ceq] = subject(l)
 c = [];
 work  =load('work.mat');
 H = work.H;
 N0 = work.n0;
 Pmax = work.Pmax;
 ceq = sum( max(1/l - N0./abs(H).^2,0) )- Pmax;
end