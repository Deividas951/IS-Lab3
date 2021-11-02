x = (0.1:1/22:1);
d = (1 + 1.5*sin(3*pi*x/2)) + 2.5*sin(2.8*pi*x)/0.5;
figure(1), plot(x, d, 'r*'), grid on;
% hidden layer weigths
w1 = randn(1);
w2 = randn(1);
b = randn(1);
r1 = 0.152;
r2 = 0.126;
c1 = 0.190909;   %0.192
c2 = 0.872727;  %0.877
% step
n = 0.2;
%% Network's response
cycleCount = 0;
for ind = 1:10000;
for indx = 1:length(x)
   f1_1 = exp(-(x(indx)-c1)^2/(2*r1^2));
   f2_1 = exp(-(x(indx)-c2)^2/(2*r2^2));
   y = f1_1*w1+f2_1*w2+b;
   e = d(indx) - y;
   % Weight update
   w1 = w1 + n*e*f1_1;
   w2 = w2 + n*e*f2_1;
   b = b+n*e;
   cycleCount = cycleCount+1;
end
end

%% Test

for indx = 1:length(x)
   f1_1 = exp(-(x(indx)-c1)^2/(2*r1^2));
   f2_1 = exp(-(x(indx)-c2)^2/(2*r2^2));
   y = f1_1*w1+f2_1*w2+b;
   yfin = y;
   Yt(indx) = yfin;
   
end

figure
ylim([0,5])
plot(x, d, 'r*',x, Yt);