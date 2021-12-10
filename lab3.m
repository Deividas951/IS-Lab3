x = (0.1:1/22:1);
d = (1 + 1.5*sin(3*pi*x/2)) + 2.5*sin(2.8*pi*x)/0.5;
figure(1), plot(x, d, 'r*'), grid on;
% hidden layer weigths
w1 = randn(1);
w2 = randn(1);
b = randn(1);
% r1 = 0.152; % Good hand picked r values
% r2 = 0.126;
c1 = 0.190909;   %0.192
c2 = 0.872727;  %0.877
% step
n = 0.2;
errorStep = 0.000001
%% finding centers automatically
[pks,locs,w,p] = findpeaks(d)
pksCount = length(pks);
c = zeros(1,pksCount)
rc = zeros(1,pksCount) % Indexes of peaks
peaks = 1;
for indx = 1:length(d)
    if pksCount+1 == peaks
        break;
    elseif d(indx) == pks(peaks)
            c(peaks) = x(indx)
            rc(peaks) = indx;
            peaks = peaks + 1;
    end
end
%% Network's response
cycleCount = 0;
       nc = 1;
       r1 = x(rc(1)+nc)-c(1);
       r2 = c(2) - x(rc(2)-nc);
       nc = nc+1;
%% Multiple smaller iterations of traning using different r1 and r2 values if error is above 0.4
for ind = 1:100;
for indx = 1:length(x)
   f1_1 = exp(-(x(indx)-c(1))^2/(2*r1^2));
   f2_1 = exp(-(x(indx)-c(2))^2/(2*r2^2));
   y = f1_1*w1+f2_1*w2+b;
   e = d(indx) - y;
   % Weight update
   w1 = w1 + n*e*f1_1;
   w2 = w2 + n*e*f2_1;
   b = b+n*e;
   % center update
   cycleCount = cycleCount+1;
end
end
if e > 0.35
       r1 = x(rc(1)+nc)-c(1);
       r2 = c(2) - x(rc(2)-nc);
       nc = nc+1;
for ind = 1:100;
for indx = 1:length(x)
   f1_1 = exp(-(x(indx)-c(1))^2/(2*r1^2));
   f2_1 = exp(-(x(indx)-c(2))^2/(2*r2^2));
   y = f1_1*w1+f2_1*w2+b;
   e = d(indx) - y;
   % Weight update
   w1 = w1 + n*e*f1_1;
   w2 = w2 + n*e*f2_1;
   b = b+n*e;
   % center update
   cycleCount = cycleCount+1;
end
end
end
if e > 0.35
       r1 = x(rc(1)+nc)-c(1);
       r2 = c(2) - x(rc(2)-nc);
       nc = nc+1;
for ind = 1:100;
for indx = 1:length(x)
   f1_1 = exp(-(x(indx)-c(1))^2/(2*r1^2));
   f2_1 = exp(-(x(indx)-c(2))^2/(2*r2^2));
   y = f1_1*w1+f2_1*w2+b;
   e = d(indx) - y;
   % Weight update
   w1 = w1 + n*e*f1_1;
   w2 = w2 + n*e*f2_1;
   b = b+n*e;
   % center update
   cycleCount = cycleCount+1;
end
end
end
if e > 0.35
       r1 = x(rc(1)+nc)-c(1);
       r2 = c(2) - x(rc(2)-nc);
       nc = nc+1;
for ind = 1:100;
for indx = 1:length(x)
   f1_1 = exp(-(x(indx)-c(1))^2/(2*r1^2));
   f2_1 = exp(-(x(indx)-c(2))^2/(2*r2^2));
   y = f1_1*w1+f2_1*w2+b;
   e = d(indx) - y;
   % Weight update
   w1 = w1 + n*e*f1_1;
   w2 = w2 + n*e*f2_1;
   b = b+n*e;
   % center update
   cycleCount = cycleCount+1;
end
end
end


% very gradual increases of r1 and r2 to improve training result even
% further, improvement is very minimal.
for inda = 1:10000;
    if e > 0.355
       r1 = r1 - errorStep;
       r2 = r2 - errorStep;

for indx = 1:length(x)
   f1_1 = exp(-(x(indx)-c(1))^2/(2*r1^2));
   f2_1 = exp(-(x(indx)-c(2))^2/(2*r2^2));
   y = f1_1*w1+f2_1*w2+b;
   e = d(indx) - y;
   % Weight update
   w1 = w1 + n*e*f1_1;
   w2 = w2 + n*e*f2_1;
   b = b+n*e;
   % center update
   cycleCount = cycleCount+1;
end  
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
plot(x, d, 'r*',x, Yt), grid on;