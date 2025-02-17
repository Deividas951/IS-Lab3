x = (0.1:1/22:1);
d = (1 + 1.5*sin(3*pi*x/2)) + 2.5*sin(2.8*pi*x)/0.5;
figure(1), plot(x, d, 'r*'), grid on;
% hidden layer weigths
w1 = randn(1);
w2 = randn(1);
b = randn(1);
% r1 = 0.152; % Good hand picked r values
% r2 = 0.126;
% c1 = 0.190909; Good hand picked c values
% c2 = 0.872727;

% step
n = 0.2;
%% finding centers automatically
[pks,locs,w,p] = findpeaks(d);
pksCount = length(pks);
c = zeros(1,pksCount);
rc = zeros(1,pksCount); % Indexes of peaks
peaks = 1;
for indx = 1:length(d)
    if pksCount+1 == peaks
        break;
    elseif d(indx) == pks(peaks)
            c(peaks) = x(indx);
            rc(peaks) = indx;
            peaks = peaks + 1;
    end
end
%% Network's response
% initializing r1 and r2 values, by simply calculating distance between one
% of the peaks to the nearest point below that peak.
cycleCount = 0;
       nc = 1;
       r1 = x(rc(1)+nc)-c(1);
       r2 = c(2) - x(rc(2)-nc);
       nc = nc+1;
%% First 100 iteration using initilized r1 and r2 values
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
%% More iteration, but modyfing r1 and r2 values until error is below 0.4
for ind = 1:10;
if e > 0.4
% r1 and r2 are being recalculated using points that are further from the
% peak, until error below defined threshold will be met
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
end

%% further decrease of error would be possible by increasing/decreasing r1 and r2 values by a very small number, but improvements are very small/

%% Test

for indx = 1:length(x)
   f1_1 = exp(-(x(indx)-c(1))^2/(2*r1^2));
   f2_1 = exp(-(x(indx)-c(2))^2/(2*r2^2));
   y = f1_1*w1+f2_1*w2+b;
   yfin = y;
   Yt(indx) = yfin;
   
end

figure
ylim([0,5])
plot(x, d, 'r*',x, Yt), grid on;