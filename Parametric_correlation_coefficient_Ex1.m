% Exercise 1
winglength = [10.4	10.8	11.1	10.2	10.3	10.2	10.7	10.5	10.8	11.2	10.6	11.4]; 
% Values of wing lengths; x variable; cm
taillength = [7.4	7.6	7.9	7.2	7.4	7.1	7.4	7.2	7.8	7.7	7.8	8.3];
% Values of tail lengths; y variable; cm
% 1. Plot
plot(winglength,taillength, 'ko') % To plot wing length and tail length
% Use of 'ko' was taken from Josh's code-- not using it results in
% connected lines
xlabel('Wing length (cm)')
ylabel('Tail length (cm)')
% By eye the variables do look positively correlated

% 2. Calculate r(x,y) and r(y,x)
% Using equations
% Much of this code is based off Josh's answers
n = 12; % Could use n = length(winglength) instead
samplemeanx = sum(winglength)/n;
samplemeany = sum(taillength)/n;
ssex = sum((winglength - samplemeanx).^2);
ssey = sum((taillength - samplemeany).^2);
scovxy = sum((winglength - samplemeanx).*(taillength - samplemeany)); 
% Above is covariance between x and y
rxy = scovxy/(sqrt(ssex)*sqrt(ssey))
ryx = scovxy/(sqrt(ssey)*sqrt(ssex))

% Using matlab function
corrcoef(winglength, taillength) % Returns matrices of values for covariance
% Between each variable with itself as well as with the other
% Both methods return the same r value

% 3.
% Standard error of rxy
serxy = sqrt((1-rxy^2)/(n-2)) % This is the standard error of rxy
% 95% confidence intervals
z = 0.5*log((1+rxy)/(1-rxy)); % Fisher's z transformation of r; use log function for ln
sd = sqrt(1/(n-3)); % Standard deviation of z
zs = z+[1 -1].*norminv(0.025).*sd; % From Josh's code; the z criterion
CI95 = (exp(2.*zs)-1)./(exp(2.*zs)+1); % From Josh's code
disp(sprintf('sem=%.4f, 95 pct CI = [%.4f %.4f]', serxy, CI95(1), CI95(2))) % From Josh's code

% 4. 
% Does rxy have p<0.5 if Hnull is r=0?
% This code is from Josh's answers
tval = rxy/serxy;           % Compute t statistic
prob = 1-tcdf(tval,n-2); % Compute p, using n-2 degrees of freedom
disp(sprintf('p=%.4f for H0: r=0', prob))
% Thus it is significant

% 5.
% Is Yale's value significantly different from the one we measured?
% His r = 0.75
zme = z; % Calculated the z score for the measured data in part 3
zh = 0.5*log((1+0.75)/(1-0.75)); % z score for Yale's r
zhvsme = (zme-zh)/(sqrt(1/(n-3)));
prob2 = 1-tcdf(zhvsme/2,inf); % From Josh's code
disp(sprintf('p=%.4f for H0: r=0.75', prob2)) % From Josh's code
% It is the same-- not significantly different

% 6. Calculate power and n needed
% This part is taken from Josh's code
v = n-2;
zm = 0.5*log((1+rxy)/(1-rxy)); % Same as original z calculated

tcrit = tinv(1-0.05/2,n-2); % t statistic
rcrit = sqrt(tcrit^2/(tcrit^2+(n-2)));
zr = 0.5*log((1+rcrit)/(1-rcrit));

Zb = (zm-zr)*sqrt(n-3);
power = normcdf(Zb);
disp(['The power of the test of H0:r=0 is ' num2str(power)])

%calculate the n needed to ensure that H0 (r=0) is rejected 99% of the time
%when |r|>= 0.5 at a 0.05 level of significance
Often = 0.01; %we want to reject H0 (1-often)% of the time [here, 99%]
rho = 0.5; % when r >=0.5 -- I believe this is the effect size decided upon?
AlphaValue = 0.05; % and the hypothesis is tested at AlphaValue of significance
Zb = tinv(1-Often,inf);
Za = tinv(1-AlphaValue/2,inf);
zeta = 0.5*log((1+rho)/(1-rho));
SampleSize = round(((Zb+Za)/zeta)^2+3);
disp(['To reject H0:r=0 99% of the time, when r=0.5 and alpha=0.05, we need n>=' num2str(SampleSize)])

