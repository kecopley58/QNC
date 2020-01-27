% Exercise 1
n = 10; % number quanta available
p = 0.2; % probability of release
binopdf(0,n,p) % probabilityof releasing zero
binopdf(1,n,p)
binopdf(3,n,p)
binopdf(4,n,p)
binopdf(5,n,p)
binopdf(6,n,p)
binopdf(7,n,p)
binopdf(8,n,p)
binopdf(9,n,p)
binopdf(10,n,p)

% Exercise 2
n = 14; % number quanta available
p = 0.1; % probability of release
binopdf(8,n,p) % probability of k=8 if p=0.1
p = 0.1:0.1:1.0; % range of possible p
k = 8; % 8 released
binopdf(k,n,p) % probability of 8 released for range of p

% Exercise 3
% Most of the code for this exercise is from the answer set
k1 = 8; % previous k
k2 = 5; % 5 released second time
n = 14; % number quanta available
ps = (0:0.01:1.0)';           % Array of possible release probabilities -- compute
                              % at a resolution of 0.01
nps = length(ps);             % Get length of array of values of p
probs = binopdf( ....         % Get value of the binomial distribution for each
   repmat([k1 k2], nps, 1), ...  % combination of k, n, p. The output is a matrix
   repmat([n1 n2], nps, 1), ...  % with two rows: 1) n1, k1  2) n2, k2
   repmat(ps, 1, 2));            % columns are different values of p
% Since assuming the experiments are independent, will just multiply
% to get total likelihood; for log will just sum
subplot(2,1,1); cla reset; hold on;
ylabel('Likelihood');
likelihoodFcn = prod(probs,2);   % Compute the product
plot(ps, likelihoodFcn);         % Plot it
maxLikelihood = max(likelihoodFcn); % Get the maximum likelihood
plot(ps(likelihoodFcn==maxLikelihood), maxLikelihood, 'ko');
subplot(2,1,2); cla reset; hold on;
ylabel('Log-likelihood');
logLikelihoodFcn = sum(log(probs),2); % Compute the sum
plot(ps, logLikelihoodFcn);      % Plot it
maxLogLikelihood = max(logLikelihoodFcn); % Get the maximum likelihood
plot(ps(logLikelihoodFcn==maxLogLikelihood), maxLogLikelihood, 'ko');

%Looking at how estimate improves with increased sample size
n = 14;                       % number of available quanta
pRelease = 0.3;               % assumed probability of release
ks = 0:n;                     % possible values of k
ps = 0:0.01:1.0;              % possible release probabilities
pdat = zeros(length(ps), 2);  % pre-allocate matrix to hold likelihoods per p
TINY = 0.0001;                % to avoid multiplying/taking logs of really small numbers
for sampleSize = round(logspace(0,3,30))  % try different sample sizes
   
   % Simulate experiments -- get simulated counts for the given n,
   % pRelease, and number of experiments
   sCounts = binornd(n,pRelease,sampleSize,1);
   
   % Plot experiment, theoretical binomial pdf
   subplot(3,1,1); cla reset; hold on;
   title(sprintf('Sample size = %d', sampleSize))
   ylabel('Probability');
   xlabel('Release count');
   xlim([ks(1) ks(end)]);
   
   % Plot normalized histogram of simulated counts
   sCountHistogram = hist(sCounts, ks);
   bar(ks, sCountHistogram./sum(sCountHistogram));
   
   % Plot theoretical pdf
   plot(ks, binopdf(ks, n, pRelease), 'ro-');
   
   % compute (log) lik for each p
   pdat(:,1) = 1; % initialize so we can keep track of product of likelihoods
   pdat(:,2) = 0; % initialize so we can keep track of sum of log-likelihoods
   
   % Loop through each possible value of release probability
   for pp = 1:length(ps)
      
      % Compute the probabilities of obtaining the data, given the assumed
      % release probabilty
      probs = binopdf(ks, n, ps(pp));
      
      % Loop through each possible value of k we could have obtained in the
      % simulated experiment
      for cc = 1:length(ks)
         
         % Did we actually obtain that particular value of k?
         if sCountHistogram(cc) > 0
            
            % Avoid really small numbers
            if probs(cc) < TINY
               lprob = log(TINY);
               pprob = TINY;
            else
               lprob = log(probs(cc));
               pprob = probs(cc);
            end
            
            % Product of likelihoods
            pdat(pp,1) = pdat(pp,1).*pprob.^sCountHistogram(cc);
            
            % Sum of log likelihoods
            pdat(pp,2) = pdat(pp,2)+lprob.*sCountHistogram(cc);
         end
      end
   end
   
   % Use binofit to estimate p from the simulated data. This uses a trick
   % that assumes all of the measurements are independent and lumps them
   % together as if they were one big experiment.
   [phat, pci] = binofit(sum(sCounts),sampleSize*14);
   
   % Plot product of likelihoods
   %
   subplot(3,1,2); cla reset; hold on;
   ylabel('likelihood');

   % Plot the likelihood function (product of likelihoods)
   plot(ps, pdat(:,1));
   
   % Find the maximum
   maxp = max(pdat(:,1));
   
   % Show the actual pRelease value as a dashed line
   plot(pRelease.*[1 1], [0 maxp], 'r--');
   
   % Show the values obtained from binofit + CI
   plot(pci, maxp.*[1 1], 'm-', 'LineWidth', 2);
   plot(phat, maxp, 'm*');
   
   % Show the maximum value of our computed likelihood function
   plot(ps(pdat(:,1)==maxp), maxp, 'ko', 'MarkerSize', 12); 
   
   % plot sum of log-likelihoods
   %
   subplot(3,1,3); cla reset; hold on;
   ylabel('log likelihood');
   xlabel('Release probability');
   axis([0 1 log(TINY)*1000 0]);

   % Plot the likelihood function (sum of log-likelihoods)
   plot(ps, pdat(:,2));
   
   % Find the maximum
   maxp = max(pdat(:,2));
   
   % Show the actual pRelease value as a dashed line
   plot(pRelease.*[1 1], [0 min(pdat(:,2))], 'r--');
   
   % Show the values obtained from binofit + CI
   plot(pci, maxp.*[1 1], 'm-', 'LineWidth', 2);
   plot(phat, maxp, 'm*');

   % Show the maximum value of our computed likelihood function
   plot(ps(pdat(:,2)==maxp), maxp, 'g*'); 
   
   % Wait
   pause(0.1);
end

% Exercise 4
% Again, this is code from the answer set as I am still inexperiencing 
% with coding in matlab; I believe I mostly follow the code though
counts = [0	0 3 10 19 26 16 16 5 5 0 0 0 0 0]; % The experimental outcomes
n = length(counts)-1;   % Number of available quanta in each experiment
ks = 0:n;               % Possible values of k, as a row vector
nks = length(ks);       % Length of k
ps = (0:0.01:1.0)';     % Possible values of release probability, as a column vector
nps = length(ps);       % Length of p

% Compute the value of the binomial distribution for each possible value of 
%  k, n, p. Make a matrix in which:
%     - columns correspond to different values of p
%     - rows correspond to different values of k
probs = binopdf( ...
   repmat(ks, nps, 1), ...    % Repeat ks-row vector for each p
   n, ...                     % Always assuming the same 'n'
   repmat(ps, 1, nks));       % Repeat ps-column vector for each k

% Make a matrix of outcomes (in rows) that are repeated along the columns so we can
% use them to compute likelihoods for each possible value of release probability (p)
countsMatrix = repmat(counts, nps, 1);

% Compute likelihood function, which takes the product of all likelihoods
% associated with each measured outcome. 
likelihoodFcn = prod(probs.^countsMatrix,2);
pHat_fromLiklihood = ps(likelihoodFcn==max(likelihoodFcn))

% Compute log-likelihood function, which takes the sum of all log-likelihoods
% associated with each measured outcome. 
logLikelihoodFcn = sum(log(probs).*countsMatrix,2);
pHat_fromLogLikelihood = ps(logLikelihoodFcn==max(logLikelihoodFcn))

% use binofit
pHat = binofit(sum(counts.*ks),sum(counts)*n)

% Exercise 5
n = 14; % number quanta available
k = 7; % 7 released
pHat = binofit(k,n) % calculating max likelihoods
pNull = 0.3; % true release p
binopdf(k,n,pNull) % did T have an effect -- is p<0.05? 

% Bonus exercise
% The code here is from the answer set
% The data table
data = [ ...
   4.0 615 206 33 2 0 0; ...
   3.5 604 339 94 11 2 0; ...
   0.0 332 126 21 1 0 0; ...
   2.0 573 443 154 28 2 0; ...
   6.5 172 176 89 12 1 0; ...
   3.0 80 224 200 32 4 0];

xs = 0:5; % x-axis

% For each session
num_sessions = size(data, 1);
for ii = 1:num_sessions 
   
   % Compute relevant variables
   nx = data(ii,2:end); % the count data
   N  = sum(nx); % the total number of trials
   m  = sum(nx(2:end).*xs(2:end))/N;  % mean
   v  = sum((xs-m).^2.*nx)/N;  % variance
   p  = 1 - (v/m); % release probabilty
   n  = m/p; % available quanta per trial

   % Set up the plot
   subplot(1, num_sessions, ii); cla reset; hold on;
   bar(xs+0.5, nx./N, 'FaceColor', 'k');

   % Compute the binomial probabilities according to the equations at the
   % top of p. 762
   binomialCounts = zeros(size(xs));
   binomialCounts(1) = sum(nx).*(1-p).^n;   
   for jj = 2:length(binomialCounts)
       binomialCounts(jj) = binomialCounts(jj-1).*(m-p.*(jj-2))/((jj-1).*(1-p));
   end
   binomialCounts = round(binomialCounts);

   % Normalize for pdf and plot
   plot(xs+0.5, binomialCounts./sum(binomialCounts), 'ro-', 'MarkerFaceColor', 'r', 'LineWidth', 2);

   % If you want to compute Chi-2 goodness-of-fit,  k-1 degrees of freedom
   % A little bit of a cheat -- assume all bins contribute even when
   % binomialCounts=0 (because nx is always zero then, too)
   % pb = 1-chi2cdf(nansum((nx-binomialCounts).^2./binomialCounts), length(binomialCounts)-1);   

   % Get Possion pdf
   pps = poisspdf(xs, m);
   plot(xs+0.5, pps, 'bo-', 'MarkerFaceColor', 'b', 'LineWidth', 2);

   % If you want to compute Chi-2 goodness-of-fit, k-1 degrees of freedom
   % poissonCounts = round(pps.*N);
   % pp = 1-chi2cdf(nansum((nx-poissonCounts).^2./poissonCounts), length(poissonCounts)-1);

   % Show titles,labels, legend
   axis([0 5 0 1]);
   set(gca, 'FontSize', 12);
   h=title({sprintf('Temp=%.1f', data(ii,1)); sprintf('p=%.2f', p)});
   set(h, 'FontWeight', 'normal');
   if ii == 1
      xlabel('Number of Events')
      ylabel('Probability')
      legend('data', 'binomial', 'poisson');
   else
      set(gca, 'YTickLabel', '');
   end
end
        
