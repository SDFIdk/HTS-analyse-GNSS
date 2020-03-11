function beta = styrkefunktion(delta,z,sigma)
  %Calculates the statistical 2-tailed power test of the z-test.
  % Called in the function StyrkeFunktionPlot.m on line 18.
  %inputs:
  %     delta: linspace(-5,5,500);
  %     z    : stats.z_crit (Z-test, calculated in linreg on line 172).
  %     sigma: stats.std_known ( Standard error of slope (known), calculated in linreg on line 166)
  
%CV*SE + mu = mean needed to reject the null (CV = Critical value).
   
beta = 1 - normcdf((z.*sigma - delta)./sigma) + normcdf( (-z.*sigma - delta)./sigma); %Two-tailed power test

%Calculates the probability of rejecting the null hypothesis when the null hypothesis is false for a range of proposed sample means (delta, linspace -5,5).

end

%On statistical power of z-test:
% power is the probability of rejecting null hypothesis when it is false. Note that a type 2 error is the probability of not rejecting the null hyp. when it is false,
% so power is simply 1 - the prob. of a type 2 error (1 - Beta) so power depends on mu, sigma, n and alpha.

#When a researcher designs a study to test a hypothesis, he/she should compute the power of the test (i.e., the likelihood of avoiding a Type II error).

#To compute the power of a hypothesis test, use the following three-step procedure.

#1) Define the region of acceptance for a hypothesis test.

#2) Specify the critical parameter value. The critical parameter value is an alternative to the value specified in the null hypothesis. 
#The difference between the critical parameter value and the value from the null hypothesis is called the effect size. That is, the effect size is equal to 
#the critical parameter value minus the value from the null hypothesis.

#3) Compute power. Assume that the true population parameter is equal to the critical parameter value, rather than the value specified in the null hypothesis.
# Based on that assumption, compute the probability that the sample estimate of the population parameter will fall outside the region of acceptance. 
# That probability is the power of the test.
  
  % The normcdf(X,MU,SIGMA) function computes the cumulative distribution function (CDF) at X of the normal distribution
  % with mean MU and standard deviation SIGMA. Default values are MU = 0 and SIGMA = 1.
  
  %The power of a statistical test is the probability of rejecting the null hypothesis when the null hypothesis is false.
  %Low power corresponds to little evidence, so high power is desirable (power lies between 0 and 1).