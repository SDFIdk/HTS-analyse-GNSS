function [b, stats] = linreg(t, x, varargin);
  % linear regression and student t test
  % Nov 2017, oldjo@sdfe.dk
  
  %Default values
  alpha = 0.05; %95% confidence interval
  Beta_0 = 0;   %Test against null hypothesis
      
  if length(varargin) > 0
    alpha = varargin{1};
    if length(varargin) > 1
      Beta_0 = varargin{2};
      if isnan(Beta_0)
        Beta_0 = 0;
      end
      if length(varargin) > 2
        sigma_0_all = varargin{3};
      end
    end
  end
  
  stats.alpha = alpha;
  stats.gamma = 1-alpha;
  
  % If inputs are row vectors, transpose them:
  if size(t,1) < size(t,2)
    t = t';
  end
    if size(x,1) < size(x,2)
    x = x';
  end
  
  T = [ones(length(t),1) t];
  
  % Linear regression: x = b(1) + b(2)*t + epsilon
  b = T\x;
   
  % Calcualte R^2 stat
  stats.Rsqr = 1 - sum((x - (b(1)+b(2)*t)).^2)/sum((x - mean(x)).^2);
  
  N = length(t); %Num of points
  stats.N = N;
  stats.df = N - 2; %Degrees of freedom (N - (one for slope, one for intercept)
  t_years = t.\365.25;
  
  %Residuals in
  resid = (x - (b(1)+b(2).*t));
  stats.resid = resid;
  
  %Standard deviation of residuals
  stats.resid_std = std(resid);
  
  % Redundant code to confirm what "std" does.
  % This is the "corrected sample standard deviation":
  stats.s_n = sqrt((1/(N-2))*sum(resid.^2));
  
  % Aslak's mail (Nov. 23) identical to s_n
  stats.sigma_0 = sqrt((resid'*eye(length(resid))*resid)./(length(resid)-2));
  stats.sigma_0_sqrd_top = resid'*eye(length(resid))*resid;
  stats.sigma_0_sqrd_bottom = length(resid)-2;
  stats.sigma_hat = sqrt((1/(N-1.5))*sum(resid.^2));
  
  %Vandermonde matrix
  V = horzcat(ones(length(t),1),t);
  stats.V = V;
  
  % Sum of Squares of Residuals  
  stats.SSR = sum(resid.^2);
  MSR = stats.SSR/stats.df;
  stats.MSR = MSR;
  stats.sigma_b = (stats.s_n)/sqrt(sum( ( t_years - mean(t_years) ).^2 ));
  
  % Chapter 5: Regression models, p. 127 bottom 
  stats.sigma_B_hat = sqrt((MSR*inv(V'*V))(2,2));
     
  % https://en.wikipedia.org/wiki/Student's_t-test#Slope_of_a_regression_line
  stats.SE_beta_hat = sqrt( (1/(N-2) * sum( (x - ( b(1) + b(2)*t )).^2) ))...
              / sqrt(sum( ( t - mean(t) ).^2 ));
              
  %Critical value of Student's t distribution for two tails. (CORRECT)
  stats.t_crit = tinv(1-(alpha)/2,N-2);
  
  % Calculate a t_score agains either the null hypothesis or a supplied Beta_0  
  stats.t_score = ( (b(2) - Beta_0)*sqrt(N-2) ) / ...
                  sqrt( stats.SSR / sum( (t - mean(t)).^2 ) );                  
  stats.t_score2 = (b(2) - Beta_0)/stats.SE_beta_hat;
  stats.t_score3 = (b(2) - Beta_0)/stats.sigma_B_hat;
  
  % Calculate confidence interval
  stats.conf = stats.t_crit * stats.sigma_B_hat;
  stats.bmax = b(2) + stats.conf;
  stats.bmin = b(2) - stats.conf;
  stats.confinterval_estimated = tinv([alpha/2 1-(alpha/2)],N-2)*stats.sigma_B_hat+b(2);
  
  % Store the Beta_0 in stats
  stats.Beta_0 = Beta_0;
  
  t2 = linspace(min(t),max(t),100);
  stats.t2 = t2;
  stats.x2 = polyval([b(2),b(1)],t2);
    
  for i = 1:length(t2)
    stats.conf_mean_response(i) = stats.t_crit*sqrt((MSR.^2)*[1,t2(i)]*inv(T'*T)*[1;t2(i)]);
    stats.prediction(i) = stats.t_crit*sqrt((MSR.^2)*(1 + [1,t2(i)]*inv(T'*T)*[1;t2(i)]));
  end
  
  %Covarience matrix SIGMA
  stats.SIGMA3 = stats.sigma_0^2 .* (V'*eye(length(x))*V)^(-1);

  %Only do the below if sigma_0_all is defined
  if length(varargin) > 2
        
    stats.SIGMA4 = (V'*(sigma_0_all.^(-2).*eye(length(x)))*V)^(-1);
    stats.SIGMA5 = sigma_0_all^2 .* (V'*eye(length(x))*V)^(-1);
    
    stats.confinterval_known = norminv([0.025 0.975])*sqrt(stats.SIGMA5(2,2))+b(2);
    
     for i = 1:length(t2)
       stats.conf_mean_response(i) = stats.t_crit*sqrt((sigma_0_all^2)*[1,t2(i)]*inv(T'*T)*[1;t2(i)]);
       stats.prediction(i) = stats.t_crit*sqrt((sigma_0_all.^2)*(1 + [1,t2(i)]*inv(T'*T)*[1;t2(i)]));
    end
  end
end