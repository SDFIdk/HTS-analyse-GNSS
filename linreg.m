function [b, stats] = linreg(t, x, varargin);
  % Linear regression and student t test for GNNS data
  % Nov 2017, oldjo@sdfe.dk
  %
  % outputs:
  %   b:                      polynomial coefficients
  %                           x = b(1) + b(2)*t + epsilon
  %   stats:                  structure with statistics
  %                           See code below for information on the fields
  %
  % inputs:
  %   t                       vector with time [decimal years]
  %   x                       vector with heights [mm]
  %   varargin(1) alpha       confidence interval
  %   varargin(2) Beta_0      DTU slope to compare to
  %   varargin(3) sigma_0_all all data sigma
 
  
  % Default values
  alpha = 0.05; %95% confidence interval
  Beta_0 = 0;   %Test against null hypothesis
      
  % Variable argument input    
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
  
  % Store the constants in 'stats' structure
  stats.Beta_0 = Beta_0;
  stats.alpha = alpha;
  stats.gamma = 1-alpha; %never used.
  
  % If inputs are row vectors, turn them into colums vectors:
  if size(t,1) < size(t,2)
    t = t';
  end
    if size(x,1) < size(x,2)
    x = x';
  end
  
  % Linear regression: x = b(1) + b(2)*t + epsilon
  T = [ones(length(t),1) t];
  b = T\x;
   
  % Calcualte R^2 stat
  stats.Rsqr = 1 - sum((x - (b(1)+b(2)*t)).^2)/sum((x - mean(x)).^2);
  
  N = length(t); % Num of points
  stats.N = N;
  stats.df = N - 2; % Degrees of freedom (N - (one for slope, one for intercept)
  
   
  % Residuals
  resid = (x - (b(1)+b(2).*t));
  stats.resid = resid;
  
  % Standard deviation of residuals
  stats.resid_std = std(resid);
  
  % Redundant code to confirm what "std" does.
  % This is the "corrected sample standard deviation":
  stats.s_n = sqrt((1/(N-2))*sum(resid.^2));
  
  % Aslak's mail (Nov. 23) identical to s_n
  stats.sigma_0 = sqrt((resid'*eye(length(resid))*resid)./(N-2));
  
  stats.sigma_0_sqrd_top = resid'*eye(length(resid))*resid;
  stats.sigma_0_sqrd_bottom = N-2;
  
  stats.sigma_hat = sqrt((1/(N-1.5))*sum(resid.^2));
  
  % Vandermonde matrix
  V = horzcat(ones(length(t),1),t);
  stats.V = V;
  
  % Sum of Squares of Residuals (p.127)
  stats.SSR = sum(resid.^2);
  
  % Mean Squared Residual (aka. Mean Squered Error) (p.127)
  MSR = stats.SSR/stats.df;
  stats.MSR = MSR;
  
  % Estimated Standar Error of Slope, (p.127) Chapter 5: Regression models
  stats.sigma_B_hat = sqrt((MSR*inv(V'*V))(2,2));
  
  % Equivalent:  
  % https://en.wikipedia.org/wiki/Student's_t-test#Slope_of_a_regression_line
  %stats.SE_beta_hat = sqrt( (1/(N-2) * sum( (x - ( b(1) + b(2)*t )).^2) ))...
  %                    / sqrt(sum( ( t - mean(t) ).^2 ));
              
  %Critical value of Student's t distribution for two tails.
  stats.t_crit = tinv(1-(alpha)/2,N-2);
  
  % Calculate a t_score agains either the null hypothesis or a supplied Beta_0  
  stats.t_score = ( (b(2) - Beta_0)*sqrt(N-2) ) / ...
                  sqrt( stats.SSR / sum( (t - mean(t)).^2 ) ); %possible division by zero                 
  % This were equivalent to the above.
  %stats.t_score3 = (b(2) - Beta_0)/stats.sigma_B_hat;

  % Students t-test
  if abs(stats.t_score) < stats.t_crit
    stats.t_test = 1;
  else  
    stats.t_test = 0;
  end
  
  % Calculate confidence interval, tivn is the Student's t distribution
  % NOTE: Slight confusion if 'tinv' uses the correct fractile. Ask Aslak
  stats.confinterval_estimated = tinv([alpha/2 1-(alpha/2)],N-2)*stats.sigma_B_hat+b(2);
  
  % Generate data for plots
  t2 = linspace(min(t),max(t),100);
  stats.t2 = t2;
  stats.x2 = polyval([b(2),b(1)],t2); % Best fit    
  
  stats.StdMiddleEpochEst = sqrt((MSR)*[1,mean(t)]*inv(T'*T)*[1;mean(t)]);
  
  warning("off"); % To avoid near singular warnings
  for i = 1:length(t2)
    stats.conf_mean_response(i) = stats.t_crit*sqrt((MSR)*[1,t2(i)]*inv(T'*T)*[1;t2(i)]);
    stats.prediction(i) = stats.t_crit*sqrt((MSR)*(1 + [1,t2(i)]*inv(T'*T)*[1;t2(i)]));
  end
  warning("on");
  
  % Estimated Covarience matrix SIGMA for the regression coefficients
  stats.SIGMAest = stats.sigma_0^2 .* (V'*eye(length(x))*V)^(-1);
  
  %STD students t
  stats.std_unknown = stats.t_crit*stats.sigma_B_hat/norminv(1-alpha/2);

  % Only do the below if sigma_0_all is defined
  if length(varargin) > 2
    
    stats.StdMiddleEpochPooled = sqrt((sigma_0_all^2)*[1,mean(t)]*inv(T'*T)*[1;mean(t)]);
    
    % Estimated Covarience matrix SIGMA for the regression coefficients using
    % pooled std    
    stats.SIGMAknown = (V'*(sigma_0_all.^(-2).*eye(length(x)))*V)^(-1);
        
    % Standard error of slope
    stats.std_known = sqrt(stats.SIGMAknown(2,2));
    
    % Confidence interval
    stats.confinterval_known = norminv([alpha./2 1-alpha./2])*stats.std_known+b(2);
           
    % Z test
    stats.z_crit = norminv(1-(alpha)/2);
    stats.z_score_known = (b(2) - Beta_0)/stats.std_known;
    % Hypothesis test:
    if abs(stats.z_score_known) < stats.z_crit
      stats.z_test_known = 1; % accepted
    else
      stats.z_test_known = 0; % rejected
    end
    
    % Generate data for plots
    warning("off");
    for i = 1:length(t2)
       stats.conf_mean_response_known(i) = ...
               stats.z_crit*sqrt((sigma_0_all^2)*[1,t2(i)]*inv(T'*T)*[1;t2(i)]);
               
       stats.prediction_known(i) = ...
        stats.z_crit*sqrt((sigma_0_all.^2)*(1 + [1,t2(i)]*inv(T'*T)*[1;t2(i)]));
    end
    warning("on");
   
   
  end
end