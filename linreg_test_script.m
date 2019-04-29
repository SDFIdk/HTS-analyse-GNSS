%linreg_test_script

%Generate test data
%t = 1:100;
%x = 1*t+10*randn(size(t));

% MPG example
%t = [2.560 2.300 1.975 1.915 2.020 2.815];
%x = [27.5 27.2 34.1 35.1 31.8 22.0];

% Load vs. Extension
%t = [3.072 3.154 3.238 3.322 3.403 3.487 3.569 3.652 3.737 3.816 3.902];
%x = [181.063 185.313 191.375 196.609 201.406 207.594 212.984 217.641...
%     223.609 228.203 234.422]

% Crystal Growth Example
%t = 2:2:28;
%x = [0.08 1.12 4.43 4.98 4.92 7.18 5.57 8.40 8.81 10.82 11.16 10.12 13.12 15.04];


%alpha = 0.05;
%Beta0 =  .50;
%sigma0 = 0.03;

[px,sx] = linreg(t,x,alpha,Beta0,sigma0);
%[px,sx] = linreg(t,x,alpha,Beta0);

t2 = linspace(min(t),max(t),100);

fig = figure;

h = plot(t,x,'*k;Data;',...
     t2,polyval([px(2),px(1)],t2),'b-;Fit;',...
     t2,polyval([Beta0,0],t2-mean(t))+mean(x),['r-;Beta_0 = ' num2str(Beta0) ';'],...
     sx.t2,sx.x2+sx.conf_mean_response,'-g;95% Confidence band for Estimated Mean Response;',...
     sx.t2,sx.x2-sx.conf_mean_response,'-g;;',...
     sx.t2,sx.x2+sx.prediction,'--g;95% Prediction band for Predicted Responses;',...
     sx.t2,sx.x2-sx.prediction,'--g;;');
     
     legend(fig,"location","southoutside");
hold on;
          
pkg load optim

    [px2,sx2] = polyfit(t,x,1);
    [xf,dxf] = polyconf(px2,t2,sx2,'pi');
    h = plot(t2,xf+dxf,'b--;95% Confidence Interval, polyconf;',...
      t2,xf-dxf,'b--;;');%,t,x,'*k;Data;');

      
textstring = ["Slope = " num2str(px(2)) char(177) num2str(sx.t_crit*sx.sigma_B_hat) "\n"];
textstring = [textstring "|t| = " num2str(abs(sx.t_score))];
if abs(sx.t_score) > sx.t_crit
  textstring = [textstring " > "];
  textstring2 = "Null hypothesis rejected";
else
  textstring = [textstring " < "];
  textstring2 = "Null hypothesis accepted";
end
textstring = [textstring "t_{crit} = " num2str(sx.t_crit) "\n" textstring2 "\n"];
textstring = [textstring "Est. Conf. interval = " num2str(sx.confinterval_estimated) "\n"];
textstring = [textstring "Known Conf. interval = " num2str(sx.confinterval_known) "\n"];
text(mean(t),min(x),textstring,'interpreter','tex')      
