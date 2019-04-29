function fig = linregplot(t, x, alpha, Beta0, sigma0, varargin)
%LINREGPLOT 
%
%output: figure 
%
%input:
%  t: vector of dates in decimal years
%  x: vector of heights in mm
%  alpha:
%  Beta0:
%  sigma0: 
%  varargin: Optional title string
[px,sx] = linreg(t,x,alpha,Beta0,sigma0);

t2 = linspace(min(t),max(t),100);

fig = figure;

h = plot(t,x,'*k;Data;',...
     t2,polyval([px(2),px(1)],t2),'b-;Fit;',...
     t2,polyval([Beta0,0],t2-mean(t))+mean(x),['r-;Beta_0 = ' num2str(Beta0) ';'],...
     sx.t2,sx.x2+sx.conf_mean_response,'-g;95% Confidence band for Estimated True Value;',...
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
hold off;

ax1 = gca;
xlabel(ax1,{'Date [years (GPSweek)]'});
ylabel(ax1,'UTM32 Height [mm]');

xticks = get(ax1, 'xticklabel');
for i = 1:length(xticks)
  if length(xticks{i}) == 4
    xticks2{i} = char([xticks{i} " (" num2str(gpsweek(datenum(xticks{i},'YYYY'))) ")"]);
  else
    dateNumber = datenum(xticks{i}(1:4),'YYYY') + str2num(xticks{i}(6))*36.525;
    xticks2{i} = char([xticks{i} " (" num2str(gpsweek(dateNumber)) ")"]);
  end
end
set(ax1, 'TickLabelInterpreter', 'latex');
set(ax1, 'xticklabel', xticks2);

      
textstring = ["Slope = " num2str(px(2)) char(177) num2str(sx.t_crit*sx.sigma_B_hat) " [mm/year] \n"];
textstring = [textstring "Standard deviation sigma_hat: " num2str(sx.t_crit*sx.sigma_B_hat/1.96) " [mm/year] Z\n"];
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
textstring = [textstring "alpha = " num2str(alpha) "\n"];
text(mean(t),min(x),textstring,'interpreter','tex')      

if length(varargin) > 0
  title(varargin{1})
else
  title('Linear Regression')
end
  