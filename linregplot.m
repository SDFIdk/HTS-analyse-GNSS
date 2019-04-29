function fig = linregplot(t, x, alpha, Beta0, sigma0, show_gpsweek, plot_middle_epoch, y_limit, varargin)
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

decimalsInPlots = "%.2f";

[px,sx] = linreg(t,x,alpha,Beta0,sigma0);

t2 = linspace(min(t),max(t),100);

fig = figure(1,"visible","off");
set(gcf, 'Position', get(0, 'Screensize'));

ax1 = axes('Position',[.15 .15 .7 .7]);
ax_text = axes('Position',[0 0 1 1],'Visible','off');

if (!isnan(Beta0))
  h = plot(ax1,t,x,'*k;GNSS observation;',...
     t2,polyval([px(2),px(1)],t2),['k-;Fit (R^2 = ' sprintf(decimalsInPlots,sx.Rsqr) ');'],...
     t2,polyval([Beta0,0],t2-mean(t))+mean(x),['r-;Beta_0 = ' sprintf(decimalsInPlots,Beta0)...
     ' [mm]/year (DTU Space absolute uplift model 2016);'],...
    
     sx.t2,sx.x2+sx.conf_mean_response,['-g;' num2str(100 - alpha.*100) '% Confidence interval for Estimated True Ellipsoidal Height;'],...
     sx.t2,sx.x2-sx.conf_mean_response,'-g;;',...     %sx.t2,sx.x2+sx.prediction,['--g;' num2str(100 - alpha.*100) '% Prediction interval for Predicted GNSS observation;'],... %sx.t2,sx.x2-sx.prediction,'--g;;',...
     sx.t2,sx.x2+sx.conf_mean_response_known,['-b;' num2str(100 - alpha.*100) '% Confidence interval for Estimated True Ellipsoidal Height (pooled);'],...
     sx.t2,sx.x2-sx.conf_mean_response_known,'-b;;');...
     %sx.t2,sx.x2+sx.prediction_known,['--b;' num2str(100 - alpha.*100) '% Prediction interval for Predicted GNSS observation;'],...
     %sx.t2,sx.x2-sx.prediction_known,'--b;;');
else
  h = plot(ax1,t,x,'*k;GNSS observation;',...
     t2,polyval([px(2),px(1)],t2),['k-;Fit (R^2 = ' sprintf(decimalsInPlots,sx.Rsqr) ');'],...

     sx.t2,sx.x2+sx.conf_mean_response,['-g;' num2str(100 - alpha.*100) '% Confidence interval for Estimated True Ellipsoidal Height;'],...
     sx.t2,sx.x2-sx.conf_mean_response,'-g;;',... %sx.t2,sx.x2+sx.prediction,['--g;' num2str(100 - alpha.*100) '% Prediction interval for Predicted GNSS observation;'],... %sx.t2,sx.x2-sx.prediction,'--g;;',...
     sx.t2,sx.x2+sx.conf_mean_response_known,['-b;' num2str(100 - alpha.*100) '% Confidence interval for Estimated True Ellipsoidal Height (pooled);'],...
     sx.t2,sx.x2-sx.conf_mean_response_known,'-b;;');...
     %sx.t2,sx.x2+sx.prediction_known,['--b;' num2str(100 - alpha.*100) '% Prediction interval for Predicted GNSS observation;'],...
     %sx.t2,sx.x2-sx.prediction_known,'--b;;');
end
     
legend(fig,"location","southoutside");
legend(fig,"boxoff");
hold on;

% Middle Epoch
if plot_middle_epoch
  plot(ax1,mean(t),mean(x),['ob;Middle Epoch Height: ' sprintf(decimalsInPlots,mean(x)) ' [mm] (Date:' datestr(mean(t)*365.25,'dd-mm-YYYY') ' , GPSweek: ' sprintf('%i',gpsweek(mean(t)*365.25)) ';']);
end

%Set y_limit
if (y_limit == 0)
      ylim(mean(ylim)+[-.5*y_limit .5*y_limit].*1e2);
end

set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);

xlabel(ax1,{'Date [years (GPS week)]'});
ylabel(ax1,'Ellipsoidal Height IGb08 [mm]');

% Show GPSweek
if show_gpsweek
  xticks = get(ax1, 'xticklabel');
  for i = 1:length(xticks)
    if length(xticks{i}) == 4
      xticks2{i} = char([xticks{i} " (" sprintf('%i',gpsweek(datenum(xticks{i},'YYYY'))) ")"]);
    else
      dateNumber = datenum(xticks{i}(1:4),'YYYY') + str2num(xticks{i}(6))*36.525;
      xticks2{i} = char([xticks{i} " (" sprintf('%i',gpsweek(dateNumber)) ")"]);
    end
  end
  set(ax1, 'TickLabelInterpreter', 'tex');
  set(ax1, 'xticklabel', xticks2);
end

%Estimated std
textstring = ["Estimated standard deviation on data " sprintf(decimalsInPlots,sx.sigma_0) " [mm]\n"];
textstring = [textstring "Estimated standard deviation on data from all measurements (pooled) "...
              sprintf(decimalsInPlots,sigma0) " [mm]\n \n"];
textstring = [textstring "Slope = " sprintf(decimalsInPlots,px(2)) " [mm/year] \n \n"];
textstring = [textstring "Using estimated standard deviation from observations:\n"];
textstring = [textstring "    Standard deviation {\sigma} = " sprintf(decimalsInPlots,sx.std_unknown) " [mm/year] \n"];
textstring = [textstring "    Est. Conf. interval = [ " sprintf(decimalsInPlots,sx.confinterval_estimated(1))  " - " sprintf(decimalsInPlots,sx.confinterval_estimated(2)) " ]\n"];
textstring = [textstring "    |t| = " sprintf(decimalsInPlots,abs(sx.t_score))];
if abs(sx.t_score) > sx.t_crit
  textstring = [textstring " > "];
  textstring2 = "    Null hypothesis rejected";
else
  textstring = [textstring " < "];
  textstring2 = "    Null hypothesis accepted";
end
textstring = [textstring "t_{crit} = " sprintf(decimalsInPlots,sx.t_crit) "\n" textstring2...
 " at " int2str(round(alpha*100)) "\% significance level\n \n"];

%Pooled std
textstring = [textstring "Using pooled standard deviation on observations:\n"];
textstring = [textstring "    Standard deviation \sigma = " sprintf(decimalsInPlots,sx.std_known) " [mm/year] \n"];
textstring = [textstring "    Known Conf. interval = [ " sprintf(decimalsInPlots,sx.confinterval_known(1))  " - " sprintf(decimalsInPlots,sx.confinterval_known(2)) " ]\n"];
textstring = [textstring "    |z| = " sprintf(decimalsInPlots,abs(sx.z_score_known))];
if abs(sx.z_score_known) > sx.z_crit
  textstring = [textstring " > "];
  textstring2 = "    Null hypothesis rejected";
else
  textstring = [textstring " < "];
  textstring2 = "    Null hypothesis accepted";
end
textstring = [textstring "z_{crit} = " sprintf(decimalsInPlots,sx.z_crit) "\n" textstring2...
 " at " int2str(round(alpha*100)) "\% significance level\n \n"];
textstring = [textstring "    "];
axes(ax_text);
text(0.15,0.75,textstring,'interpreter','tex')      

%If title supplied
if length(varargin) > 0
  title(ax1,varargin{1})
else %else use generic title.
  title('Linear Regression')
end
  