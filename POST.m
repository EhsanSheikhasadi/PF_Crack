%% [POST]: MATLAB Code for RUL Calculation and Results Plot
 
function rul=POST(thetaHat,degraPredi,degraTrue)
global DegraUnit ...
TimeUnit time y thres ParamName thetaTrue signiLevel ns ny nt
np=size(thetaHat,1);
perceValue=[50 signiLevel 100-signiLevel];
figure(1);                      %% Distribution of Parameters
for j=1:np; subplot(1,np,j);
 [frq,val]=hist(thetaHat(j,:),30);
 bar(val,frq/ns/(val(2)-val(1)));
 xlabel(ParamName(j,:));
end;
figure(2);                                %% Degradation Plot
degraPI=prctile(degraPredi',perceValue)';
f1(1,:)=plot(time(1:ny),y(:,1),'.k'); hold on;
f1(2,:)=plot(time(ny:end),degraPI(:,1),'--r');
f1(3:4,:)=plot(time(ny:end),degraPI(:,2:3),':r'); 
f2=plot([0 time(end)],[thres thres],'g');
legend([f1(1:3,:); f2],'Data','Median',...
              [num2str(100-2*signiLevel) '% PI'],'Threshold')
xlabel(TimeUnit); ylabel(DegraUnit);
minY = min(y(:,1));  % Minimum y value from your data
maxY =  thres;  % Maximum y value from either your degradation prediction intervals or the threshold, whichever is larger
ylim([minY * 0.9, maxY * 1.5]);  % Set y-axis limits with a little padding
i0=0;                                       %% RUL Prediction
if y(nt(1))-y(1)<0; coeff=-1; else coeff=1;end;
for i=1:ns;
 loca=find(degraPredi(:,i)*coeff>=thres*coeff,1);
 if isempty(loca); i0=i0+1;
  disp([num2str(i) 'th not reaching thres']);
 elseif loca==1; rul(i-i0)=0;
 else rul(i-i0)=...
  interp1([degraPredi(loca,i) degraPredi(loca-1,i)], ...
           [time(ny-1+loca) time(ny-2+loca)],thres)-time(ny);
 end
end;
rulPrct=prctile(rul,perceValue);
figure(3);                             %% RUL Results Display
[frq,val]=hist(rul,30);
bar(val,frq/ns/(val(2)-val(1))); hold on;
xlabel(['RUL' ' (' TimeUnit ')']);
titleName=['at ' num2str(time(ny)) ' ' TimeUnit];
title(titleName)
fprintf(['\n Percentiles of RUL at %g ' TimeUnit], time(ny))
fprintf('\n  %gth: %g,  50th (median): %g,  %gth: %g \n', ...
perceValue(2),rulPrct(2),rulPrct(1),perceValue(3),rulPrct(3))
if ~isempty(degraTrue);                  %% True Results Plot
 figure(1); % parameters
 for j=1:np; subplot(1,np,j); hold on;
  plot(thetaTrue(j),0,'kp','markersize',18);
 end;
 figure(2); % degradation
 sl=0; if ~isempty(nt); sl=ny-nt(end); end
 f3=plot(time(sl+1:sl+length(degraTrue)),degraTrue,'k');
 legend([f1(1:3,:); f2; f3],'Data','Median', ...
       [num2str(100-2*signiLevel) '% PI'],'Threshold','True')
 figure(3); % RUL
 loca=find(degraTrue*coeff>=thres*coeff,1);
 rulTrue=interp1([degraTrue(loca) degraTrue(loca-1)], ...
             [time(sl+loca) time(sl+loca-1)],thres)-time(ny);
 plot(rulTrue,0,'kp','markersize',18);
end
end