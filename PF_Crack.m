function [thetaHat,rul]=PF_Crack
clear global; global DegraUnit initDisPar TimeUnit time y thres ParamName thetaTrue signiLevel ns ny nt np
    %====== PROBLEM DEFINITION 1 (Required Variables) ===========
WorkName='PF_Crack';            % work results are saved by WorkName
DegraUnit='Crack size (m)';     % degradation unit
TimeUnit='cycles';              % time unit (Cycles, Weeks, etc.)
time=[0:100:3500]';             % time at both measurement and prediction.
dt= 20;                         % time interval for degradation propagation
y=[0.0100 0.0109 0.0101 0.0107 0.0110 0.0123 0.0099 0.0113...
0.0132 0.0138 0.0148 0.0156 0.0155 0.0141 0.0169 0.0168]'; %[nyx1]: measured data
thres= 0.05 ;                   % threshold (critical value)
ParamName=['m'; 'C'; 's'; 'a'];                        %[npx1]: parameters' name to be estimated
initDisPar=[4.0 0.2; -23 1.1; 0.001 0; 0.01 0];        %[npx2]: prob. parameters of initial/prior dist
thetaTrue=[3.8;log(1.5e-10); 0.001; 0.01];             %[npx1]: true values of parameters
signiLevel= 5;                                         % significance level for C.I. and P.I.
ns= 5e3 ;                                              % ns: number of samples(particles)
    %=========== PROGNOSIS using PF =================
ny=length(y); 
nt=ny;
np=size(ParamName,1);  % np: number of parameters
for j=1:np     %% Initial Distribution
    param(j,:,1)=normrnd(initDisPar(j,1),initDisPar(j,2),1,ns);
end
for k=2:length(time)   %% Update Process or Prognosis
    % step1. prediction (prior)
    paramPredi=param(:,:,k-1);  %the parameters at the previous step are copied to the parameters at the current step
    nk=(time(k)-time(k-1))/dt;
    for k0=1:nk 
        paramPredi(np,:)=MODEL(paramPredi,dt); %the degradation levels, paramPredi(np,:), need to be propagated from the previous step using MODEL function
    end                                        %this propagation must be repeated nk times with dt interval
    if k<=ny   % (Update Process)
        % step2. update (likelihood)
        mu=paramPredi(np,:); s=paramPredi(np-1,:);
        zeta=sqrt(log(1+(s./mu).^2)); eta=log(mu)-0.5*zeta.^2;
        likel=lognpdf(y(k),eta,zeta); 
        % step3. resampling based on the inverse CDF method
        cdf=cumsum(likel)./sum(likel);  % the cumulative CDF of likelihood samples 
        for i=1:ns
            u=rand;  %uniform distribution between zero and one 
            loca=find(cdf>=u,1); param(:,i,k)=paramPredi(:,loca); % we want to to find the approximate location of likelihood sample based the inverse CDF method 
        end                                                       % The inverse CDF method is repeated by ns times to generate new ns samples
    else  % (Prognosis)
        param(:,:,k)=paramPredi;
    end
end
thetaHat=param(1:np-1,:,ny);        %% Final Sampling Results
paramRearr=permute(param,[3 2 1]);  %% Degradation Prediction
zHat=paramRearr(:,:,np);
mu=zHat(ny:end,:); s=paramRearr(ny:end,:,np-1);
zeta=sqrt(log(1+(s./mu).^2)); eta=log(mu)-0.5*zeta.^2;
degraPredi=lognrnd(eta,zeta);
%=========== POST-PROCESSING =================
degraTrue=[];     %% True Degradation
if ~isempty(thetaTrue)
    k=1;
    degraTrue0(1)=thetaTrue(np); 
    degraTrue(1)=thetaTrue(np);
    for k0=2:max(time)/dt+1
        degraTrue0(k0,1)= MODEL([thetaTrue(1:np-1); degraTrue0(k0-1)],dt);
        loca=find((k0-1)*dt==time,1);
        if ~isempty(loca) 
            k=k+1;
            degraTrue(k)=degraTrue0(k0);
        end
    end
end
rul=POST(thetaHat,degraPredi,degraTrue);   %% RUL & Result Disp
Name=[WorkName ' at ' num2str(time(ny)) '.mat']; save(Name);
end
%====================%
function z1=MODEL(param,dt)
global ParamName np
for j=1:np
    eval([ParamName(j,:) '=param(j,:);']); 
end
    %===== PROBLEM DEFINITION 2 (model equation) ===============
dsig=75;
z1=exp(C).*(dsig.*sqrt(pi*a)).^m.*dt+a;
    %===========================================================
end
