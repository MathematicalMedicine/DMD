function[Age,  standOTE]=predictAgeandOTEByOneFit(Age,   loaevent) 
% Input Age : time of event (uncensored) or last seen( censored)
%       Steroid :  1 Steroid Yes -1 Steroid No
%       censored :   1 for  censored,   0 for uncensored
%       loaevent :   0 for  censored,   1 for uncensored

 
censored =zeros(size(loaevent));
censored(loaevent==1) = 0;
censored(loaevent==0) = 1;

verboseNo=0;

curdist='wbl';

%disp(curdist)
pd = fitdist(Age,curdist,'Censoring',censored);

param= [pd.A pd.B];

sizeN = length(Age)  ;
OTE = zeros(size(Age));
% calculate MMR and OTE
for i=1:sizeN
 
        if censored(i)==1
            % t is the conditional survial time.
            t=condiSurvival_Median( param, Age(i));  % at ct
            Age(i) = t;
            [OTE(i)]= calOTE(curdist,param,t,loaevent(i), verboseNo);  
        else
            [OTE(i)]= calOTE(curdist,param,Age(i),loaevent(i), verboseNo);      
        end
   
end

 
standOTE = (OTE - mean(OTE))/std(OTE);
 

%sprintf('ps has mean and std of %4.3f (%4.3f)', mean(ps), std(ps))
end
 

function [cS]=condiSurvival_Median(wblAB, t)
   % calculate the conditional survival time t+...
   % Use the medain instead of mean
   
    St=1-wblcdf(t,wblAB(1),wblAB(2));  % S(t)
    % Find the t2 s.t. 1-wblcdf(t,wblAB(1),wblAB(2)) = St/2;
    cS = wblinv(1-St/2,wblAB(1),wblAB(2));
        
%     % check the inverse
%     [t St (1-wblcdf(cS,wblAB(1),wblAB(2))) cS] 
%     pause
end


function [  OTE, DR ]= calOTE(curdist,param,t,loaeventInd,verbose)
    % Do NOT return MR, DR and MMR to avoid the confusion with PHR MR 
    yy =1-cdf(curdist,t,param(1),param(2));
  
    CHF=-log(yy);
    
    MR = loaeventInd-CHF;
    DR = sign(MR) .* sqrt(-2*(MR + loaeventInd*log(loaeventInd- MR)));
    
    %MMR
    newDeltai= -log(0.5);
    MMR = newDeltai-CHF;
    
    %OTE
    OTEtemp = sign(MMR) .* sqrt(-2*(MMR + log(newDeltai- MMR)));
    %OTE
    OTE = sign(MMR) .* sqrt(-2*( newDeltai-newDeltai*log(newDeltai) - CHF+ newDeltai* log(CHF)));
    
    if verbose
        [CHF MR DR MMR OTE OTEtemp]
        pause
    end

end
