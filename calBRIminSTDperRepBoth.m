function [IntLRs, Mod ,mlenull] = calBRIminSTDperRepBoth(Steroid, dGenos, Age, LB,UB)
% This is a QT  steroid genotye interaction
% This calculates Max(LR) and Int(LR), where LR is comparing 2 hyps
%     Null: Steroid Yes has the same 3 mus as Steroid No,
%     Alt : Steroid Yes has the different 3 mus from Steroid No,
%     genotype     11       12       22
%    Steroid No    u11      u12      u22
%    Steroid Yes   u11+a11  u12+a12  u22+a22
%     
%     That is, null assumes a11=a12=a13
%              alt assumes a11~=a12~=a13
%     
%     First we assume all sigma are fixed at 1.
%     Then we vary them as well.
%
%     All QT values must be normalized for the analyses.   

%  We decided to use this fixed additive mu and multiplicative sigma on 12:30 pm Dec 4th, 2023
% This is now using OTE as QT values
% The range for alpha is -1 to 1 added to mus in Steroid No group.
% The range of beta is from [1/2, 2].
% We drop any genotype with counts <5.

% That is, the main difference from calBRI_2Dec23 is the range depending on
% the genotype (and rep) or not 

% This is Method 3 in the document.
% This calculates both BRI with / without GT='11'
 

global nullLod11 nullLod12 nullLod22 nullLodSum GTcounts maxLod minBound

sizeN= length(Age); 
minBound = -2*sizeN;

mus = zeros(1,3);
y1sigmas = zeros(1,3);
y2sigmas = zeros(1,3);

meanAge = mean(Age);
stdAge = std(Age);
Age = (Age- mean(Age))/std(Age);    % We changed the normalization scheme Jan 30, 2024 
                                    % Now we use Age   is ages of LOA (for uncensored) or ages of predicted LOA (for
                                    % censored)
    LOANo=Age((Steroid==-1));
    dGenoNo = dGenos(Steroid==-1);
    LOA0No=LOANo(dGenoNo==0);
    LOA1No=LOANo(dGenoNo==1);
    LOA2No=LOANo(dGenoNo==2);
    LOAYes=Age((Steroid==1));
    dGenoYes = dGenos(Steroid==1);
    LOA0Yes=LOAYes(dGenoYes==0);
    LOA1Yes=LOAYes(dGenoYes==1);
    LOA2Yes=LOAYes(dGenoYes==2);

    GTcounts = [length(LOA0No) length(LOA1No) length(LOA2No) length(LOA0Yes) length(LOA1Yes) length(LOA2Yes) ];
    GTstats = [mean(LOA0No) std(LOA0No) mean(LOA1No)  std(LOA1No) mean(LOA2No)  std(LOA2No)];
 
     
    mus(1)= mean(LOA0No);
    mus(2)= mean(LOA1No);
    mus(3)= mean(LOA2No);
    musYes(1)= mean(LOA0Yes);
    musYes(2)= mean(LOA1Yes);
    musYes(3)= mean(LOA2Yes);
    y1sigmas(1)= std(LOA0No);
    y1sigmas(2)= std(LOA1No);
    y1sigmas(3)= std(LOA2No);
    y2sigmas(1)= std(LOA0Yes);
    y2sigmas(2)= std(LOA1Yes);
    y2sigmas(3)= std(LOA2Yes);
   
 
    % Now there is no dropping
    if y1sigmas(1) <0.5
        y1sigmas(1) =   0.5;
        %mus(1)=musYes(1);
    end
    if y1sigmas(2) <0.5
        y1sigmas(2) =0.5;
    end
    if y1sigmas(3) <0.5
        y1sigmas(3) =0.5;
    end

    if y2sigmas(1) <0.5
        y2sigmas(1) =   0.5; 
    end
    if y2sigmas(2) <0.5
        y2sigmas(2) =   0.5; 
    end
    if y2sigmas(3) <0.5
        y2sigmas(3) =   0.5; 
    end

    nullLodSum=minBound; %log10(realmin)-500;
    options = optimoptions('fmincon' ,Display='none');
    fun = @(x) calNullLod(x,  mus,y1sigmas, LOA0Yes,LOA1Yes,LOA2Yes);
    x0= [0.5 1];
    A = [];
    b = [];
    Aeq =[];
    beq =[];
    lb= [-3 -3];
    ub= [3 3];
    nonlcon =[];
    %sprintf('repId=%d mus = %f %f %f', repId, mus)
    [mlenull , nullLod] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub, nonlcon, options);
  
    nullLod = - nullLod;

    maxLod = log10(realmin);
    funLR = @(x) calEpsLR(x, musYes,y2sigmas,LOA0Yes,LOA1Yes,LOA2Yes );
    ND=6;
    NF=1;
 

    MXPTS = 100000;
    EA = 1e-5;
    [IntLR, AE, NV, FL ] = adapt( ND, NF, LB, UB, MXPTS, funLR, EA); % ER, KEY );
    vol = prod(UB-LB);
    IntLR =IntLR/vol;

    thmax=[0 0 0 0 0 0];
    thMod = -calEpsLod(thmax, musYes,y2sigmas,  LOA0Yes,LOA1Yes,LOA2Yes ); % theoretical MOD

    Mod = max(maxLod, thMod); %-oppMod;
    
    IntLRs = log10(IntLR);
  
betaV=2^mlenull(2);
alphaV=mlenull(1)*stdAge  ;

mlenull=[alphaV betaV];
end
 



function [Lod] = calNullLod(x, mus,sigmas, LOA0Yes,LOA1Yes,LOA2Yes)
 
global nullLod11 nullLod12 nullLod22 nullLodSum  minBound

mus = mus +x(1);
sigmas = sigmas * 2.^x(2);

pLod0 = sum(log10(normpdf(LOA0Yes,mus(1),sigmas(1)))) ;
if ~isnan(pLod0)        
    if ~isreal(pLod0)
        pLod0= real(pLod0);
    end
    Lod=pLod0;
else
    pLod0 = NaN;
    Lod=0;
end

pLod1 = sum(log10(normpdf(LOA1Yes,mus(2),sigmas(2))));
if ~isnan(pLod1)  
    if ~isreal(pLod1)
        pLod1= real(pLod1);
    end
    Lod= Lod+pLod1;
else
    pLod1= NaN;
end

pLod2 = sum(log10(normpdf(LOA2Yes,mus(3),sigmas(3)))) ;
if ~isnan(pLod2)  
    if ~isreal(pLod2)
        pLod2= real(pLod2);
    end
    Lod= Lod+pLod2; 
else
    pLod2= NaN; 
end


if Lod > nullLodSum
    nullLodSum = Lod;
    nullLod11=pLod0;
    nullLod12=pLod1;
    nullLod22=pLod2;
 
end
if ~isreal(Lod)
    Lod= real(Lod);
end
if Lod < minBound 
    Lod = minBound; 
end

% we neee m.l.e. so, add - to find the min.
Lod = -Lod; 
end

%This is the version I calculate the altLike and nullLike separately.
function [LR] = calEpsLR(x, mus,sigmas,  LOA0Yes,LOA1Yes,LOA2Yes )
global nullLod11 nullLod12 nullLod22    maxLod

mus = mus+x(1:3)';
sigmas = sigmas .* 2.^x(4:6)'; 

if ~isnan(nullLod11) 
    pLod0 = sum(log10(normpdf(LOA0Yes,mus(1),sigmas(1))));
    Lod=pLod0-nullLod11;
     
else
    Lod=0;
end
 
if ~isnan(nullLod12)
    pLod1=sum(log10(normpdf(LOA1Yes,mus(2),sigmas(2)))); 
    Lod = Lod +pLod1 -nullLod12;
end

if ~isnan(nullLod22)
    pLod2 = sum(log10(normpdf(LOA2Yes,mus(3),sigmas(3))))  ;
    Lod = Lod +pLod2 -nullLod22;
end

LR= 10^(Lod); 
if Lod > maxLod
    maxLod = Lod;
     
end



end

%This is the version I calculate the altLike and nullLike separately.
function [Lod] = calEpsLod(x,  mus,sigmas, LOA0Yes,LOA1Yes,LOA2Yes )
global nullLod11 nullLod12 nullLod22  

mus = mus+x(1:3);
sigmas = sigmas .* 2.^x(4:6);
 

if ~isnan(nullLod11) 
    pLod0 = sum(log10(normpdf(LOA0Yes,mus(1),sigmas(1))));
    Lod=pLod0-nullLod11;
    %sprintf('using GT = 11 ')
else
    Lod=0;
end
 
if ~isnan(nullLod12)
    pLod1=sum(log10(normpdf(LOA1Yes,mus(2),sigmas(2)))); 
    Lod = Lod +pLod1 -nullLod12;
end

if ~isnan(nullLod22)
    pLod2 = sum(log10(normpdf(LOA2Yes,mus(3),sigmas(3))))  ;
    Lod = Lod +pLod2 -nullLod22;
end

Lod= -Lod; 

end 

 
 

