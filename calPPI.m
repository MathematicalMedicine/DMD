function [] = calPPI(fileName)     
% This code calculates the PPI for Gene x Environment interaction
% Input: store the followins in csv file.
%     Age: age of LOA (uncensored) or age of last seen (censored)
%     Steroid : 1 for yes -1 for no steroid use
%     dGenos : Genotypes 0 for homogeneous rare alleles
%                        1 for heterogenous alleles
%                        2 for homogeneous common alleles
%     censored : 1 for censored and 0 for uncensored
% 
% Output: BRI, Mod, PPI
%
% Usage: calPPI('repNM1.csv')

 
    Tb = readtable(fileName);     
    %dGenos,Steroid,simLOA,simCT,simLOACT,LOAandBarLOA,censored
    Age = Tb.simLOACT;
    Steroid = Tb.Steroid;
    dGenos = Tb.dGenos;

    LOA = ones(size(Age));
    LOA(Tb.censored == 1) = 0;
    [~,  OTE]=predictAgeandOTEByOneFit(Age, LOA) ;

    counts = [sum((dGenos==0)&(Steroid==-1)) sum((dGenos==1)&(Steroid==-1)) sum((dGenos==2)&(Steroid==-1)) ...
        sum((dGenos==0)&(Steroid==1)) sum((dGenos==1)&(Steroid==1)) sum((dGenos==2)&(Steroid==1))];
    if (counts(1) <3 || counts(4) <3) &&( counts(1) <2 || counts(5) <3)
        BRI=0.0;
        Mod=0.0;
    else
        if (counts(1) <3 || counts(4) <3)
            LB= [-3 -3 -1 -1];
            UB= [3 3 log(3)/log(2) log(3)/log(2)];
            [BRI, Mod] = calBRIminSTDperRepDropGTBoth( Steroid,dGenos,  OTE, LB, UB);
            
        else
            LB= [-3 -3 -3 -1 -1 -1];
            UB= [3 3 3 log(3)/log(2) log(3)/log(2) log(3)/log(2)];
            [BRI, Mod ] = calBRIminSTDperRepBoth( Steroid,dGenos,  OTE, LB, UB);
            
        end

    end
    PPI = (10^BRI)*0.0004;
    PPI = PPI/(PPI + 0.9996); 
     
    sprintf('BRI= %f Mod = %f PPI=%f', BRI, Mod, PPI);
 

end


 




