# Gene x Environment Interaction Analysis
This code calculates the PPI for Gene x Environment interaction.

Input: store the followins in csv file.
     Age: age of LOA (uncensored) or age of last seen (censored)
     Steroid : 1 for yes -1 for no steroid use
     dGenos : Genotypes 0 for homogeneous rare alleles
                        1 for heterogenous alleles
                        2 for homogeneous common alleles
     censored : 1 for censored and 0 for uncensored
 
Output: BRI, Mod, PPI

Usage: calPPI('repNM1.csv')


## Files included
calPPI.m : the main code to read the input file and calculates the PPI
calBRIminSTDperRepBoth.m 
calBRIminSTDperRepDropGTBoth.m
predictAgeandOTEByOneFit.m
repNM1.csv : a replicate in a simulated dataset

