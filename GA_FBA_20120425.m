%GA_FBA_20120425.m

%written by the Senger Research Group, 4/25/2012

%this program is called by the "Genetic_Algorithm.m" program and performs
%all FBA and calculations of unconstrained exchange fluxes.

%The input chromosome from the genetic algorithm is called "chrom1"

%The total unconstrained exchange fluxes are returned as the variable
%"result"

%this program also amends the stoichiometric matrix to update the biomass
%equation before each run.



%first, define what is in the chrom1 variable

mATP=chrom1(1,1);
glycogen=chrom1(1,2);
AA=chrom1(1,3);
lipid=chrom1(1,4);
DNA=chrom1(1,5);
RNA=chrom1(1,6);
peptidoglycan=chrom1(1,7);

%insert chromosome values into the stoichiometric matrix

%mATP
% r.S(200,150)=-mATP-0.171;  %ATP
r.S(151,150)=mATP;         %ADP
r.S(424,150)=mATP;         %H20
r.S(427,150)=-mATP-0.0005; %H+

%glycogen
r.S(406,150)=-glycogen;

%amino acids (keeping the fractions constant to the Reed model)
r.S(169,150)=-AA.*0.096;    %L-Ala
r.S(192,150)=-AA.*0.0053;   %L-Arg
r.S(195,150)=-AA.*0.0451;   %L-Asn
r.S(197,150)=-AA.*0.0451;   %L-Asp
r.S(244,150)=-AA.*0.0171;   %L-Cys
r.S(382,150)=-AA.*0.0492;   %L-Gln
r.S(388,150)=-AA.*0.0492;   %L-Glu
r.S(393,150)=-AA.*0.115;    %L-Gly
r.S(436,150)=-AA.*0.0177;   %L-His
r.S(455,150)=-AA.*0.0543;   %L-Iso
r.S(481,150)=-AA.*0.0842;   %L-Leu
r.S(490,150)=-AA.*0.0642;   %L-Lys
r.S(517,150)=-AA.*0.0287;   %L-Met
r.S(583,150)=-AA.*0.0346;   %L-Phe
r.S(611,150)=-AA.*0.0413;   %L-Pro
r.S(651,150)=-AA.*0.0403;   %L-Ser
r.S(692,150)=-AA.*0.0474;   %L-Thr
r.S(707,150)=-AA.*0.0106;   %L-Trp
r.S(714,150)=-AA.*0.0258;   %L-Tyr
r.S(750,150)=-AA.*0.0791;   %L-Val

%lipids (Reed fractions constant)
r.S(227,150)=-lipid.*0.011;  %cardiolipin
r.S(489,150)=-lipid.*0.923;  %lipopolysaccharide
r.S(580,150)=-lipid.*0.0549; %phosphatidylglycerol
r.S(614,150)=-lipid.*0.011;  %phosphatidylserine

%DNA (constant fractions)
r.S(254,150)=-DNA.*0.247;    %dATP
r.S(259,150)=-DNA.*0.253;    %dCTP
r.S(267,150)=-DNA.*0.253;    %dGTP
r.S(302,150)=-DNA.*0.247;    %dTTP

%RNA (constant fractions)
r.S(240,150)=-RNA.*0.198;    %CTP
r.S(417,150)=-RNA.*0.319;    %GTP
r.S(749,150)=-RNA.*0.214;    %UTP

%ATP
r.S(200,150)=-mATP-RNA.*0.287;

%Peptidoglycan
r.S(579,150)=-peptidoglycan;


%run FBA
perform_FBA_20120425

%calculate the unconstrained exchange reactions
unconstrained_exchange_calculation_20120425

%generate the result variable
result=EXuncon;