function FBAsolution=perform_FBA(t1)

%function form for "perform_FBA_20120425.m"

%written by the Senger Research Group, 4/25/2012

%this program performs simple FBA that also minimizes total flux

FBAsolution=optimizeCbModel(t1,[],'one');   %performs pFBA
%FBAsolution=optimizeCbModel(t1);            %performs FBA

