function writeCbToSBMLv2(model,fileName)
%writeCbToSBML write a COBRA model to SBML
%
% writeCbToSBML(model,fileName)
%
%INPUTS
% model         COBRA model structure
% fileName      Name of xml file output
%
        sbmlModel = convertCobraToSBMLv2(model);
        OutputSBML(sbmlModel, fileName);
