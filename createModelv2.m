function model = createModelv2(rxnAbrList,rxnNameList,rxnList,revFlagList,...
    lowerBoundList,upperBoundList,subSystemList,grRuleList,geneNameList,...
    systNameList)
%createModelv2 Create a COBRA model from inputs or an empty model
%structure if no inputs are provided.
%
% model = createModelv2(rxnAbrList,rxnNameList,rxnList,revFlagList,...
%    lowerBoundList,upperBoundList,subSystemList,grRuleList,geneNameList,...
%    systNameList)
%
%INPUTS
% rxnAbrList            List of names of the new reactions
% rxnNameList           List of names of the new reactions
% rxnList               List of reactions: format: {'A -> B + 2 C'}
%                       If the compartment of a metabolite is not
%                       specified, it is assumed to be cytoplasmic, i.e. [c]
%
%OPTIONAL INPUTS
% revFlagList           List of reversibility flag (opt, default = 1)
% lowerBoundList        List of lower bound (Default = 0 or -vMax)
% upperBoundList        List of upper bound (Default = vMax)
% subSystemList         List of subsystem (Default = '')
% grRuleList            List of gene-reaction rule in boolean format (and/or allowed)
%                       (Default = '');
% geneNameList          List of gene names (used only for translation
%                       from common gene names to systematic gene names)
% systNameList          List of systematic names
%
%OUTPUT
% model                 COBRA model structure
%
% Ines Thiele 01/09

%create blank model
% the order of the fields and type of initiation of some of them were
% modified by Hadi Nazem-Bokaee on 3/13/2012
model = struct();

model.rxns=cell(''); model.rxnNames=cell(''); model.S=sparse(0,0);
model.grRules=cell(''); model.genes=cell(''); 
model.rules=cell(''); model.rxnGeneMat=sparse(0,0); 
model.proteins=cell(''); model.subSystems=cell(''); model.rev=zeros(0,1);
model.lb=zeros(0,1);model.ub=zeros(0,1); 
model.c=zeros(0,1);model.b=zeros(0,1);
model.confidenceScores=cell(''); model.rxnECNumbers=cell('');
model.rxnNotes=cell(''); model.rxnReferences=cell('');

model.mets=cell('');model.metNames=cell('');
model.metFormulasNeutral=cell(''); model.metFormulas=cell(''); 
model.metCharge=zeros(0,1); model.metCompartment=cell('');
model.metKEGGID=cell(''); model.metPubChemID=cell('');
model.metChEBIID=cell(''); model.metInChIString=cell(''); 
model.metSmiles=cell('');

if nargin < 1
    return;
end

nRxns = length(rxnNameList);
if nargin < 9
    geneNameList(1:nRxns,1) = {''};
    systNameList(1:nRxns,1) = {''};
end
if nargin < 8
    grRuleList(1:nRxns,1) = {''};
end
if nargin < 7
    subSystemList(1:nRxns,1) = {''};
end
if nargin < 5
    lowerBoundList = -1000*ones(nRxns,1);
end
if nargin < 6
    upperBoundList = 1000*ones(nRxns,1);
end
if nargin < 4
    revFlagList = ones(nRxns,1);
end
if isempty(revFlagList)
    revFlagList = zeros(nRxns,1);
    revFlagList(find(lowerBoundList)< 0) = 1;
end

for i = 1 : nRxns
    if i==nRxns
        pause(eps)
    end
    if ~isempty(grRuleList{i})
        if ~isempty(strfind(grRuleList{i},','))
          grRuleList{i}= (regexprep(grRuleList{i},',',' or ')); 
        end
        % the following 3 lines were added on 2/28/12 by Hadi Nazem-Bokaee
        if ~isempty(strfind(grRuleList{i},'|'))
          grRuleList{i}= (regexprep(grRuleList{i},'|',' or ')); 
        end
        if ~isempty(strfind(grRuleList{i},'&'))
           grRuleList{i} = (regexprep(grRuleList{i},'&',' and '));
        end
       if ~isempty(strfind(grRuleList{i},'+'))
          grRuleList{i}= (regexprep(grRuleList{i},'+',' and '));
       end
    end
    [metaboliteList,stoichCoeffList,revFlag] = rxnParser(rxnList{i});
    for q=1:length(metaboliteList)
        if length(metaboliteList{q})<=3 || ~strcmp(metaboliteList{q}(end-2),'[')
            %assuming the default compartment is cytoplasmic
            metaboliteList{q}=[metaboliteList{q},'[c]'];
        end
    end
    % the following 4 lines were added on 11/23/11 by Hadi Nazem-Bokaee
    % setting the lower bound to zero if the reaction is irreversible
    if ~(revFlag)
       lowerBoundList(i) = 0; 
    end
    model = addReactionv2(model,{rxnAbrList{i},rxnNameList{i}},metaboliteList,stoichCoeffList,...
        revFlag,lowerBoundList(i),upperBoundList(i),0,...
        subSystemList{i},grRuleList{i},geneNameList{i},systNameList{i},false);
end