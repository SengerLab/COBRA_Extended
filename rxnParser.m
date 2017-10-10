function [metaboliteList,stoichCoeffList,revFlag] = rxnParser(equation)
%rxnParser parses a reaction equation to a list of metabolites and a
%list of stoichiometric coefficients. It also returns
%reversibility/irreversibility of the reaction in form of the logical
%variable, revFlag.

%INPUT
% equation          Reaction equation
%                   Any format can be handled with this function. Some
%                   examples are shown here:
%                   
%, may contain symbols '+', '->', '<=>' in
%                   addition to stoichiometric coefficients and metabolite names
%                   examples: 
%                   '0.01 cdpdag-SC[m] + 0.01 pg-SC[m]  -> 0.01 clpn-SC[m] + cmp[m] + h[m]' (irreversible reaction)
%                   'cit[c] + icit[x]  <=> cit[x] + icit[c] ' (reversible reaction)
%                   If no stoichiometric coefficient is provided, it is assumed
%                   to be = 1
%
%OUTPUTS
% metaboliteList    Cell array with metabolite names
% stoichCoeffList   List of S coefficients
% revFlag           Indicates whether the reaction is reversible (true) or
%                   not (false)
%
% Example:
%
%  formula = '0.01 cdpdag-SC[m] + 0.01 pg-SC[m]  -> 0.01 clpn-SC[m] + cmp[m] + h[m]'
%
%  [metaboliteList,stoichCoeffList,revFlag] = parseRxnFormula(formula)
%
%  metaboliteList = 
%   'cdpdag-SC[m]'    'pg-SC[m]'    'clpn-SC[m]'    'cmp[m]'    'h[m]'
%  stoichCoeffList = 
%   -0.01 -0.01 0.01 1 1
%  revFlag =
%   false
%
% Markus Herrgard 6/1/07
%
% Richard Que 1/25/10 Modified to handle '-->' and '<==>' as arrows 
% as well as reactionsformatted as '[compartment] : A --> C'. 
%
% Hadi 11/25/11     Modified to handle '=>', and '<->' as arrows as
%                   well as stoichiometric coefficients inside paranthesis as 
%                   '[compartment] : (2) A + (3) B => (4) C'
%                   The code still cannot handle '<=' as arrow
% Hadi Nazem-Bokaee 2/18/2012
%  
 
stoichCoeffList = [];
metaboliteList = {};
revFlag = true;
compartment = '';
reactants = []; products = [];
rxnSymbols = {' -> ',' <-> ',' <- ' ; ' --> ',' <--> ',' <-- ' ; ' => ',' <=> ',' <= ' ; ' ==> ',' <==> ',' <== '};
if strcmp(equation(1,1),'[')
    [endStr, matchStr] = regexp(equation, '\[\w+\]', 'end', 'match');
    compartment = char(matchStr);
    equation = equation(1,endStr+3:end);
end
for i = 1:size(rxnSymbols,1)
    for j = 1:size(rxnSymbols,2)
        symLoc{i,j} = strfind(equation,rxnSymbols{i,j});
    end
end
[r,c] = find(~cellfun('isempty',symLoc));
if c == 1
    revFlag = false;
    reactants = strtrim(equation(1,1:symLoc{r,c}));
    products = strtrim(equation(1,symLoc{r,c}+4:end));
elseif c == 2
%     revFlag = true;
    reactants = strtrim(equation(1,1:symLoc{r,c}));
    products = strtrim(equation(1,symLoc{r,c}+5:end));
elseif c == 3
    revFlag = false;
    reactants = strtrim(equation(1,symLoc{r,c}+4:end));
    products = strtrim(equation(1,1:symLoc{r,c}));
end
if ~isempty(reactants)
    rMets = regexp(reactants, '\s[+]', 'split');
else
    rMets = '';
end
if ~isempty(products)
    pMets = regexp(products, '\s[+]', 'split');
else
    pMets = '';
end
if ~isempty(rMets)
    for r = 1:length(rMets)
        temp_rMet = regexp(strtrim(rMets{1,r}), '\s', 'split');
        if length(temp_rMet) > 1
            if ~isnan(str2double(temp_rMet{1,1}))
                stoichCoeffList(end+1) = -1 * str2double(temp_rMet{1,1});
                k = 2;
            elseif strcmp(temp_rMet{1,1}(1,1),'(')
                prant = findstr(')', temp_rMet{1,1});
                stoichCoeffList(end+1) = -1 * str2double(temp_rMet{1,1}(2:prant-1));
                k = 2; prant = [];
            else
                stoichCoeffList(end+1) = -1;
                k = 1;
            end
            temp_name = [];
            for i = k:length(temp_rMet)
                temp_name = [temp_name ' ' char(temp_rMet{1,i})];
            end
        else
            stoichCoeffList(end+1) = -1;
            temp_name = char(temp_rMet);
        end
        if ~isempty(compartment)
            metaboliteList{end+1} = strtrim([temp_name compartment]);
        else % assuming the default compartment is cytosolic
            if strcmp(char(temp_name(1,end)),']')
                metaboliteList{end+1} = strtrim(char(temp_name));
            else
                metaboliteList{end+1} = strtrim([char(temp_name) '[c]']);
            end
        end
    end
end
if ~isempty(pMets)
    for p = 1:length(pMets)
        temp_pMet = regexp(strtrim(pMets{1,p}), '\s', 'split');
        if length(temp_pMet) > 1
            if ~isnan(str2double(temp_pMet{1,1}))
                stoichCoeffList(end+1) = str2double(temp_pMet{1,1});
                k = 2;
            elseif strcmp(temp_pMet{1,1}(1,1),'(')
                prant = findstr(')', temp_pMet{1,1});
                stoichCoeffList(end+1) = str2double(temp_pMet{1,1}(2:prant-1));
                k = 2; prant = [];
            else
                stoichCoeffList(end+1) = 1;
                k = 1;
            end
            temp_name = [];
            for i = k:length(temp_pMet)
                temp_name = [temp_name ' ' char(temp_pMet{1,i})];
            end
        else
            stoichCoeffList(end+1) = 1;
            temp_name = char(temp_pMet);
        end
        if ~isempty(compartment)
            metaboliteList{end+1} = strtrim([temp_name compartment]);
        else % assuming the default compartment is cytosolic
            if strcmp(char(temp_name(1,end)),']')
                metaboliteList{end+1} = strtrim(char(temp_name));
            else
                metaboliteList{end+1} = strtrim([char(temp_name) '[c]']);
            end
        end
    end
end