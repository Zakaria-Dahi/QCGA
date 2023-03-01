function [pop_ini] = initialisation(indiv,dimension)

%   to prepare Matlab call to Python

% make sure to put any initialisation after "clear classes"

% not really necessary unless there is a change in the function itself
% clear classes

char_pop = []; % contains the generated population in forms of char array
pop_ini = zeros(indiv,dimension); % contains the generated population as normal integer array

% load the function: "initialisation_Py_Func" form  " initialisation.py" file
mod = py.importlib.import_module('initialisation');
py.importlib.reload(mod);
measured = py.dict(py.initialisation.initialisation_Py_Func(indiv,dimension));

% this is done to convert python dictionary into a Matlab structure
measured = py.json.dumps(measured);
measured = char(measured);
measured = jsondecode(measured);
classical_indiv = cell2mat(replace(fieldnames(measured),'x',''));
occurence = cell2mat(struct2cell(measured));

% create the population
% the flip is used to ensure the correspondance between the quantum and the classical registers
for i =1: length(occurence)
    char_pop = [char_pop; repmat(flip(classical_indiv(i,:)),occurence(i),1)];
end

% convert the char population to integer population
for i =1:indiv
    ones_index = find(all(ismember(char_pop(i,:),'1'),1));
    pop_ini(i,ones_index)= 1;
    % make sure that no individual equals [000000...]
    if isequal(pop_ini(i,:),zeros(1,dimension)) == 1
        pop_ini(i,randi(dimension))= 1;
    end
end
end