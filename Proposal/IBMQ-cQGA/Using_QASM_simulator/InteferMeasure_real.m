function [offspring] = InteferMeasure_real(candidate,dimension,exe)

%%   this function creates the quantum_circuit
%%   make sure to put any initialisation after "clear classes"
     %     not really necessary unless there is a change in the function itself
     %   clear classes


%% Start the fucntion
% char_offspring contains the generated population in forms of char array
offspring = zeros(1,dimension); % contains the generated population as normal integer array


% load the function: "initialisation_Py_Func" form  " initialisation.py" file
mod = py.importlib.import_module('InteferMeasure_real');
py.importlib.reload(mod);
measured = py.dict(py.InteferMeasure_real.InteferMeasure_real_Py_Func(candidate,dimension,exe));

% this is done to convert python dictionary into a Matlab structure
measured = py.json.dumps(measured);
measured = char(measured);
measured = jsondecode(measured);
classical_indiv = cell2mat(replace(fieldnames(measured),'x',''));
%occurence = cell2mat(struct2cell(measured)); % no need occurence is always 1

% the flip is used to ensure the correspondance between the quantum and the classical bits
char_offspring = flip(classical_indiv);


% convert the char population to integer population
ones_index = find(all(ismember(char_offspring,'1'),1));
offspring(ones_index)= 1;
% make sure that no individual equals [000000...]
if isequal(offspring,zeros(1,dimension)) == 1
    offspring(randi(dimension))= 1;
end
end



