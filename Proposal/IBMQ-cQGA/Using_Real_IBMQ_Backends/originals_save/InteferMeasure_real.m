function [offspring] = InteferMeasure_real(candidate,dimension,rows,exe)

%%   this function creates the quantum_circuit
%%   make sure to put any initialisation after "clear classes"
     %     not really necessary unless there is a change in the function itself
     %   clear classes

% convert the matlab matrix into a strcuture accepted by python: numpy.array
% we can only pass 1xN vectors from Matlab to Python
candidate= py.numpy.array(candidate(:)');
%% Start the fucntion
char_offspring = []; % contains the generated population in forms of char array
% char_offspring contains the generated population in forms of char array
offspring = zeros(rows,dimension); % contains the generated population as normal integer array


% load the function: "initialisation_Py_Func" form  " initialisation.py" file
mod = py.importlib.import_module('InteferMeasure_real');
py.importlib.reload(mod);
measured = py.dict(py.InteferMeasure_real.InteferMeasure_real_Py_Func(candidate,dimension,rows,exe));

% this is done to convert python dictionary into a Matlab structure
measured = py.json.dumps(measured);
measured = char(measured);
measured = jsondecode(measured);
classical_indiv = cell2mat(replace(fieldnames(measured),'x',''));
occurence = cell2mat(struct2cell(measured)); 

% create the offspring
% the flip is used to ensure the correspondance between the quantum and the classical bits
for i =1: length(occurence)
    char_offspring = [char_offspring; repmat(flip(classical_indiv(i,:)),occurence(i),1)];
end


% convert the char population to integer population
for i =1:rows
    ones_index = find(all(ismember(char_offspring(i,:),'1'),1));
    offspring(i,ones_index)= 1;
    % make sure that no individual equals [000000...]
    if isequal(offspring(i,:),zeros(1,dimension)) == 1
        offspring(i,randi(dimension))= 1;
    end
end
end



