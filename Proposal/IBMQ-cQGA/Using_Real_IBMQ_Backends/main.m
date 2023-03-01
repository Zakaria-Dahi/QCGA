% -- Synchronous Cellular "Real" Quantum Genetic Algorithm --------
% Author: Dr. Zakaria Dahi
%Last time checked: 28-06-2020
% Initialisation: biased 1/sqrt(2)
% Interference: standard rotation
% crossover: two-point
% mutation: bit-flip
% selection: binary tournament
% grid: square
% Neighbourhood: Moore
% additional info: records the gbest (binary and probabilistic + variance of the population)
% experiments made using the "real mode": performed only on 15 cells instances.
% settings: except for the value of theta, we use the same parametring as cQLGA

clear all % clear in Matlab
clear classes % clear in Python


%% -=-=-= Uncomment to execute on the cluster -=-=-=-=-=
% read and set python interpreter: according to Andrés method
   % python = './zakvenv/bin/python';
   % if exist(python, 'file')
   %     pyversion(python);
   % end

%% -------- Initialisation of POI Libs -------------------------------------
% Add Java POI Libs to matlab javapath
javaaddpath('Jar/poi-3.8-20120326.jar');
javaaddpath('Jar/poi-ooxml-3.8-20120326.jar');
javaaddpath('Jar/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('Jar/xmlbeans-2.3.0.jar');
javaaddpath('Jar/dom4j-1.6.1.jar');
javaaddpath('Jar/stax-api-1.0.1.jar');
%% -------------------- Starting the execution of the program --------------

% Mode: 1 for real mode
% Mode: 0 for simulation mode

mode =1;

if mode == 1
    inst_tackle = 28; % Tackle only 15 cells (max 15 qubits in real mode IBMQ_16 melbourne)
	% later place a loop of execution exe=1:30
	% the index of the execution, will beused to load a different api token each time
    exe = 1; 
    % to import the API token once for all: avoid issues of "already loaded credentials"
    mod = py.importlib.import_module('import_api');
    py.importlib.reload(mod);
    response = py.dict(py.import_api.import_api_func());
else
    inst_tackle = 12; % tackle all the networks
end

for ind=inst_tackle:inst_tackle
%% -------------   Initialize the parameters of the experiments ------------
%% Save the best fitnesses reported by the state-of-the-art algorithm 
    % the 0 are for the instances higher after the classical 4x4, 6x6, 8x8 and 10x10 cells
global state_of_the_art 
       state_of_the_art = [98535  97156 95038 173701 182331  174519 307695 287149  264204  385927 357368  370868 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%%%% -------------- The number of execution ------------------------------
global execution 
       execution = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
%% ------ Number of individuals contained in the population --------------------
% Note:  
	% always remember the IBMQ restriction in real mode 8192 shots  = individuals (max)
global  indiv
        indiv = 400;
%% Save the dimension of the individual (it corresponds to the size of the network)
global Dimension
       Network_size = size(Instance(ind));
       Dimension = Network_size(1);
%% Number of itterations 
global iter
if mode == 1
    iter = ceil((175000 - indiv)/(2*indiv));
else
     %% No limits in simulation mode
    iter = ceil((175000 - indiv)/(2*indiv));
end   
%% Number of execution 
global Nexe
       Nexe = execution(ind);
%% Instances and neighbourhood of the network
global network 
       network = Instance(ind);
global neighbourhoud
       if ind == 1 || ind == 2 || ind == 3 || ind == 13 || ind == 16 ||  ind == 17
           neighbourhoud = neighbour(4);
       end
       if ind == 4 || ind == 5 || ind == 6 || ind == 14
           neighbourhoud = neighbour(6);
       end
       if ind == 7 || ind == 8 || ind == 9 || ind == 15
           neighbourhoud = neighbour(8);
       end
       if ind == 10 || ind == 11 || ind == 12
           neighbourhoud = neighbour(10);
       end
       if ind == 18 
          neighbourhoud =   neighbour(19);
       end
       if ind == 19
          neighbourhoud =   neighbour(63);
       end
       if ind == 20 
          neighbourhoud =   neighbour(99);
       end
       if ind == 21 
          neighbourhoud =   neighbour(144);
       end
       if ind == 22 
          neighbourhoud =   neighbour(196);
       end
       if ind == 23 
          neighbourhoud =   neighbour(256);
       end
       if ind == 24
          neighbourhoud =   neighbour(324);
       end
       if ind == 25
          neighbourhoud =   neighbour(400);
       end
       if ind == 26 
          neighbourhoud =   neighbour(900);
       end
       if ind == 27
          neighbourhoud =   neighbour(900);
       end
       if ind == 28 % this instance is used only for real mode
          neighbourhoud =   neighbour(15);
       end
%% ---- Saving the type and the dimension of the network used -------------
global instance 
      if ind == 1
          instance =  'Network_1_4x4';
      end
      if ind == 2
          instance =  'Network_2_4x4';       
      end
      if ind == 3
          instance =  'Network_3_4x4';         
      end
      if ind == 4
          instance =  'Network_1_6x6';         
      end
      if ind == 5
          instance =  'Network_2_6x6';         
      end
      if ind == 6
          instance =  'Network_3_6x6';         
      end
      if ind == 7
          instance =  'Network_1_8x8';         
      end
      if ind == 8
          instance =  'Network_2_8x8';
      end
      if ind == 9
          instance =  'Network_3_8x8';          
      end
      if ind == 10
          instance =  'Network_1_10x10';         
      end
      if ind == 11
          instance =  'Network_2_10x10';          
      end      
      if ind == 12
          instance =  'Network_3_10x10';
      end
      if ind == 13
          instance =  'Network_4_4x4';         
      end
      if ind == 14
          instance =  'Network_4_6x6';          
      end      
      if ind == 15
          instance =  'Network_4_8x8';
      end
      if ind == 16
          instance =  'Network_5_4x4';
      end
      if ind == 17
          instance =  'Network_6_4x4';
      end
      if ind == 18
          instance =  'Network_19_Cells';
      end
      if ind == 19
          instance =  'Network_7x9_Cells';
      end
      if ind == 20
          instance =  'Network_9x11_Cells';
      end
      if ind == 21
          instance =  'Network_12x12_Cells';
      end
      if ind == 22
          instance =  'Network_14x14_Cells';
      end
      if ind == 23
          instance =  'Network_16x16_Cells';
      end
      if ind == 24
          instance =  'Network_18x18_Cells';
      end
      if ind == 25
          instance =  'Network_20x20_Cells';
      end
      if ind == 26
          instance =  'Network_1_30x30_Cells';
      end
      if ind == 27
          instance =  'Network_2_30x30_Cells';
      end
      if ind == 28 % this instance is used only for real mode
          instance =  'Network_15x15_Cells';
      end
%% -------- Prameters  of the dynamic cellular Genetic Algorithm  ------
   pcr = 1;  % In all my works I found that Cr is better be equal to 1
   pm = 0.1;  % generally SOTA sets the Pm to be low.
   theta = 180/iter; % 180/218 = 0.8257, This sis done to give the theta a chance to reach 180 at the end of the iterations
   % pay attention: convert theta to radians because the rotation in IBMQ is in radians
   global grid_x % I set the grid to be square to show that this is a basic squeleton workable without any tuning
          grid_x = 20;
   global grid_y 
          grid_y = 20;
%% variable to save for each execution  ----------------------
     % Fitness of the best individual reach after the end of all the ietteration 
     global ALL_EXECUTION
            ALL_EXECUTION = [];
     % Save all the execution Time per run 
     global ALL_TIME
            ALL_TIME = [];
     % save all the best fitnesses through the generations
     global ALL_FITNESSES
            ALL_FITNESSES = [];
     % save all the itterations where the fitness value where extracted
     global ALL_ITTERATIONS
            ALL_ITTERATIONS = [];
     % Mean of all the fitness reached  
     global moy
            moy =0;
     % Standard deviation of the fitness reached 
     global ecart 
            ecart = 0;
     % save all the information relative to the experiment 
     global vect_one 
            vect_one = [];
     global vect_two
            vect_two = [];
     % save to best individual all along the NEXE execution 
     global FITNESS_GBEST
            FITNESS_GBEST = 1000000000000000000000000000000000000;
     % save the best individual all along the Nexe exveution 
     global GBEST_ALL
            GBEST_ALL = [];
     %%% To save the number of execution needed to obtain the best result obtained by the state-of-the-art algorithms 
     global iter_needed 
            iter_needed = zeros(1,Nexe);
     %%% To save the time needed to obtain the best result obtained by the state-of-the-art algorithms 
     global time_needed 
            time_needed = zeros(1,Nexe);
%% ------------- Start of the execution loop -------------------------------
for exe=1:Nexe
%% ----- Change the seed of the Mersenne Twister Random Generator ----------
rng shuffle
% ------------------------------------------------------------------------
ti = 0;
tic
%% --- To reset The number Of Fitness Evolution ----------
it = 0;
%% variable used during the calcul ---------------------------------------
     % matrix to record the best inidividual (binary)
     gbest_bin_mat = [];
     % matrix to record the best inidividual (probabilistic)
     gbest_prob_mat = [];
     % matrix to record the variance of the population 
     var_element = [];
     var_mat = [];     
     % Vector conaiting the fitness of all individuals of the population 
     global fit
            fit = [];
     % Vector containing the best individual so far 
     global gbest
            gbest = [];
     % The population 
     global population
            population = [];
     % To save the itteration and the corresponding fitness value of this itteration (fitness evauation)
     ALL_FITNESSES = [];
     ALL_ITTERATIONS = [];
%% ------------ Generate the intial population ---------------------
if mode == 1 % the real mode
    population_bin = initialisation_real(indiv,Dimension,exe);
else % the simulation mode
    population_bin = initialisation(indiv,Dimension);
end
% at the begining al the theta angles are set to 0
population = zeros(indiv,Dimension);
%% ------ Compute the initial fitnesses for the initial population -------
for w=1:indiv
    result  = RC_Function(population_bin(w,:),Dimension,network,neighbourhoud);
    fit = [fit , result];
end
%% -- calculate the best individual (it's index , it's fitness and ..) ---
[fitness_best,indbest] = min(fit);
gbest = population_bin(indbest,:);
ALL_FITNESSES = [ALL_FITNESSES fitness_best]; 
ALL_ITTERATIONS = [ALL_ITTERATIONS indiv];
% ----- record the best individual --------------
gbest_bin_mat = [gbest_bin_mat;gbest];
gbest_prob_mat = [gbest_prob_mat; population(indbest,:)];
% ----- record the variance of the population ---
for k=1:Dimension
    var_element = [var_element var(population(:,k))];
end
var_mat = [var_mat;var_element];

for it=1:iter % --------- itteration loop  --------------------------------
% ------  Initialise the parameters for info recording: variance ----------
var_element = [];
gbest_prob = [];
% --------------- Create the auxiliary population  ------------------------
futur_population = population;
futur_fit = fit;
% ----------- Stock the produced offspring when browsing ech node ---------
produced_cx = []; % stock the realvalued offspring (in terms of angles)
produced_cx_bin = []; %  stock the binary offspring (in terms of 0 and 1)
mat_chosen = []; % sttock the sequencing of the processes parents.
% -------------------------------------------------------------------------
for i=0:(grid_y - 1)
    for j=1:grid_x
        % ----------- Select the First Parentat position i,j --------------
        parent_1 = (i * grid_x) + j;
        chosen = parent_1;
        mat_chosen = [mat_chosen chosen];
        parent_1 = population(chosen,:);
        fit_parent_1 = fit(chosen);
        % ----------- To ensure that during the last execution we brose only 200 individuals -------
        % PS: during the last loop we are allowed to use 0.25 of 800
        % evaluations = 200. So, we browse only half this number as
        % individuals. Because, each individual iplies 2 FEs.
   
        if (it == iter) && (chosen > (0.25 * 2 * indiv / 2))
            break;
        end
        % ----------- Select the second parent -----------------------
        y = floor(chosen/grid_x);
        x = chosen - (y * grid_x);
        neighbourhood = [];
        if y == 0
            %% ------ First element of the grid ---------
            if x == 1 
                x1 = x + 1;
                y1 = y;
                x2 = x;
                y2 = y+1;
                x3 = x;
                y3 = grid_y - 1;
                x4 = grid_x;
                y4 = y;
                % Moore
                x5 = x + 1;
                y5 = y + 1;
                x6 = x + 1;
                y6 = grid_y - 1;
                x7 = grid_x;
                y7 = grid_y - 1;
                x8 = grid_x;
                y8 = y+1;
                candidate1 = (y1 * grid_x) + x1;
                candidate2 = (y2 * grid_x) + x2;
                candidate3 = (y3 * grid_x) + x3;
                candidate4 = (y4 * grid_x) + x4;
                candidate5 = (y5 * grid_x) + x5;
                candidate6 = (y6 * grid_x) + x6;
                candidate7 = (y7 * grid_x) + x7;
                candidate8 = (y8 * grid_x) + x8;
                neighbourhood = [candidate2 candidate1 candidate3 candidate4 candidate5 candidate6 candidate7 candidate8];
            else 
               %% ----------- First line of the grid -----------------------
                x1 = x + 1;
                y1 = y;
                x2 = x;
                y2 = y+1;
                x3 = x;
                y3 = grid_y - 1;
                x4 = x - 1;
                y4 = y;
                % Moore
                x5 = x + 1;
                y5 = y + 1;
                x6 = x - 1;
                y6 = y + 1;
                x7 = x + 1;
                y7 = grid_y - 1;
                x8 = x - 1;
                y8 = grid_y - 1;
                candidate1 = (y1 * grid_x) + x1;
                candidate2 = (y2 * grid_x) + x2;
                candidate3 = (y3 * grid_x) + x3;
                candidate4 = (y4 * grid_x) + x4;
                candidate5 = (y5 * grid_x) + x5;
                candidate6 = (y6 * grid_x) + x6;
                candidate7 = (y7 * grid_x) + x7;
                candidate8 = (y8 * grid_x) + x8;
                neighbourhood = [candidate4 candidate1 candidate2  candidate3 candidate5 candidate6 candidate7 candidate8];               
            end
        else
             if x == 0
                if y == 1
                   %% ------------ the last element of the first line -------------
                   x1 = grid_x - 1;
                   y1 = 0;
                   x2 = 0;
                   y2 = 2;
                   x3 = 1;
                   y3 = 0;
                   x4 = 0;
                   y4 = grid_y;
                   % Moore
                   x5 = grid_x - 1;
                   y5 = 1;
                   x6 = 1;
                   y6 = 1;
                   x7 = 1;
                   y7 = grid_y - 1;
                   x8 = grid_x - 1;
                   y8 = grid_y - 1;
                   candidate1 = (y1 * grid_x) + x1;
                   candidate2 = (y2 * grid_x) + x2;
                   candidate3 = (y3 * grid_x) + x3;
                   candidate4 = (y4 * grid_x) + x4;
                   candidate5 = (y5 * grid_x) + x5;
                   candidate6 = (y6 * grid_x) + x6;
                   candidate7 = (y7 * grid_x) + x7;
                   candidate8 = (y8 * grid_x) + x8;
                   neighbourhood = [candidate3 candidate1 candidate2 candidate4 candidate5 candidate6 candidate7 candidate8];
                else
                   if y == grid_y
                      %% ------------ the last element of the last line -------------
                      x1 = 1;
                      y1 = grid_y - 1;
                      x2 = 0;
                      y2 = grid_y - 1;
                      x3 = grid_x - 1;
                      y3 = grid_y - 1;
                      x4 = 0;
                      y4 = 1;
                      % Moore
                      x5 =  1;
                      y5 =  grid_y - 2;
                      x6 =  grid_x - 1;
                      y6 =  grid_y - 2;
                      x7 =  1;
                      y7 =  0;
                      x8 =  grid_x - 1;
                      y8 =  0;
                      candidate1 = (y1 * grid_x) + x1;
                      candidate2 = (y2 * grid_x) + x2;
                      candidate3 = (y3 * grid_x) + x3;
                      candidate4 = (y4 * grid_x) + x4;
                      candidate5 = (y5 * grid_x) + x5;
                      candidate6 = (y6 * grid_x) + x6;
                      candidate7 = (y7 * grid_x) + x7;
                      candidate8 = (y8 * grid_x) + x8;
                      neighbourhood = [candidate4 candidate1 candidate2 candidate3 candidate5 candidate6 candidate7 candidate8];
                   else
                      %% ------------ the last element of the rest line -------------
                      x1 = 0;
                      y1 = y - 1;
                      x2 = 1;
                      y2 = y - 1;
                      x3 = grid_x - 1;
                      y3 = y - 1;
                      x4 = 0;
                      y4 = y + 1;
                      % Moore
                      x5 =  1;
                      y5 =  y - 2;
                      x6 =  grid_x - 1;
                      y6 =  y - 2;
                      x7 =  1;
                      y7 =  y;
                      x8 =  grid_x - 1;
                      y8 =  y;
                      candidate1 = (y1 * grid_x) + x1;
                      candidate2 = (y2 * grid_x) + x2;
                      candidate3 = (y3 * grid_x) + x3;
                      candidate4 = (y4 * grid_x) + x4;
                      candidate5 = (y5 * grid_x) + x5;
                      candidate6 = (y6 * grid_x) + x6;
                      candidate7 = (y7 * grid_x) + x7;
                      candidate8 = (y8 * grid_x) + x8;
                      neighbourhood = [candidate1 candidate2 candidate3 candidate4 candidate5 candidate6 candidate7 candidate8];                      
                   end
                end
             else
                 if x == 1
                    if y == (grid_y - 1)
                       %% --------------- first element of the last line -----------
                        x1 = 1;
                        y1 = 0;
                        x2 = 1;
                        y2 = y - 1;
                        x3 = x + 1;
                        y3 = y;
                        x4 = 0;
                        y4 = y + 1;
                        % Moore
                        x5 =  2;
                        y5 =  0;
                        x6 =  grid_x;
                        y6 =  0;
                        x7 =  2;
                        y7 =  y - 1;
                        x8 =  grid_x;
                        y8 =  y - 1;
                        candidate1 = (y1 * grid_x) + x1;
                        candidate2 = (y2 * grid_x) + x2;
                        candidate3 = (y3 * grid_x) + x3;
                        candidate4 = (y4 * grid_x) + x4;
                        candidate5 = (y5 * grid_x) + x5;
                        candidate6 = (y6 * grid_x) + x6;
                        candidate7 = (y7 * grid_x) + x7;
                        candidate8 = (y8 * grid_x) + x8;
                        neighbourhood = [candidate1 candidate2 candidate3 candidate4 candidate5 candidate6 candidate7 candidate8];                 
                    else
                       %% --------------- first element of the remaining lines ----------- 
                        x1 = 1;
                        y1 = y - 1;
                        x2 = 1;
                        y2 = y + 1;
                        x3 = x + 1;
                        y3 = y;
                        x4 = 0;
                        y4 = y + 1;
                        % Moore
                        x5 =  2;
                        y5 =  y - 1;
                        x6 =  grid_x;
                        y6 =  y - 1;
                        x7 =  2;
                        y7 =  y + 1;
                        x8 =  grid_x;
                        y8 =  y + 1;
                        candidate1 = (y1 * grid_x) + x1;
                        candidate2 = (y2 * grid_x) + x2;
                        candidate3 = (y3 * grid_x) + x3;
                        candidate4 = (y4 * grid_x) + x4;
                        candidate5 = (y5 * grid_x) + x5;
                        candidate6 = (y6 * grid_x) + x6;
                        candidate7 = (y7 * grid_x) + x7;
                        candidate8 = (y8 * grid_x) + x8;
                        neighbourhood = [candidate1 candidate2 candidate3 candidate4 candidate5 candidate6 candidate7 candidate8];                          
                    end
                 else
                    if  y  == (grid_y - 1)
                       %% --------------- other element of the last line ----------- 
                        x1 = x;
                        y1 = y - 1;
                        x2 = x - 1;
                        y2 = y;
                        x3 = x + 1;
                        y3 = y;
                        x4 = x;
                        y4 = 0;
                        % Moore
                        x5 =  x + 1;
                        y5 =  y - 1;
                        x6 =  x - 1;
                        y6 =  y - 1;
                        x7 =  x + 1;
                        y7 =  0;
                        x8 =  x - 1;
                        y8 =  0;
                        candidate1 = (y1 * grid_x) + x1;
                        candidate2 = (y2 * grid_x) + x2;
                        candidate3 = (y3 * grid_x) + x3;
                        candidate4 = (y4 * grid_x) + x4;
                        candidate5 = (y5 * grid_x) + x5;
                        candidate6 = (y6 * grid_x) + x6;
                        candidate7 = (y7 * grid_x) + x7;
                        candidate8 = (y8 * grid_x) + x8;
                        neighbourhood = [candidate1 candidate2 candidate3 candidate4 candidate5 candidate6 candidate7 candidate8];                     
                    else
                       %% --------------- other element of the remaining lines ----------- 
                        x1 = x - 1;
                        y1 = y ;
                        x2 = x + 1;
                        y2 = y;
                        x3 = x;
                        y3 = y + 1;
                        x4 = x;
                        y4 = y - 1;
                        % Moore
                        x5 =  x + 1;
                        y5 =  y + 1;
                        x6 =  x - 1;
                        y6 =  y + 1;
                        x7 =  x + 1;
                        y7 =  y - 1;
                        x8 =  x - 1;
                        y8 =  y - 1;
                        candidate1 = (y1 * grid_x) + x1;
                        candidate2 = (y2 * grid_x) + x2;
                        candidate3 = (y3 * grid_x) + x3;
                        candidate4 = (y4 * grid_x) + x4;
                        candidate5 = (y5 * grid_x) + x5;
                        candidate6 = (y6 * grid_x) + x6;
                        candidate7 = (y7 * grid_x) + x7;
                        candidate8 = (y8 * grid_x) + x8;
                        neighbourhood = [candidate4 candidate1 candidate2 candidate3 candidate5 candidate6 candidate7 candidate8];                        
                    end
                 end
             end
        end
        % ------ Perform the interference: rotation gate --------------------
        % we are not really applying the rotation gate here but we are
        % modifying the theta that will be used later to apply the rotation around Z axis
        parent_processed = parent_1 + theta;
        parent_1 = parent_processed;
        % check if there is actually elements < 0 or > 180
        index_inf_zero = find(parent_1 < 0);
        index_sup_oneeighty = find(parent_1 > 180);
        if isempty(index_inf_zero) == 0
            parent_1(index_inf_zero) = 0;
        end
        if isempty(index_sup_oneeighty) == 0
            parent_1(index_sup_oneeighty) = 180;
        end
        % ------ Perform the binary tournament on the neighbourhoud --------- 
        fit_list = []; 
        for k=1:length(neighbourhood) 
            fit_list = [fit_list fit(neighbourhood(k))]; 
        end
        % ------ extract the best individual from the neighbourhood----------
        [chosen_value,ind_chosen] = min(fit_list); 
        indice = neighbourhood(ind_chosen); 
        parent_2 = population(indice,:); 
        %% ----------------- Perform a two-point crossover -------------------------
        if rand <= pcr 
            cross_point_one = randi(Dimension); 
            cross_point_two = randi(Dimension);  
            % -- condition to avoid the first point is after the second -----
            if cross_point_one > cross_point_two 
                a = cross_point_one; 
                cross_point_one = cross_point_two; 
                cross_point_two = a; 
            end
            c1 = [parent_1(1:cross_point_one) parent_2((cross_point_one+1):cross_point_two) parent_1((cross_point_two + 1):Dimension)];
            c2 = [parent_2(1:cross_point_one) parent_1((cross_point_one+1):cross_point_two) parent_2((cross_point_two + 1):Dimension)];
        else
            c1 = parent_1;
            c2 = parent_2;
        end
        %% --------------- Perform a bit-flip mutation ----------------------------
          % We do not really apply the mutation here
          % we compute the new theta to be used for the mutation later
          % to apply a bit flip we just need to apply 180 - theta
        c1_mutated = c1;
        c2_mutated = c2;
        c1_pm_con = rand(1,Dimension);
        c2_pm_con = rand(1,Dimension); 
        c1_mut_pos = find(c1_pm_con <= pm);
        c2_mut_pos = find(c2_pm_con <= pm);
        c1_mutated(c1_mut_pos) = 180 - c1(c1_mut_pos); 
        c2_mutated(c2_mut_pos) = 180 - c2(c2_mut_pos);
        %% ------------ Stock the produced offspring ------------------------
        produced_cx = [produced_cx;c1_mutated;c2_mutated];
    end
end
%% -------------- extract the size of the produced offspring -----------------
[lim_i,lim_j] = size(produced_cx);
% --- Send all the produced offspring for Interference and Measurement -------
	% previsouly we send each offspring by its own. This implies that we
	% create a quantum circuit for each offspring and send each circuit in its own
	% job. This implies that the fairshare queuing algorithm is applied to
	% queue our 800 jobs. This makes the exeperiments last VERY VERY long.
	% So, we decided to create all the circuits of all the offsprings once
	% Once for all we will send all the circuits in ONE job. This decrease
	% tremendously the time of experiments.
if mode == 1 % the real mode
   index_mat = 1;
   while index_mat <= lim_i
        % used to devide the offspring to blocks of 75 individual max
        % (max allowed circuit on IBMQ_melbourne backend is 75 quantum circuits)
        concat_produced_cx_bin = [];
        left = lim_i - index_mat;
        if left >= 74
           lb = index_mat;
           ub = index_mat + 74;
           index_mat = ub + 1;
           size_concat = 75;
        else
           lb = index_mat;
           ub = index_mat + left;
           index_mat = ub + 1;
           size_concat = left + 1;
        end
        % -=-=-=-=-=-=-=-=-=--=- Quantum Measurement -=-=-=-=-=-=-=-=-=-=-=
        concat_produced_cx_bin = InteferMeasure_real(produced_cx(lb:ub,:),Dimension,size_concat,exe);
        produced_cx_bin = [produced_cx_bin;concat_produced_cx_bin];
   end
end

iterated_chosen = 0; % variable to stock an index to browse the array of the "mat_chosen"
for i=1:2:(lim_i-1)
        iterated_chosen = iterated_chosen + 1; % iterate the index of the browser for "mat_chosen"
        if mode ==  1 % the real mode
            % ---- recover the binary-valued offspirng after measuement ----
            c1_bin = produced_cx_bin(i,:);
            c2_bin = produced_cx_bin(i+1,:);
        else % the simulation mode
            %  ------------ Quantum measurement ---------
				% we do not only apply the measurement here
				% we create the whole quantum circuit
				% it contains: interference and measurement          
            c1_bin = InteferMeasure(produced_cx(i,:),Dimension);
            c2_bin = InteferMeasure(produced_cx(i+1,:),Dimension);
        end
        %% ------------ Test if offs are feasible solutions -----------------
        if isequal(c1_bin,zeros(1,Dimension)) == 1
            c1_bin(randi(Dimension)) = 1;
        end
        if isequal(c2_bin,zeros(1,Dimension)) == 1
            c2_bin(randi(Dimension)) = 1;
        end
        %% ---------- evaluate the two new offspring ----------------------
        fit_off_1 = RC_Function(c1_bin,Dimension,network,neighbourhoud);
        fit_off_2 = RC_Function(c2_bin,Dimension,network,neighbourhoud);  
        %% ------------- extract the processed parent ---------------------
        chosen = mat_chosen(iterated_chosen);
        parent_1 = population(chosen,:);
        fit_parent_1 = fit(chosen);
        %% --- compare the two offspring and the parent being processed ---
        if fit_off_1 <= fit_parent_1
            elected_child = c1_mutated;
            elected_fit = fit_off_1;
        end
        if fit_off_2 <= fit_parent_1
            elected_child = c2_mutated;
            elected_fit = fit_off_2;
        end
        if (fit_off_1 > fit_parent_1) && (fit_off_2 > fit_parent_1)
            elected_child = parent_1;
            elected_fit = fit_parent_1;
        end        
        %% ---------- Update the population using binary tournament --------
        futur_population(chosen,:) = elected_child;
        futur_fit(chosen) = elected_fit;
        %% ---------- Check if one of the produced offspring is better than the Gbest ---
        if fit_off_1 <= fitness_best
            gbest = c1_bin;
            gbest_prob = c1_mutated;
            fitness_best = fit_off_1;
        end   
        if fit_off_2 <= fitness_best
            gbest = c2_bin;
            gbest_prob = c2_mutated;
            fitness_best = fit_off_2;
        end  
end
%% ------------ Update the population --------------------------------
population = futur_population;
fit = futur_fit;
%% Recording the number of itterations needed to obtain results as good or better than State-of-the-art Algorithms
if fitness_best <= state_of_the_art(ind)
   if iter_needed(exe) == 0 
       iter_needed(exe)  = it;
	   time_needed(exe)  = toc;
   end
end
% ----- record the best individual --------------
gbest_bin_mat = [gbest_bin_mat;gbest];
gbest_prob_mat = [gbest_prob_mat;gbest_prob];
% ----- record the variance of the population ---
for k=1:Dimension
    var_element = [var_element var(population(:,k))];
end
var_mat = [var_mat;var_element];
%% ------------------------------------------------------------------------
ALL_FITNESSES = [ALL_FITNESSES fitness_best]; 
if it == iter
    ALL_ITTERATIONS = [ALL_ITTERATIONS ((((it - 1)* 2 *indiv) + (0.25 * 2 * indiv)) + indiv)];
else
    ALL_ITTERATIONS = [ALL_ITTERATIONS ((it * 2 *indiv) + indiv)];
end
end %-------------- end of the loop ---------------------------------------
if FITNESS_GBEST > fitness_best
	FITNESS_GBEST = fitness_best;
	GBEST_ALL = gbest;
end
ALL_EXECUTION = [ALL_EXECUTION fitness_best];
oo = toc;
ti = ti + oo;
ALL_TIME = [ALL_TIME ti];
end
%% ----  Writing the result -----------------------------------------------
moy =  mean(ALL_EXECUTION);
ecart = std (ALL_EXECUTION);
datte = mat2cell(date,1);
timee =  num2str(sum(ALL_TIME));
timme =  mat2cell(timee,1);
%% ---- Writing on Excel File ---------------------------------------------------------------------------------------------
vect_one = [datte Nexe it indiv Dimension instance timme mean(ALL_TIME) std(ALL_TIME)  min(ALL_EXECUTION) max(ALL_EXECUTION) moy ecart];

headers = {'date', 'Nbr_Execution','Nbr_Fitness_Evaluations','Nbr_Individual','Dimension (cells)','Instance','Execution_Time','Mean','Std','Best','Worst','Mean','Std'};
name = strcat(instance,'.xls');
xlwrite(name,[headers;vect_one],1);
%% --- Recording the results obtained by the best individual  -------------------------------
datee =  num2str(date);
result = Gbest_Show(GBEST_ALL,Dimension,network,neighbourhoud);
vect_two = [datte result(1) result(2) result(3)];
headers = {'date','Fitness','Update Location Cost','Paging Cost'};
xlwrite(name,[headers;vect_two],2);

%% ---- Recording The Reporting Cells ID ------------------------------------------------------
% I commented this block because when using large-sized networks, writing the IDs of cells would 
% be too much for the excel sheet: there is a limited number of cell you can wirte in excel.
% vect_two = [cell2mat(result(4))'];
% headers = {'Reporting Cell'};
% xlwrite(name,vect_two',3);
%% ---- Recording the ID of Non Reporting Cells -----------------------------------------------
% vect_two = [cell2mat(result(5))];
% headers = {'Non Reporting Cell'};
% xlwrite(name,vect_two',4);

xlwrite(name,ALL_EXECUTION,6);
xlwrite(name,ALL_TIME,7);
xlwrite(name,transpose(ALL_FITNESSES),8);
xlwrite(name,transpose(ALL_ITTERATIONS),9);
xlwrite(name,transpose(fit),10);
xlwrite(name,transpose(iter_needed),11);
xlwrite(name,transpose(time_needed),12);

% write info on the population
xlwrite(name,gbest_bin_mat,13);
xlwrite(name,gbest_prob_mat,14);
xlwrite(name,var_mat,15);
end
%% I added this command because Daniel told me that if i don't add it it will ot escape and display the results of the run 
exit;
