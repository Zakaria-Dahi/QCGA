% -- Popupatio-based Incremental Learning Algorithm (PBIL) --------
% Source: Solving the reporting cells problem by using a parallel team of evolutionary algorithms
% Note: it applies two updates on the probability vector
% date of implementation: 10/01/2019
% Author: Dr Zakaria DAHI

clear all
%% --- Initialisation of POI Libs: Add Java POI Libs to matlab javapath ----
javaaddpath('Jar/poi-3.8-20120326.jar');
javaaddpath('Jar/poi-ooxml-3.8-20120326.jar');
javaaddpath('Jar/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('Jar/xmlbeans-2.3.0.jar');
javaaddpath('Jar/dom4j-1.6.1.jar');
javaaddpath('Jar/stax-api-1.0.1.jar');
%% -------------------- Starting the execution of the PBIL -----------------
for ind=1:12
    
%% -------------   Initialize the parameters of the experiments ------------

%% Save the best fitnesses reported by the state-of-the-art algorithm 
    % the 0 are for the instances higher after the classical 4x4, 6x6, 8x8 and 10x10 cells
       state_of_the_art = [98535  97156 95038 173701 182331  174519 307695 287149  264204  385927 357368  370868 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%% The number of execution 
       execution = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
%% Number of individuals contained in the population 
        indiv = 175;
%% Save the dimension of the individual (it corresponds to the network's size)
       Network_size = size(Instance(ind));
       Dimension = Network_size(1);
%% Number of itterations 
       iter = 1000;
%% Number of execution 
       Nexe = execution(ind);
%% Instances and neighbourhood of the network
       network = Instance(ind);
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
%% Saving the type and the dimension of the network used 
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
      
%% -------- Prameters  of the : PBIL  ------
    Lr = 0.1; % learning rate % value taken from the above-cited article
    new_lr = Lr/2; % see email clarification PBIL to understand.
    Pm = 0.075; % mutation probability % value taken from the above-cited article
    Ms = 0.075; % mutation shift % value taken from the above-cited article
%% variable to save for each execution  ----------------------
     % Fitness of the best individual reached at the end of each execution 
     ALL_EXECUTION = [];
     % Save all the execution's time of each execution
     ALL_TIME = [];
     % Mean of ALL_EXECUTION  
     moy =0;
     % Standard deviation of ALL_EXECUTION
     ecart = 0;
     % save all the information relative to the experiment 
     vect_one = [];
     vect_two = [];
     % save to best individual ever in that execution
     FITNESS_GBEST = 1000000000000000000000000000000000000;
     % save the best individual ever in that execution
     GBEST_ALL = [];
     %%% To save the number of execution needed to obtain the best result reached by SOTA 
     iter_needed = zeros(1,Nexe);
     %%% To save the time needed to obtain the best result reached by SOTA 
     time_needed = zeros(1,Nexe);
%% ------------- Start of the execution loop -------------------------------
for exe=1:Nexe
%% ----- Change the seed of the Mersenne Twister Random Generator ----------
rng shuffle
%% ------------------------ Start recording time ---------------------------
ti = 0;
tic
%% ------------  variable used during the calcul --------------------------
     % The population 
     population = [];
     % Vector conaiting the fitness of all individuals of the population 
     fit = [];
     % Vector containing the best individual in the population
     gbest = [];
     % variable containing the fitnes of best individual 
     fitness_best = 1000000000000000000000000000000000000000000000000000;
     % Vector containing the second best individual in the population
     gbest_sec = [];
     % variable containing the fitnes of the second best individual 
     fitness_best_sec = fitness_best;  % I assigned the value of gbest just to avoid wasting space
     % Vector containing the  previous best solution 
     prev_gbest = [];
     % variable containing the fitness of previous best solution
     prev_fitness_best = fitness_best; % I assigned the value of gbest just to avoid wasting space  
     % Vector containing the  previous second best solution 
     prev_gbest_sec = [];
     % variable containing the fitness of previous second best solution
     prev_fitness_best_sec = fitness_best_sec; % I assigned the value of gbest just to avoid wasting space
     % save all the best fitnesses obtained through the iterations
     ALL_FITNESSES = [];
     % save the iterations number in each execution
     ALL_ITTERATIONS = [];
     % The probability vector |0.5|0.5|0.5|0.5|0.5|
     probability_vector = ones(1,Dimension)*0.5;
     
     
%% --------------------------- itteration loops ---------------------------------------
for it=1:iter 
    %% ------------ Ininitialising variables: population, fitness vector--------------
    population = rand(indiv,Dimension);
    fit = zeros(1,indiv);
    %% ------------ Start the PBIL process -------------------------------------------
    for  w=1:indiv
        %% ------------------ creating the population --------------------------------
        population(w,:) =  population(w,:) < probability_vector;
        %% ------------------ Evaluate the population --------------------------------
        fit(w)  = RC_Function(population(w,:),Dimension,network,neighbourhoud);
    end
    %% -- find the best individual (its index, fitness, etc) ----------
    [fitness_best,indbest] = min(fit);
    gbest = population(indbest,:);
    %% -- find the second best individual (it's index, fitness, etc) --
    [fitness_best_sec,indbest_sec] = min(fit(fit>min(fit)));
    gbest_sec = population(indbest_sec,:);
    %% -- make sure we are using the best indiviudal of all iterations ---
    % this is not mentioned in the article but if we do not do it, we will
    % end up having an up-down plot of the fitness_best fitness (not normal
    % for a metaheuristic ref: #Alba do not reinvent the wheel)
    if fitness_best > prev_fitness_best
        fitness_best = prev_fitness_best;
        gbest = prev_gbest;
    else
        prev_fitness_best = fitness_best;
        prev_gbest = gbest;
    end
    %% -- make sure we are using the second best indiviudal of itterations ---
    % this is not mentioned in the article but if we do not do it, we will
    % end up having an up-down plot of the fitness_best_sec fitness (not normal 
    % for a metaheuristic ref: #Alba do not reinvent the wheel)
    if fitness_best_sec > prev_fitness_best_sec
        fitness_best_sec = prev_fitness_best_sec;
        gbest_sec = prev_gbest_sec;
    else
        prev_fitness_best_sec = fitness_best_sec;
        prev_gbest_sec = gbest_sec;
    end
    %% ------------- updating the probability centeral vector ----------------
    % new_lr = Lr/2; we use this strategy as the one used in the original source (see email clarification.pdf)
    probability_vector = (probability_vector * (1 - new_lr)) + (gbest * new_lr); 
    probability_vector = (probability_vector * (1 - new_lr)) + (gbest_sec * new_lr); 
    %% ------------- Mutate the probability central evctor --------------------
    rand_vect = rand(1,Dimension);
    pm_vect = ones(1,Dimension) * Pm;
    cond_vect = rand_vect < pm_vect;
    mut_ind = find (cond_vect == 1);
    probability_vector(mut_ind) =  (probability_vector(mut_ind)*(1 - Ms)) + (rand * Ms);
    %% Recording the number of itterations needed to obtain results better/== SOTA
    if fitness_best <= state_of_the_art(ind)
        if iter_needed(exe) == 0
            iter_needed(exe)  = it;
            time_needed(exe)  = toc;
        end
    end
    %% ------------------------------------------------------------------------
    ALL_FITNESSES = [ALL_FITNESSES fitness_best];
    ALL_ITTERATIONS = [ALL_ITTERATIONS (it*indiv)];
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
end
%% I added this command because Daniel told me that if i don't add it it will ot escape and display the results of the run 
exit;
