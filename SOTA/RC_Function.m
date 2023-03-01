function [result] = RC_Function(individual,Dimension,network,neighbourhoud)
%% --------  Computing the Location updates related to every reporting cells of the network ------
% --------- Variable that stock the constant introduced in most of the Papers ------
Beta = 10;
% -- Variable to stock the location update cost ------------------------
Lu = 0;
% -- Variable to Stock The Reporting Cells ID --------------------------
Reporting_cells = [];
% -- Variable to stock the Reporting Cell's Vinicity ------------------
Vinicity_vector = [];
Non_reporting_cells_vector = [];
Paging_Vinicity_Data = [];
%% -------------- Compute The location Update and Extract The RC ID ------
for cell=1:Dimension
    if individual(cell) == 1
        %% ----- Compute The Location Update For the Reporting Cells -----
        Lu = Lu + network(cell,1);
        %% ----- Extract The Reporting Cells ID --------------------------
        Reporting_cells = [Reporting_cells cell];
    end
end
%% ------- Computing the vinicity for the reporting cell of the network ---------------------
% -------- Variable to stock the number of the RCs Found ----------------
Size_RC = length(Reporting_cells);
for RC=1:Size_RC % --- For Every Reporting Cell Do : ---------------------------
    % -- When processing each reporting cell this vector is used to save temporarly the ID of non reporting cells related to the neghbourhood of the actual reporting cell -------------
    Non_reporting_cells = [];
    NRC_ID = [];
    % ------ I've Incremented the vinicity by 1 from the start to not forget to count the cell itself --- 
    vinicity = 1;
    %% ----- Extract the ID of the first level of neighbouring cells ---------------
    for i=2:7
      % ------ This condition here is to exclude counting the 999 which is in reality not a cell ---- 
      if neighbourhoud(Reporting_cells(RC),i) ~= 1000
          % --  Test if the actual neighbouring cell is already a reporting cell ------
          Response =  isempty(find((neighbourhoud(Reporting_cells(RC),i) ==  Reporting_cells) == 1));
          if Response == 1
              %% ----- Recording the Non reporting cells related to each reporting cell --------            
              Non_reporting_cells = [Non_reporting_cells  neighbourhoud(Reporting_cells(RC),i)];
              %% ----- Recording the Non reporting cells related to each reporting cell --------
              NRC_ID = [NRC_ID neighbourhoud(Reporting_cells(RC),i)];
              %% ----- Increasing the Vinicity Value Each Time  A non reporting Cell belong to the neighbourhood of the actual Reporting Cell ----
              vinicity = vinicity + 1;
          end
      end
    end
    %% ---- Extracting the 1,2,3 ... N level of non reporting cells ----------------------------
    % ----  Variable to stock the Number of NRC of the actual RC  -------------------------
    if isempty(Non_reporting_cells) == 0
     Size_NRC = length(Non_reporting_cells);
        for j=1:Size_NRC
             result = 0;
             result = Vinicity_NRC(Non_reporting_cells(j),Reporting_cells,neighbourhoud,NRC_ID);
             vinicity = vinicity +  result(1); 
             NRC_ID = result(2:length(result));
        end
    end  
% ---- This vector is used to stock the vinicity of the non reporting cells related to the treated reporting cell at this time 
Vinicity_vector = [Vinicity_vector vinicity];
%% -- Record The Non-Reporting  Cells Related to Each one Of the Reporting-Cells
%  -- The vector is with the format x, NRC where X is the number of non-reporting-cells and NRC are the ID of the NRC
% -- I added recovering the size because the 2017 version uses an additional input ---
[rr,cc] = size(NRC_ID);
Non_reporting_cells_vector = [Non_reporting_cells_vector mat2cell(NRC_ID,rr)];
end
Vinicity_Reporting_Cells = [Reporting_cells ;Vinicity_vector];
%% ----------- Extracting the Vinicity of the non-reporting cells -------------------
Non_reporting_cell_ID = [];
Non_reporting_cell_Vinicity = [];
for cell=1:Dimension
    if individual(cell) == 0
        %% ----- Extract The ID of the Non-Reporting Cell -----
           Non_reporting_cell_ID = [Non_reporting_cell_ID cell];
    end
end
%% -- Comparing the number of the non-reporting cell to each one of the neighbourhood of the reporting ones -------
for k=1:length(Non_reporting_cell_ID)
    vinicity_record = [];
    for z=1:length(Non_reporting_cells_vector)
        Neighbour_RC = cell2mat(Non_reporting_cells_vector(z));
        Response =  (isempty(find((Non_reporting_cell_ID(k) == Neighbour_RC) == 1)));
        if Response == 0
           vinicity_record = [vinicity_record Vinicity_Reporting_Cells(2,z)];
        end
    end
    Non_reporting_cell_Vinicity = [Non_reporting_cell_Vinicity max(vinicity_record)]; 
end
Vinicity_Non_Reporting_Cells = [Non_reporting_cell_ID;Non_reporting_cell_Vinicity];
%% --------- Computing The paging cost * Vinicity -----------------------------------
Paging_Vinicity_Data = [Vinicity_Non_Reporting_Cells Vinicity_Reporting_Cells];
Pc = 0;
for id_one=1:Dimension
    for id_two=1:Dimension
        if Paging_Vinicity_Data(1,id_two) == id_one
           Pc = (network(id_one,2) * Paging_Vinicity_Data(2,id_two)) + Pc;
        end
    end
end
%% ------- Affecting The results ----------------------------------------------------
result = (Beta * Lu) + Pc;
end


