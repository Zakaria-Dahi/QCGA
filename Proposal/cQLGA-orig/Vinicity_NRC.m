function [Final_Result] = Vinicity_NRC(Actual_NRC,Actual_RC,Actual_neighbourhoud,History_NRC_ID)
%% ---- Compute the vinicity of the non reporting cells ------------------------------
%% ----- Extract the first level of neighbouring cells ---------------
Actual_vinicity = 0;
Actual_NRC_ID   = [];
Non_reporting_cells_Neighbourhoud = [];
%% ----- Extract the ID of the first level of neighbouring cells ---------------
% --- This loop is used to brows the 6 canonical neighbouring cell of each Cell in a network ----
for i=2:7 
    % -- Variable used to stock the response of the following check -----------------
    check = 0;
    % Verify if the current neighbour is not already visited ---------------
    for w=1:length(History_NRC_ID)
        if Actual_neighbourhoud(Actual_NRC,i) == History_NRC_ID(w)
            check = 1;
        end
    end
    % --- Record the cell as a neighbour only if it is not recorded yet --------
    if check == 0 
        % ------ This condition here is to exclude counting the 999 which is in reality not a cell ----    
        if (Actual_neighbourhoud(Actual_NRC,i) ~= 1000) && (Actual_neighbourhoud(Actual_NRC,i) ~= Actual_NRC) 
            % --  Test if the actual neighbouring cell is already a reporting cell ------
            Response =  isempty(find((Actual_neighbourhoud(Actual_NRC,i) ==  Actual_RC) == 1));
            if Response == 1
                        %% ----- Recording the Non reporting cells related to each reporting cell --------            
                        Non_reporting_cells_Neighbourhoud = [Non_reporting_cells_Neighbourhoud  Actual_neighbourhoud(Actual_NRC,i)];
                        %% ----- Recording the Non reporting cells related to each reporting cell --------
                        Actual_NRC_ID = [Actual_NRC_ID Actual_neighbourhoud(Actual_NRC,i)];
                        %% ----- Increasing the Vinicity Value Each Time  A non reporting Cell belong to the neighbourhood of the actual Reporting Cell ----
                        Actual_vinicity = Actual_vinicity + 1;
                        % --- Add the new visited cells to the one already visite  -----------------------------
                        History_NRC_ID = [History_NRC_ID Actual_neighbourhoud(Actual_NRC,i)];
            end
        end
    end
end
%% ---- Extracting the 1,2,3 ... N level of non reporting cells ----------------------------
% ----  Variable to stock the Number of NRC of the actual RC  -------------------------
if isempty(Non_reporting_cells_Neighbourhoud) == 0
    Size_NRC = length(Non_reporting_cells_Neighbourhoud);
    for j=1:Size_NRC
        result = 0;
        result = Vinicity_NRC(Non_reporting_cells_Neighbourhoud(j),Actual_RC,Actual_neighbourhoud,History_NRC_ID);
        Actual_vinicity = Actual_vinicity +  result(1); 
        History_NRC_ID = result(2:length(result));
    end
end
Final_Result = [Actual_vinicity History_NRC_ID];
end

