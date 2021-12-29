% Date: 15-12-2021
% About the code:
	% - Contains the neighbourhood of each cell of the 15 cells network. 
	% - The neighbourhood is organised as follows: {ID of the cell, ID of neighbourhood}. 
	% - The ID "1000" is used for padding when the processed cell has no neighbourhing cell.


 Neighbour_matrix = [1 2 6 1000 1000 1000 1000;    
 2 1 6 7 8 3 1000;
 3 2 8 4 1000 1000 1000;
 4 3 8 9 10 5 1000;
 5 4 10 1000 1000 1000 1000 ;
 6 1 2 7 11 1000 1000;
 7 2 6 11 12 13 8;
 8 3 2 7 13 9 4;
 9 4 8 13 14 15 10;
 10 5 4 9 15 1000 1000;
 11 6 7 12 1000 1000 1000 ;
 12 11 7 13 1000 1000 1000 ;
 13 14 9 8 7 12 1000 ;
 14 15 9 13 1000 1000 1000 ;
 15 10 9 14 1000 1000 1000];
