This folder contains the files responsible for the interference and measurement.
I had a problem which was that the number of entries returned was not equal to the number of quantum circuits I was sending.
The problem was that I was appending the result of each circuit to a dictionary that already exists. 
The problem with the predefined option "append" is that it represents each entry only once. Meaning: if we have two times a measurement with "10001", it will be represented only once.