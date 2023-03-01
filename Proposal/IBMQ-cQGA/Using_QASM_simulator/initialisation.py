# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 20:27:12 2019

@author: Zakaria Dahi
@Insitution: Zakaria Dahi

this file perform the initialisation in simulation mode 

The intialisation of a  Quantum Circuit of "qbits" Qbits and executed  "shots" times
qbits: the # of qbits to be generated
shots: the number of individual in the population
"""

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, Aer, execute
   
def initialisation_Py_Func(shots,qbits):
    # the integer when coming from matlab is considered float: need to convert it to "int"
    shots = int(shots)
    qbits = int(qbits)
    #create a quantum circuit(score): with quantum register nd a classical register    
    q = QuantumRegister(qbits)
    c = ClassicalRegister(qbits)
    qc = QuantumCircuit(q, c)
    #Apply hadamard gate and measurement on each qbit in the quantum register
    #1/sqrt(2) initialisation: 1/2 probability to get 0 or 1
    # create the superposition of states
    for i in range(qbits):
        #hadamard gate: create superposition
        qc.h(q[i])
        # measurement
        qc.measure(q[i],c[i])
    # enable the QASM simulator    
    backend = Aer.get_backend('qasm_simulator')
    # execute the quantum circuit
    job_sim = execute(qc, backend, shots=shots)
    sim_result = job_sim.result()
    # get the results of the simulation
    results_m = sim_result.get_counts(qc)
    # return the result of the simulation
    return(results_m);