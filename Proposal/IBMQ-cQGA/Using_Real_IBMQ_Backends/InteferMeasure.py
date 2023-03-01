# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 08:26:04 2019

@author: Zakaria Dahi
@Insitution: Zakaria Dahi

This creates a quantum circuit and apply : interference and measurement
candid: represents the angles thetat used to perform the rotation around "Z" axis
dim: represents the number of qubits
"""

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, Aer, execute
from math import pi
   
def InteferMeasure_Py_Func(cand,dim):
    # the integer when coming from matlab is considered float: need to convert it to "int"
    dim = int(dim)
    #create a quantum circuit(score): with quantum register nd a classical register    
    q = QuantumRegister(dim)
    c = ClassicalRegister(dim)
    qc = QuantumCircuit(q, c)
    #Apply on each qbit in the quantum register:
            # 1) hadamard gate,
            # 2) rotation arounz "z" using the angle theta
            # 3) hadamard gate
            # 4) measurement 
    #1/sqrt(2) initialisation: 1/2 probability to get 0 or 1
    # create the superposition of states
    for i in range(dim):
        #hadamard gate: create superposition
        qc.h(q[i])
        # convert angle in degrees to radians
        angle = (cand[i]*pi)/180
        # apply rotation around "z" axis
        qc.rz(angle,q[i])
        #hadamard gate: create superposition
        qc.h(q[i])
        # measurement
        qc.measure(q[i],c[i])
    # enable the QASM simulator    
    backend = Aer.get_backend('qasm_simulator')
    # execute the quantum circuit
    job_sim = execute(qc, backend, shots=1)
    sim_result = job_sim.result()
    # get the results of the simulation
    results_m = sim_result.get_counts(qc)
    # return the result of the simulation
    return(results_m);