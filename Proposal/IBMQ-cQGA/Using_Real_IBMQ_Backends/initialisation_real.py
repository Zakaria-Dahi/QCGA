# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 20:27:12 2019
last time checked: 28-06-2020

@author: Zakaria Dahi
@Insitution: Zakaria Dahi

this file perform the initialisation in real mode 

The intialisation of a  Quantum Circuit of "qbits" Qbits and executed  "shots" times
qbits: the # of qbits to be generated
shots: the number of individual in the population
exe: correspond to the number of the API token to recover
"""

""" Before I used to import compile also but it seems it is deprecated """

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, IBMQ
from qiskit.compiler import transpile, assemble

   
def initialisation_real_Py_Func(shots,qbits,exe):
    # the integer when coming from matlab is considered float: need to convert it to "int"
    exe = int(exe)
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
    # recover the 16 qubit melbourne
    provider = IBMQ.get_provider(group='open')
    backend = provider.get_backend('ibmq_16_melbourne')
    
    """ -=-=- this is a deprectaed block of code to execute a job -=-=-=-
    
    # execute the quantum circuit
    qobj = compile(qc, backend=backend, shots=shots)
    job = backend.run(qobj)
    # get the results of the execution
    result = job.result()
    counts = result.get_counts(qc)
    """
    
    """ -=-=- this is a new block of code to execute a job -=-=-=-"""

    mapped_circuit = transpile(qc, backend=backend)
    qobj = assemble(mapped_circuit, backend=backend, shots=shots)
    job = backend.run(qobj)
    result_sim = job.result()
    counts = result_sim.get_counts()  
    
    # return the result of the simulation
    return(counts);
