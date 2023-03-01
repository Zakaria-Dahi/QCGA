# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 08:26:04 2019
Last time checked: 28-06-2020


@author: Zakaria Dahi
@Insitution: Zakaria Dahi

all of this is done in real mode

This creates a quantum circuit and apply : interference and measurement
candid: represents the angles thetat used to perform the rotation around "Z" axis
dim: represents the number of qubits
exe: the id of the execution to retrieve the correponsing api token
"""

""" Before I used to import compile also but it seems it is deprecated """

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, IBMQ
from math import pi
import numpy as np
from qiskit.compiler import transpile, assemble

def InteferMeasure_real_Py_Func(cand,dim,rows,exe):
    # declare an empty list that will contain the produced quantum circuits 
    list_qc = [];
    # declare an empty dictionary that will contain the prodcued counts
    dict_qc = {}
    # the integer when coming from matlab is considered float: need to convert it to "int"
    dim = int(dim)
    exe = int(exe)
    rows = int(rows)
    # recover the original matrix of offspring after conversion in matlab
    cand = np.reshape(np.ravel(cand, order='F'), (rows,dim), order='F')
    # Browse the offspring and create a quantum circuit for each offspring
    for x in range(rows):
        #create a quantum circuit(score): with quantum register and a classical register    
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
            angle = (cand[x][i]*pi)/180
            # apply rotation around "z" axis
            qc.rz(angle,q[i])
            #hadamard gate: create superposition
            qc.h(q[i])
            # measurement
            qc.measure(q[i],c[i])
        # append the produced quantum circuit to the list of quantum circuit
        list_qc.append(qc)
        
    # recover the 16 qubit melbourne backed
    provider = IBMQ.get_provider(group='open')
    backend = provider.get_backend('ibmq_16_melbourne')
    
    """ -=-=- this is a deprectaed block of code to execute a job -=-=-=
    
    # execute the quantum circuit
    qobj = compile(qc, backend=backend, shots=1)
    job = backend.run(qobj)
    # get the results of the execution
    result = job.result()
    counts = result.get_counts(qc)
    """

    """ -=-=- this is a new block of code to execute a job -=-=-=-"""

    mapped_circuit = transpile(list_qc, backend=backend)
    qobj = assemble(mapped_circuit, backend=backend, shots=1)
    job = backend.run(qobj)
    result_sim = job.result()
    
    # recover the counts
    for z in range(rows):
        counts = result_sim.get_counts(z)  
        # verify if the key already exists in the dictionary
        values_view = counts.keys()
        value_iterator = iter(values_view)
        first_value = next(value_iterator)
        # if the key already exist then increment its value
        if first_value in dict_qc:
            dict_qc[first_value] = dict_qc[first_value] + 1
        else:
            dict_qc.update(counts)
    
    # return the result of the simulation
    return(dict_qc);
