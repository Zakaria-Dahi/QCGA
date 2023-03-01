# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 08:26:04 2019

@author: Zakaria Dahi
@Insitution: Zakaria Dahi

all of this is done in real mode

This creates a quantum circuit and apply : interference and measurement
candid: represents the angles thetat used to perform the rotation around "Z" axis
dim: represents the number of qubits
exe: the id of the execution to retrieve the correponsing api token
"""

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, compile, IBMQ
from math import pi
   
def InteferMeasure_real_Py_Func(cand,dim,exe):
    # delete the current accounts
    #IBMQ.delete_accounts()
    # the integer when coming from matlab is considered float: need to convert it to "int"
    dim = int(dim)
    exe = int(exe)
    #create a quantum circuit(score): with quantum register nd a classical register    
    q = QuantumRegister(dim)
    c = ClassicalRegister(dim)
    qc = QuantumCircuit(q, c)
    # import the API token
    import api1
    IBMQ.enable_account(api1.apitoken)
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
    # recover the 16 qubit melbourne
    backend = IBMQ.get_backend('ibmq_16_melbourne')
    # execute the quantum circuit
    qobj = compile(qc, backend=backend, shots=1)
    job = backend.run(qobj)
    # get the results of the execution
    result = job.result()
    counts = result.get_counts(qc)
    # return the result of the simulation
    return(counts);
