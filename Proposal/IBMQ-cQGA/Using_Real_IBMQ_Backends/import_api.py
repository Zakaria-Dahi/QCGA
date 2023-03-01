# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 10:50:46 2020
last time checked: 28-06-2020

@author: Zakaria Abdelmoiz Dahi
Before creating this file, I was importing and loading the api token every time I call a python function: initialisation and interference
solution: I made this file to load the credentials once for all at the begining of the execution
Avoided issue: this prevent having bugs of "credentials already loaded"
"""

#import the API token
from qiskit import IBMQ
import api

def import_api_func():
    # desable active account : in case we are sure there is another account in use
    #IBMQ.disable_account()
    # enable the wanted account
    IBMQ.enable_account(api.apitoken)
    return();