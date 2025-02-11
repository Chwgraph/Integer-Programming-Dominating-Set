# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:41:18 2025

@author: howiech
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import DOMIP

# Function to read a matrix from user input
def read_matrix():
    matrix = []
    rows = int(input("Enter the number of rows: "))
    for i in range(rows):
        row = input(f"Enter row {i + 1} (space-separated values): ")
        # Split the input string into a list of integers
        matrix.append(list(map(int, row.split())))
    return matrix, rows

#function to generate a random adjacency matrix
def generate_random_adjacency_matrix():
    # Generate a random matrix with values 0 or 1
    n=np.random.randint(5,100)
    matrix = np.random.randint(0, 2, size=(n,n))

    # Make it symmetric to ensure it's undirected
    matrix = np.tril(matrix) + np.tril(matrix, -1).T

    # Set diagonal to 0 (no self-loops)
    np.fill_diagonal(matrix, 0)

    return matrix, n

def generate_Z(k):
    return np.random.randint(0,2,k)

# Read the matrix or generate a random adjacency matrix
# matrix, rows = read_matrix()
matrix, rows = generate_random_adjacency_matrix()
Z=generate_Z(rows)
print(Z)

params = {
"WLSACCESSID": 'xxxxxxxxxx',
"WLSSECRET": 'xxxxxxxxx',
"LICENSEID": xxxxxxxx,
}

d1=DOMIP.dom(matrix, rows,'CP',params,10000)
# var=d1.MDS(5, 's')
var=d1.CMDS('MAR')
# var=d1.PMDS(Z)
if var!=[]:
    d1.draw(var)
#var=d1.DRDP(4)
#if var!=[]:
 #   d1.draw_DRDP(var)
