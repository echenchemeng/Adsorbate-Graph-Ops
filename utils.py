# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 13:08:02 2022

@author: eccn3
"""

import numpy as np

def min_img_dist(cell, pos1, pos2):
    periodics = np.array([-1, 0, 1])
    distances = np.empty(27)
    
    idx = 0
    
    for i in periodics:
        for j in periodics:
            for k in periodics:
                distances[idx] = np.linalg.norm(np.sum(pos1 + np.array([i,j,k]) * cell.T, axis = 1))
                idx +=1
    return min(distances)


