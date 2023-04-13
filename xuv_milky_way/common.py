import numpy as np
import pandas as pd
import pickle
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
from scipy import stats
import time
import os

#Find nearest neighbour in array
def find_nearest(array, value):
    """
    This function finds the array value nearest to the input value.

    array: array with some values from which you want to find the nearest to the input value
    value: your input value
    return: returns the nearest value from the array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

#Find nearest 2 neighbours in array
def find_2_nearest(array, value):
    """
    This function finds the 2 array values nearest to the input value.

    array: array with some values from which you want to find the 2 nearest to the input value
    value: your input value
    return: returns the 2 nearest values from the array
    """
    array = np.asarray(array)
    dif = np.abs(array - value)
    sort = np.argsort(dif)
    idx = sort[0] #Nearest neighbour
    idx2 = sort[1] #2nd nearest neighbour
    return array[idx], array[idx2]

#Interpolation between two data points
def interpolation(a, b, x):
    """
    This function gives you the 2D interpolated result between two points.

    a: array with the two x points you want to interpolate between.
    b: array with the two y points you want to interpolate between.
    x: point on the x axis at which we want to interpolate.
    return: returns the y value found from the 2D interpolation.
    """
    output = b[0] + (x - a[0]) * ((b[1] - b[0])/(a[1] - a[0]))
    return output

# Function that finds a file and returns its path
def find(name, path):
    """
    This function finds the document stated.

    name: file name (string).
    path: file path (string)
    return: if the document is found it returns the whole path name, if it is not found it returns None.
    """
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)