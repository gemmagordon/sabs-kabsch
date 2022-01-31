import matplotlib.pyplot as plt
import scipy as sp
import numpy as np


# https://zpl.fi/aligning-point-patterns-with-kabsch-umeyama-algorithm/
# https://pymolwiki.org/index.php/Kabsch#:~:text=The%20Kabsch%20algorithm%20uses%20linear,implementation%20of%20the%20Kabsch%20algorithm.
# https://stackoverflow.com/questions/60877274/optimal-rotation-in-3d-with-kabsch-algorithm
# https://towardsdatascience.com/the-definitive-procedure-for-aligning-two-sets-of-3d-points-with-the-kabsch-algorithm-a7ec2126c87e

# STEP 1  - 2 sets of coordinates/points to be aligned
# STEP 2 - find centroids 
# STEP 3 - align by centroids
# STEP 4 - compute rotation matrix
# STEP 5 - calculate scale factor 
# STEP 6 - translate back to original position

# input 2 paired sets of n points with 3 co-ordinates i.e. (Nx3 matrix)
# find centroids
# find variance to calculate scaling
# calculate covariance matrix
# compute SVD - H = UDV^T
# detect and prevent reflection with another matrix
# calculate optimal rotation matrix R (using SVD) and scale factor 


# 1 initialise the two matrices as np.arrays
def import_matrices():

    A = [[]]
    B = [[]]

    return A, B


# 2 translate A and B so their centroids are on the origin (0,0,0) - do this by subtracting from each element the average of the whole column
def translate_to_origin(A, B):


    return A_translated, B_translated


# 3 find optimal rotation matrix that turns A as close as possible to B (minimise RMSD) & apply that rotation to A
def compute_covariance_matrix(A_translated, B_translated):
    # multiply Bt (transposed) by A (gives H)
    H = np.matmul(B_translated.T, A_translated)

    return H 


def compute_optimal_rotation_matrix(H):
    # find SVD of H 
    U, S, V = np.linalg.svd(H) # np.linalg.svd returns transpose of V
    Vt = V
    # keep matrices 1 and 3 from the decomposition (U and V) - will be Vt if python svd gives transposed V
    # rotation matrix R is Vt x Ut (transposed U) (transpose V if not given as transposed)
    R = np.matmul(Vt, U.T)

    return R


def apply_rotation(A_translated, R):

    # to get rotated A, multiply A x R
    A_rotated = np.matmul(A_translated, R)

    return A_rotated


# 4 translate A back to where B originally was centered (add back averages to matrix columns)
def revert_translation():

    return A