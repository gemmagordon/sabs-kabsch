import numpy as np

# https://zpl.fi/aligning-point-patterns-with-kabsch-umeyama-algorithm/
# https://pymolwiki.org/index.php/Kabsch#:~:text=The%20Kabsch%20algorithm%20uses%20linear,implementation%20of%20the%20Kabsch%20algorithm.
# https://stackoverflow.com/questions/60877274/optimal-rotation-in-3d-with-kabsch-algorithm
# https://towardsdatascience.com/the-definitive-procedure-for-aligning-two-sets-of-3d-points-with-the-kabsch-algorithm-a7ec2126c87e


# NOTE import own data? set import_matrices arg with default option if no data
# NOTE docstrings
# NOTE calculate RMSD? 
# NOTE generate visualisation of data? 


# 1 initialise the two matrices as np.arrays
def import_matrices():

    # create two random matrices
    A = np.random.random((3,10))
    B = np.random.random((3,10))

    return A, B


# 2 translate A and B so their centroids are on the origin (0,0,0) - do this by subtracting from each element the average of the whole column
def translate_to_origin(A, B):

    # find and store mean value for each column (x, y, z)
    Ax_mean, Ay_mean, Az_mean = np.mean(A, axis=1)
    Bx_mean, By_mean, Bz_mean = np.mean(B, axis=1)

    # A and B can be split into 3 cols of x,y,z co-ords
    Ax, Ay, Az = A
    Bx, By, Bz = B
    
    # set up lists of cols and means to iterate through
    columns = [Ax, Ay, Az, Bx, By, Bz] 
    means = [Ax_mean, Ay_mean, Az_mean, Bx_mean, By_mean, Bz_mean]

    # set up empty matrices to hold updated values
    translated_matrices = np.empty(shape=(6,10))
    translated_column = np.empty(shape=(1,10))

    # iterate through column values and substract column mean for each col in each matrix
    for column, mean in zip(columns, means):
        for coord in column:
            translated_coord = coord - mean
            np.append(translated_column, translated_coord)


    np.append(translated_matrices, translated_column)
    
    # separate results back into the two respective matrices A and B 
    A_translated = translated_matrices[0:3, 0:10]
    B_translated = translated_matrices[3:6, 0:10]


    return A_translated, B_translated, means



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

def revert_translation(means, B_translated, A_rotated):

    # A and B can be split into 3 cols of x,y,z co-ords
    Ax, Ay, Az = A_rotated
    Bx, By, Bz = B_translated

    columns = [Ax, Ay, Az, Bx, By, Bz]

    # set up empty matrices to hold updated values
    reverted_matrices = np.empty(shape=(6,10))
    reverted_column = np.empty(shape=(1,10))

    # iterate through column values and add column mean for each col in each matrix
    for column, mean in zip(columns, means):
        for coord in column:
            reverted_coord = coord + mean
            np.append(reverted_column, reverted_coord)

    np.append(reverted_matrices, reverted_column)
    
    # separate results back into the two respective matrices A and B 
    A_reverted = reverted_matrices[0:3, 0:10]
    B_reverted = reverted_matrices[3:6, 0:10]

    return A_reverted, B_reverted


def run_kabsch():

    A, B = import_matrices()
    print('A shape: ', A.shape)
    print('B shape: ', B.shape)
    A_translated, B_translated, means = translate_to_origin(A, B)
    H = compute_covariance_matrix(A_translated, B_translated)
    print('H shape: ', H.shape)
    R = compute_optimal_rotation_matrix(H)
    print('R shape: ', R.shape)
    A_rotated = apply_rotation(A_translated, R)
    print('A_rotated shape: ', A_rotated.shape)
    A_reverted, B_reverted = revert_translation(means=means, B_translated=B_translated, A_rotated=A_rotated)
    print('A_reverted shape: ', A_reverted.shape)
    print('B_reverted shape: ', B_reverted.shape)

    return A_reverted, B_reverted

run_kabsch()


