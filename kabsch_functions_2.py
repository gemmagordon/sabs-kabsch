import numpy as np
import matplotlib.pyplot as plt


def generate_matrices():

    # generate x, y, z co-ords and concat into matrix
    theta = np.linspace(10, 20, 50)

    A_z = np.linspace(10, 20, 50)
    r = A_z**2 + 1
    A_x = r * np.sin(theta)
    A_y = r * np.cos(theta)

    B_z = np.linspace(10, 20, 50)
    q = B_z**2 + 1
    B_x = q * np.cos(theta) # switch sin and cos for x and y 
    B_y = q * np.sin(theta)

    # extract points into matrix
    A = np.array([A_x, A_y, A_z], dtype='float32')
    B = np.array([B_x, B_y, B_z], dtype='float32')

    return A, B


# 2 translate A and B so their centroids are on the origin (0,0,0) - do this by subtracting from each element the average of the whole column
def translate_to_origin(A, B):
    
    # find and store mean value for each column (x, y, z)
    Ax_mean, Ay_mean, Az_mean = np.mean(A, axis=1)
    Bx_mean, By_mean, Bz_mean = np.mean(B, axis=1)

    # A and B can be split into 3 cols of x,y,z co-ords
    A_x, A_y, A_z = A
    B_x, B_y, B_z = B
    
    # set up lists of cols and means to iterate through
    columns = [A_x, A_y, A_z, B_x, B_y, B_z] 
    means = [Ax_mean, Ay_mean, Az_mean, Bx_mean, By_mean, Bz_mean]

    # set up empty matrices to hold updated values
    translated_matrices = [] # NOTE THESE ARE NOT EMPTY if using np.empty! 
    print('EMPTY MATRIX: ', translated_matrices)

    # iterate through column values and substract column mean for each col in each matrix
    subtract = lambda column: column - mean
    subtract_func = np.vectorize(subtract)

    count = 0
    for column, mean in zip(columns, means):
        count += 1 
        print(count)
        print(mean)
        print('ORIGINAL COLUMN:', column[0:3])
        translated_column = subtract_func(column)
        print('TRANSLATED COLUMN: ', translated_column[0:3])
        translated_matrices.append(translated_column) # NOTE these cols are not being appended to the new matrix 

    translated_matrices = np.column_stack((translated_matrices[0], translated_matrices[1], translated_matrices[2], 
                                            translated_matrices[3], translated_matrices[4], translated_matrices[5]))
    print('FILLED IN MATRIX: ', translated_matrices)

    # separate results back into the two respective matrices A and B 
    A_translated = translated_matrices[0:3, 0:50]
    B_translated = translated_matrices[3:6, 0:50]

    return A_translated, B_translated, means


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
    Ax_rot, Ay_rot, Az_rot = A_rotated
    Bx_rot, By_rot, Bz_rot = B_translated

    rotated_columns = [Ax_rot, Ay_rot, Az_rot, Bx_rot, By_rot, Bz_rot]

    # set up empty matrices to hold updated values
    reverted_matrices = []
    print('EMPTY REVERTED MATRIX: ', reverted_matrices)

    # iterate through column values and substract column mean for each col in each matrix
    add = lambda rotated_column: rotated_column + mean
    add_func = np.vectorize(add)

    for rotated_column, mean in zip(rotated_columns, means):
        print('ROTATED COLUMN: ', rotated_column[0:3])
        reverted_column = add_func(rotated_column)
        print('REVERTED COLUMN: ', reverted_column[0:3]) 
        reverted_matrices.append(reverted_column)
    
    reverted_matrices = np.column_stack(( reverted_matrices[0], reverted_matrices[1], reverted_matrices[2], 
                                            reverted_matrices[3], reverted_matrices[4], reverted_matrices[5]))
    print('FILLED IN MATRIX: ', reverted_matrices)

    # separate results back into the two respective matrices A and B 
    A_reverted = reverted_matrices[0:3, 0:50] 
    B_reverted = reverted_matrices[3:6, 0:50]

    return A_reverted, B_reverted


def run_kabsch():

    A, B = generate_matrices()
    A_translated, B_translated, means = translate_to_origin(A=A, B=B)
    H = compute_covariance_matrix(A_translated, B_translated)
    R = compute_optimal_rotation_matrix(H)
    A_rotated = apply_rotation(A_translated, R)
    A_reverted, B_reverted = revert_translation(means=means, B_translated=B_translated, A_rotated=A_rotated)

    # plot A and B
    ax = plt.figure(figsize=(15,15)).add_subplot(projection='3d')
    ax.plot(A[0], A[1], A[2], label='A')
    ax.plot(B[0], B[1], B[2], label='B')
    ax.legend(loc='upper left', fontsize=20)
    ax.set_xlabel('X', fontsize=20)
    ax.set_ylabel('Y', fontsize=20)
    ax.set_zlabel('Z', fontsize=20)
    plt.savefig('A_B_original.png')
    plt.close()

    # plot A and B translated
    ax2 = plt.figure(figsize=(15,15)).add_subplot(projection='3d')
    ax2.plot(A_translated[0], A_translated[1], A_translated[2], label='A translated')
    ax2.plot(B_translated[0], B_translated[1], B_translated[2], label='B translated')
    ax2.legend(loc='upper left', fontsize=20)
    ax2.set_xlabel('X', fontsize=20)
    ax2.set_ylabel('Y', fontsize=20)
    ax2.set_zlabel('Z', fontsize=20)
    plt.savefig('A_B_translated.png')
    plt.close()

    # plot reverted A and B 
    ax3 = plt.figure(figsize=(15,15)).add_subplot(projection='3d')
    ax3.plot(A_reverted[0], A_reverted[1], A_reverted[2], label='A_reverted')
    ax3.plot(B_reverted[0], B_reverted[1], B_reverted[2], label='B_reverted')
    ax3.legend(loc='upper left', fontsize=20)
    ax3.set_xlabel('X', fontsize=20)
    ax3.set_ylabel('Y', fontsize=20)
    ax3.set_zlabel('Z', fontsize=20)
    plt.savefig('A_B_reverted.png')
    plt.close()

    # plot A vs A translated vs A reverted
    ax5 = plt.figure(figsize=(15,15)).add_subplot(projection='3d')
    ax5.plot(A[0], A[1], A[2], label='A original')
    ax5.plot(A_translated[0], A_translated[1], A_translated[2], label='A translated')
    ax5.plot(A_rotated[0], A_rotated[1], A_rotated[2], label='A rotated')
    ax5.plot(A_reverted[0], A_reverted[1], A_reverted[2], label='A reverted')
    ax5.legend(loc='upper left', fontsize=20)
    ax5.set_xlabel('X', fontsize=20)
    ax5.set_ylabel('Y', fontsize=20)
    ax5.set_zlabel('Z', fontsize=20)
    plt.savefig('A_comparison.png')
    plt.close()

    # plot B vs B translated vs B reverted
    ax5 = plt.figure(figsize=(15,15)).add_subplot(projection='3d')
    ax5.plot(B[0], B[1], B[2], label='B original')
    ax5.plot(B_translated[0], B_translated[1], B_translated[2], label='B translated')
    ax5.plot(B_reverted[0], B_reverted[1], B_reverted[2], label='B reverted')
    ax5.legend(loc='upper left', fontsize=20)
    ax5.set_xlabel('X', fontsize=20)
    ax5.set_ylabel('Y', fontsize=20)
    ax5.set_zlabel('Z', fontsize=20)
    plt.savefig('B_comparison.png')
    plt.close()

    return 

run_kabsch()