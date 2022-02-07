import numpy as np
import matplotlib.pyplot as plt

# TODO
# import own matrices (coords) option
# try surface plots
# function for plotting comparisons
# README - how to use, description of algorithm


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
def translate_to_origin(matrix):
    

    matrix_means = np.mean(matrix, axis=1)
    matrix_translated = matrix.copy()
    for axis in range(len(matrix_means)):
        matrix_translated[axis] -= matrix_means[axis]

    return matrix_translated, matrix_means


def compute_covariance_matrix(A_translated, B_translated):
    
    # multiply Bt (transposed) by A (gives H)
    H = np.matmul(B_translated.T, A_translated)

    return H 


def compute_optimal_rotation_matrix(H):
    
    # find SVD of H 
    U, S, V = np.linalg.svd(H) # np.linalg.svd does not return transpose of V
    Vt = V.T
    # keep matrices 1 and 3 from the decomposition (U and V) - will be Vt if python svd gives transposed V
    # rotation matrix R is Vt x Ut (transposed U) (transpose V if not given as transposed)
    R = np.matmul(Vt, U.T)
    
    return R


def apply_rotation(A_translated, R):

    # to get rotated A, multiply A x R
    A_rotated = np.matmul(A_translated, R)

    return A_rotated


# 4 translate A back to where B originally was centered (add back averages to matrix columns)
def revert_translation(A_rotated, B_translated, A_means, B_means):
    
    A_reverted = A_rotated.copy()
    for axis in range(len(A_means)):
        A_reverted[axis] += A_means[axis]

    B_reverted = B_translated.copy()
    for axis in range(len(B_means)):
        B_reverted[axis] += B_means[axis]

    return A_reverted, B_reverted


def run_kabsch():

    A, B = generate_matrices()
    A_translated, A_means = translate_to_origin(A)
    B_translated, B_means = translate_to_origin(B)
    H = compute_covariance_matrix(A_translated, B_translated)
    R = compute_optimal_rotation_matrix(H)
    A_rotated = apply_rotation(A_translated, R)
    A_reverted, B_reverted = revert_translation(A_rotated, B_translated, A_means, B_means)


    # plot A and B
    ax = plt.figure(figsize=(15,15)).add_subplot(projection='3d')
    ax.plot(A[0], A[1], A[2], label='A', color='r')
    ax.plot(B[0], B[1], B[2], label='B', color='b')
    ax.plot(A_translated[0], A_translated[1], A_translated[2], label='A translated', color='r', linestyle='dashed')
    ax.plot(B_translated[0], B_translated[1], B_translated[2], label='B translated', color='b', linestyle='dashed')
    ax.legend(loc='upper left', fontsize=20)
    ax.set_xlabel('X', fontsize=20)
    ax.set_ylabel('Y', fontsize=20)
    ax.set_zlabel('Z', fontsize=20)
    plt.savefig('A_B_original_vs_translated.png')
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

    # plot A vs A translated vs A rotated vs A reverted
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
    ax5.plot(B[0], B[1], B[2], label='B original', color='b')
    ax5.plot(B_translated[0], B_translated[1], B_translated[2], label='B translated', color='r')
    ax5.plot(B_reverted[0], B_reverted[1], B_reverted[2], label='B reverted', linestyle='dashed', color='y')
    ax5.legend(loc='upper left', fontsize=20)
    ax5.set_xlabel('X', fontsize=20)
    ax5.set_ylabel('Y', fontsize=20)
    ax5.set_zlabel('Z', fontsize=20)
    plt.savefig('B_comparison.png')
    plt.close()

    # plot A vs A reverted
    ax = plt.figure(figsize=(15,15)).add_subplot(projection='3d')
    ax.plot(A_reverted[0], A_reverted[1], A_reverted[2], label='A reverted')
    ax.plot(A[0], A[1], A[2], label='A original')
    ax.legend(loc='upper left', fontsize=20)
    ax.set_xlabel('X', fontsize=20)
    ax.set_ylabel('Y', fontsize=20)
    ax.set_zlabel('Z', fontsize=20)
    plt.savefig('A_vs_A_reverted.png')
    plt.close()

    # plot A reverted vs B 
    ax = plt.figure(figsize=(15,15)).add_subplot(projection='3d')
    ax.plot(A_reverted[0], A_reverted[1], A_reverted[2], label='A reverted')
    ax.plot(B[0], B[1], B[2], label='B original')
    ax.legend(loc='upper left', fontsize=20)
    ax.set_xlabel('X', fontsize=20)
    ax.set_ylabel('Y', fontsize=20)
    ax.set_zlabel('Z', fontsize=20)
    plt.savefig('B_vs_A_reverted.png')
    plt.close()

    # plot A vs B 
    ax = plt.figure(figsize=(15,15)).add_subplot(projection='3d')
    ax.plot(A[0], A[1], A[2], label='A original')
    ax.plot(B[0], B[1], B[2], label='B original')
    ax.legend(loc='upper left', fontsize=20)
    ax.set_xlabel('X', fontsize=20)
    ax.set_ylabel('Y', fontsize=20)
    ax.set_zlabel('Z', fontsize=20)
    plt.savefig('B_vs_A_original.png')
    plt.close()

    return 

run_kabsch()