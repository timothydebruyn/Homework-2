
# region imports
import copy as cp  # a quick way to access deepcopy through an alias
# endregion

# region functions
def FirstNonZero_Index(R):
    """
    Finds pivot for a row (i.e., first non-zero number in a row reading from left to right)
    :param R: a row vector
    :return: the column index (start counting at zero) of the first non-zero number
    """
    for ColumnIndex in range(len(R)):
        if R[ColumnIndex] != 0.0:
            return ColumnIndex  # found a non-zero element in the row, so return the column index
    return -1  # if the whole row is zeros, returns -1

def MakeDiagDom(A):
    """
    Reorders the rows of matrix A to put the largest absolute values along the diagonal.
    :param A: The matrix to sort
    :return: The sorted matrix
    """
    m = len(A)
    for i in range(m):
        max_row = max(range(i, m), key=lambda r: abs(A[r][i]))  # Find row with the largest absolute value in column i
        if max_row != i:
            A = SwapRows(A, i, max_row)  # Swap rows
    return A

# region row operations
def SwapRows(A, r1, r2):
    '''
    One of the elementary row operations in Gaussian elimination.
    :param A: A matrix
    :param r1: index of row 1
    :param r2: index of row 2
    :return: The A matrix after the row swap is done.
    '''
    rmax = max(r1, r2)  # the larger index
    rmin = min(r1, r2)  # the smaller index
    RMax = A[rmax]  # temporarily store the row vector at rmax in a variable
    RMin = A.pop(rmin)  # pop function removes this row from the matrix and shifts all larger indices by -1
    A.insert(rmin, RMax)  # insert the row RMax at the location rmin. All higher index rows increase index by +1
    A[rmax] = RMin  # now, replace row rmax with RMin
    return A  # done

def MultRow(R, s=1):
    '''
    Used to multiply a row vector by a scalar value
    :param R: the row vector
    :param s: the scalar with default value = 1
    :return: a new row vector multiplied by the scalar (s*R)
    '''
    return [r * s for r in R]

def AddRows(R1, R2, s=1.0):
    '''
    Adds a scalar multiple of row vector R2 to row vector R1.
    R2 and R1 must be the same length
    :param R1: a row vector
    :param R2: another row vector
    :param s: a scalar
    :return: a new row vector (R1+s*R2)
    '''
    return [r1 + s * r2 for r1, r2 in zip(R1, R2)]

# endregion

# the Echelon form of a matrix is when I produce an upper triangular matrix by Gaussian elimination
def EchelonForm(A):
    '''
    I'm expecting a Matrix of m rows by n columns.
    This function performs row operations (Gauss elimination) to produce echelon form matrix.
    :param Matrix: the matrix
    :return: the echelon form of the matrix
    '''
    m, n = len(A), len(A[0])  # number of rows and columns of A
    Ech = cp.deepcopy(A)  # make a deep copy of A so that I don't actually change A

    # order the rows by first non-zero in each column
    for i in range(m):  # iterate through all rows
        for r in range(i, m):  # iterate through all rows below row i
            p = FirstNonZero_Index(Ech[r])  # find column index in row r that is non-zero
            if p == i:  # found a row with non-zero in ith position
                Ech = SwapRows(Ech, r, i)  # move this row to the ith row
                break  # stops iterating through rows below i since I found a suitable row to put in position i
        if Ech[i][i] != 0.0:  # if I found a non-zero value for the [i][i] pivot
            for r in range(i + 1, m):  # now add multiples of row i to rows i+1 to m in order to make column i values zero below row i
                p = FirstNonZero_Index(Ech[r])
                if p == i:  # found row p has a nonzero element in column i
                    s = -Ech[r][p] / Ech[i][i]
                    Ech[r] = AddRows(Ech[r], Ech[i], s)
    return Ech

def ReducedEchelonForm(A):
    """
    This functions first creates an echelon form matrix from A and then calculates a reduced echelon form of A
    by subsequent row operations.
    :param A: The matrix to work on
    :return: The reduced echelon form of the matrix A
    """
    REF = EchelonForm(A)  # first reduce to echelon form
    for i in range(len(A) - 1, -1, -1):  # iterate from last row to row 0
        R = REF[i]
        j = FirstNonZero_Index(R)  # find the first non-zero column in r
        R = MultRow(R, 1.0 / R[j])  # make jth value equal to 1.0
        REF[i] = R
        for ii in range(i - 1, -1, -1):  # remember, end index in range is non-inclusive
            RR = REF[ii]
            if RR[j] != 0:
                RR = AddRows(RR, R, -RR[j])
                REF[ii] = RR
    return REF

# produce and identity matrix of the same size as A
def IDMatrix(A):
    '''
    Create and return an identity matrix of same dimensions as A
    :param A:
    :return:
    '''
    m, n = len(A), len(A[0])  # number of rows and columns
    return [[1 if j == i else 0 for j in range(n)] for i in range(m)]

# produce an augmented matrix from A and B
def AugmentMatrix(A, B):
    '''
    Create an augmented matrix from two matrices
    :param A: a matrix
    :param B: another matrix
    :return:
    '''
    C = cp.deepcopy(A)
    for i in range(len(C)):
        C[i] += B[i]  # this is called concatenating a list
    return C

# remove the jth column from matrix A
def popColumn(A, j):
    '''
    I want to remove column j from matrix A. I'm using slicing to cut out the column j
    :param A: The matrix
    :param j: Index of the column I want to remove
    :return: The matrix with column j removed
    '''
    AA = cp.deepcopy(A)
    for row in AA:
        row.pop(j)
    return AA

def insertColumn(A, b, i):
    '''
    This should insert column vector b into matrix A at index i. All columns to the right of i should move right by 1
    :param A: a matrix
    :param b: a column vector
    :param i: the index where to insert b
    :return: the new matrix with b inserted
    '''
    ANew = cp.deepcopy(A)
    for r in range(len(ANew)):
        ANew[r].insert(i, b[r])
    return ANew

def replaceColumn(A, b, i):
    '''
    This replaces a column of A with column vector b at column index i
    :param A: a matrix
    :param b: a column vector
    :param i: the column index of column to replace
    :return: a new matrix with the new column
    '''
    ANew = cp.deepcopy(A)
    ANew = popColumn(ANew, i)
    ANew = insertColumn(ANew, b, i)
    return ANew


def InvertMatrix(A):
    """
    Finds the inverse of matrix A by forming the augment matrix AI and using Gauss elimination
    to move the identity matrix to the left yielding IAinv, where Ainv is the inverse matrix
    :param A: the matrix to invert
    :return: the inverted matrix
    """
    ID = IDMatrix(A)
    Ainv = AugmentMatrix(A, ID)
    IAinv = ReducedEchelonForm(Ainv)
    return [row[-len(ID[0]):] for row in IAinv]  # extract the inverse from the augmented matrix

# use this to multiply matrices of correct dimensions
def MatrixMultiply(A, B):
    '''
    For multiplication of matrices, I need mXn * nXp to give a mXp matrix.
    So, must first check number of cols of A equals number of rows of B.
    Then, do matrix multiplication.
    :param A: A mxn matrix
    :param B: A nxp matrix
    :return: A matrix of shape mxp
    '''
    m, n = len(A), len(A[0])
    nn, p = len(B), len(B[0])
    if n != nn:
        raise ValueError(f"Cannot multiply: A's columns {n} must match B's rows {nn}")

    C = [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(p)] for i in range(m)]
    return C


def main():
    M = [[4, -1, -1, 3], [-2, -3, 1, 9], [-1, 1, 7, -6]]
    print("Original matrix:")
    for r in M:
        print(r)

    E = EchelonForm(M)
    print("\nEchelon form:")
    for r in E:
        print(r)

    RREF = ReducedEchelonForm(M)
    print("\nReduced Echelon Form:")
    for r in RREF:
        print(r)

    # for solving [A][x] = [b]
    A = popColumn(M, len(M[0]) - 1)  # remove last column of augmented matrix M

    MI = InvertMatrix(A)
    print("\nInverted Matrix:")
    for r in MI:
        print(r)

    B = MatrixMultiply(A, MI)
    print("\nA^-1 * A:")
    for r in B:
        print(r)


if __name__ == "__main__":
    main()
