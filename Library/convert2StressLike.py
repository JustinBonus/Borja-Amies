def convert2StressLike(A):
    #stress = normS(A) .. returns a stress like second order tesnor A
    #                  (in stress type storage)
    if A.shape == (3,3):
        #tensor represented as 3x3 matrix
        hold = A
    elif A.shape == (6,1) or A.shape == (6,) or A.shape == (1,6):
        #tensor represented as 6x1 or 1x6
        A = A.reshape(6,1)
        hold = np.zeros((6,1))
        hold[0] = A[0]
        hold[1] = A[1]
        hold[2] = A[2]
        hold[3] = 0.5*A[3]
        hold[4] = 0.5*A[4]
        hold[5] = 0.5*A[5]
    else:
        print('Unsupported representation of second order tensor in convert2StressLike()')
    return hold
