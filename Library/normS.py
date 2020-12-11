def normS(A):
    #val = normS(A) returns the norm of a second order tensor A (stress type)
    if A.shape == (3,3):
        #tensor represented as 3x3
        A2 = np.zeros((3,3))
        A2 = sum(A[i]*A[i] for i in range(0,len(A)))
        val = np.sqrt(sum(A2))
    elif A.shape == (6,1) or A.shape == (6,) or A.shape == (1,6):
        #tensor represented as 6x1 or 1x6
        A = A.reshape(6,1)
        A2 = np.zeros(6).reshape(6,1)
        for i in range(0,len(A)):
            A2[i] = A[i]*A[i] 
        val = np.sqrt(sum(A2[i] for i in range(0,3)) + 2*sum(A2[i] for i in range(3,6)))
    elif A.shape == (3,1) or A.shape == (3,) or A.shape == (1,3):
        #tensor represented as 3x1 or 1x3
        A = A.reshape(3,1)
        A2 = np.zeros(3).reshape(3,1)
        for i in range(0,len(A)):
            A2[i] = A[i]*A[i]
        val = np.sqrt(sum(A2[i] for i in range(0,3)))
    else:
        print('Unsupported representation of second order tensor in normS()')
        val = 0
    return val
