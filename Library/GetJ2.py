def GetJ2(A): 
    if A.shape == (3,3):
        #This is just to remain consistent with the MATLAB file, will revamp
        m = 1/3 * (A[0,0]+A[1,0]+A[2,0])
        hold = np.zeros((6,1))
        hold[0] = A[0,0] - m 
        hold[1] = A[1,0] - m
        hold[2] = A[2,0] - m
        out = 0.5 * (hold[0]*hold[0] + hold[1]*hold[1] + hold[2]*hold[2] + 2*(A[0,1]*A[0,1] + A[1,1]*A[1,1] + A[2,1]*A[2,1]))
        out = out[0]
    elif A.shape == (6,) or A.shape == (6,1) or A.shape == (1,6):    
        #J2 of 6x1 and 1x6 inputs
        A = A.reshape(6,1) #force 6x1 dimension for calculation
        m = 1/3 * (A[0] + A[1] + A[2])
        hold = np.zeros((6,1))
        hold[0] = A[0] - m 
        hold[1] = A[1] - m
        hold[2] = A[2] - m
        out = 0.5 * (hold[0]*hold[0] + hold[1]*hold[1] + hold[2]*hold[2] + 2*(A[3]*A[3] + A[4]*A[4] + A[5]*A[5]))
        out = out[0]
    return out

