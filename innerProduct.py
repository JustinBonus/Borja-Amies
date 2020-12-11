def innerProduct( X, Y, type):
    #inner = X.*Y for type = stress, strain, or neutral storage
    inner = 0.0
    if X.shape == (3, 3):
        #tensor represented as 3x3
        inner = np.multiply(X,Y) # Outputs (3,3) array
    elif X.shape == (6,1) or X.shape == (6,) or X.shape == (1,6):
        X = X.reshape(6,1)
        Y = Y.reshape(6,1)
        #Pseudo-switch
        switcher = {
            1: 2.0, #Stress
            2: 0.5, #Strain
            3: 1.0  #stress:strain
        }
        modifier = switcher.get(type, "Invalid type in innerProduct()")
        for i in range(0, len(X)):
            inner = inner + X[i]*Y[i]
            inner = inner[0] # Output scalar
            if i > 2:
                inner = inner + (modifier - 1.0)*X[i]*Y[i]
                inner = inner[0] # Output scalar
    elif X.shape == (3,1) or X.shape == (3,) or X.shape == (1,3):
        X = X.reshape(3,1)
        Y = Y.reshape(3,1)
        for i in range(0, len(X)):
            inner = inner + X[i]*Y[i]
            inner = inner[0]
    else:
        print('Unsupported representation of second order tensor in innerProduct()')
    return inner
