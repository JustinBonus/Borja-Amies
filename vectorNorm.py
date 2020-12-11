def vectorNorm(X, type):
    #norm(X) in type = stress, strain or neutral storage format
    vNorm = np.sqrt(innerProduct(X,X,type))
    return vNorm
