def GetCe(mStatic):
    # Create the Elastic 4th Order Modulus
    # elastic stiffness tensor (matrix_hh)
    a = 2*mStatic.G 
    b = mStatic.K - a/3 # K + 4/3G
    #b = mStatic.E  / (3 - 6*mStatic.v) - a / 3;
    #mC = np.array([[np.multiply(a, np.eye(3)) + np.multiply(b, np.ones((3,3))),  np.zeros((3,3))],
    #   [np.zeros((3,3)) , np.multiply(mStatic.G, np.eye(3))]])
    mC = 3*mStatic.K*GetIvol() + 2*mStatic.G*GetIdev()     
    return mC
