def vectorAngle(vec1, vec2, type):
    # First take the dev() and hydroProjector() of the vectors if measuring
    # angle from deviatoric view
    #
    #if vec1.shape == (6,1) or vec1.shape(6,) or vec1.shape == (1,6):
    #theta = cos-1((u dot v)/(||u|| dot ||v||))
    switch_type = {
        1: 2, #Stress
        2: 0.5, #Strain
        3: 1 #Stress Strain
    }
    factor = switch_type.get(type, 'Invalid type')
    angle = float(np.arccos(innerProduct(vec1, vec2,1)/(normS(vec1)*normS(vec2)))*(180/np.pi)) #Degrees     
    #else:
    #angle = 'Incorrect input vector shape'
    return angle
