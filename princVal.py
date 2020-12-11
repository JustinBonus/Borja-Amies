def princVal(vec, type):
    #Based on the code of Robert Siegwart, 2019
    #Edited by Justin Bonus, July 2019
    #type = 0 returns normal principals
    #type = 1 returns shear principals
    import numpy as np
    
    if vec.shape == (3,3):
        S = vec
    elif vec.shape == (6,1) or vec.shape == (6,) or vec.shape == (1,6):  
        vec = vec.reshape(6,1)
        S = np.array([[vec[0], vec[3], vec[5]],
                       [vec[3], vec[1], vec[4]],
                       [vec[5], vec[4], vec[2]] ])
        S = S.reshape(3,3)
    else:
        print('Unsupported vector shape in principal()')
    
    e_val, e_vec = np.linalg.eig(S) #Solve for eigenvalues and eigenvectors
    p3, p2, p1 = np.sort(e_val)   #Sort smallest to largest
    #p1, p2, p3 = e_val
    
    if type == 1:
        p1 = (p1+p3)/2
        p2 = (p1+p2)/2
        p3 = (p2+p3)/2

    return p1, p2, p3
