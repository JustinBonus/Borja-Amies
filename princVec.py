def princVec(vec, type):
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
                       [vec[5], vec[4], vec[2]] ]).reshape(3,3)

    else:
        print('Unsupported vector shape in principal()')
    
    e_val, e_vec = np.linalg.eig(S) #Solve for eigenvalues and eigenvectors
    p3, p2, p1 = np.sort(e_val)   #Sort smallest to largest
    #p1, p2, p3 = e_val
    e_val_l = e_val.tolist() #Python list
    p1_index, p2_index, p3_index = e_val_l.index(p1), e_val_l.index(p2), e_val_l.index(p3)
    p1_vec, p2_vec, p3_vec = e_vec[:,p1_index], e_vec[:,p2_index], e_vec[:,p3_index]
    
    if type == 1:
        tau1 = (p1+p3)/2
        tau2 = (p1+p2)/2
        tau3 = (p2+p3)/2
        p1_vec = (p1_vec + p3_vec)/np.linalg.norm(p1_vec+p3_vec)
        p2_vec = (p1_vec + p2_vec)/np.linalg.norm(p1_vec+p2_vec)
        p3_vec = (p2_vec + p3_vec)/np.linalg.norm(p2_vec+p3_vec)

    return p1_vec, p2_vec, p3_vec
