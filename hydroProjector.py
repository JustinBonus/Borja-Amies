def hydroProjector(stress, type):
    # Justin Bonus, July 2019
    # Project forms of stress vectors onto a plane orthogonal to the hydrostatic axis, centered at the origin
    # U_{proj-plane} = U - ((U dot n)/(||n||^2))n
    # U = [s_11, s_22, s_33, s_12, s_23, s_13]; n = Vector normal to \pi-plane
    # U = [s_1, s_2, s_3]
    #===================================
    #Type 0 for for top 3 cells of 6x1, 1 for bottom 3 cells of 6x1, 2 for 3x1 principal form
    #===================================
    #Be careful that the correct unit vector is being used or an offset to the projected plane will occur 
    #Not a major issue when viewing figure in 'ortho' 
    import numpy as np
    #stress = stress.reshape(1,6)
    def zero(stress):
        #Projecting 3x1 normal principal vectors onto the \Pi-plane
        #Converting to 6x1 in order to use other functions that are called
        zero.vec = np.array([1,1,1,0,0,0])
        zero.norm_vec = normS(zero.vec)
        zero.unit_vec = zero.vec/zero.norm_vec
        zero.norm_unit_vec = normS(zero.unit_vec)
        zero.proj_stress = stress - ((innerProduct(stress,zero.unit_vec,1))/(zero.norm_unit_vec**2))*zero.unit_vec
        return zero.proj_stress
    def one(stress):
        one.vec = np.array([0,0,0,1,1,1])
        one.norm_vec = normS(one.vec)
        one.unit_vec = one.vec/one.norm_vec
        one.norm_unit_vec = normS(one.unit_vec)
        one.proj_stress = stress - ((innerProduct(stress,one.unit_vec,1))/(one.norm_unit_vec**2))*one.unit_vec
        return one.proj_stress
    def two(stress):
        two.vec = np.array([1,1,1])
        two.norm_vec = normS(two.vec)
        two.unit_vec = two.vec/two.norm_vec
        two.norm_unit_vec = normS(two.unit_vec)
        two.proj_stress = stress - ((innerProduct(stress,two.unit_vec,1))/(two.norm_unit_vec**2))*two.unit_vec
        return two.proj_stress        
    
    switch_type = {
        1: zero,
        2: one,
        3: two
    }
    #Pseudo-switch
    case = switch_type.get(type+1, lambda: "Invalid projection type")
    case(stress) #Runs the projection function.
    proj_stress = case.proj_stress
    return proj_stress
