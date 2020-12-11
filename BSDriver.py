def BSDriver(LoadCase):
    # BoundingSurface J2 with kinematic hardening 
    # Written by Pedro Arduino, Mar. 22 2019
    # Copyright Arduino Computational Geomechanics Group
    # Ported into Python/Jupyter Notebook by Justin Bonus, Jul. 2019
    #
    #
    # LoadCase:
    #    1 ... proportionally increasing strain
    #    2 ... cyclic strain
    #    3 ... proportionally increasing stress
    #    4 ... cyclic stress
    #
    # ======  LOADING CASES ==================================================
    
    import numpy as np
    from collections import namedtuple
    
    nPoints = 200

    ## Switch for LoadCases:
    ## Pseudo-switch created by using python dictionary to hold LoadCase functions
    def case_one():
        case_one.time   = np.linspace(0,1,nPoints+1)
        case_one.strain = np.array([ 0.05, -0.015, -0.015, 0.000, 0.000, 0.000 ]).reshape(6,1) * case_one.time
        case_one.StressDriven = 0
        return case_one
    def case_two():
        nCycles = 3
        omega   = 0.15
        case_two.time   = np.linspace(0,nCycles*2*np.pi/omega,nCycles*nPoints+1);
        case_two.strain = np.array([ 0.00, -0.000, -0.000, 0.045, 0.000, 0.000 ]).reshape(6,1) * np.sin( omega*case_two.time )      
        case_two.StressDriven = 0 
        return case_two
    def case_three():
        case_three.time   = np.linspace(0,1,nPoints+1)       
        case_three.stress = np.array([[0.100],
                           [0.000],
                           [0.000],
                           [0.000],
                           [0.000],
                           [0.000]])*case_three.time + 0.0*np.array([1,1,1,0,0,0]).reshape(6,1)*np.ones( case_three.time.shape )            
        case_three.StressDriven = 1    
        return case_three
    def case_four():
        nCycles = 3
        omega   = 0.15
        case_four.time   = np.linspace(0, nCycles*2*np.pi/omega, nCycles*nPoints+1)
        case_four.stress = np.array([[0.000],
                           [0.000],
                           [0.000], #.01, .03, -.01, .05, 0, -.02
                           [0.050],
                           [0.000],
                           [0.000]])*np.sin( omega*case_four.time ) + 0.0*np.array([1,1,1,0,0,0]).reshape(6,1)*np.ones( case_four.time.shape )            
        case_four.StressDriven = 1          
        return case_four

    case_switcher = {
        1: case_one,
        2: case_two,
        3: case_three,
        4: case_four
    }    

    case = case_switcher.get(LoadCase, lambda: "Invalid LoadCase")
    case() #Runs the LoadCase function. Creates: case.time, case.strain | case.stress, case.StressDriven
    time, StressDriven = case.time, case.StressDriven 
    if StressDriven:
        stress = case.stress
        strain = np.zeros((6,1)) #initialize empty 6x1 strain numpy array for stress-driven scenario
    else:
        strain = case.strain
        stress = np.zeros((6,1)) #initialize empty 6x1 stress numpy array for strain-driven scenario
    
    Stress0 = np.zeros((6,1)) #Initialize first 'unloading' point
    StrainDriven = int(not StressDriven)

    # ========================================================================
    # ---- MATERIAL PARAMETERS
    # Static Parameters

    # Static Parameters
    E = 20 #Elastic Modulus  MPa
    v= 0.49 #Poissons ratio, less than 0.5 to allow compresibility
    G = E/(2*(1+v)) #Shear modulus
    K = E/(3*(1-2*v)) #Bulk modulus
    Kmod = 0 #Isotropic Hardening
    Su = 0.061 #Yield stress in 1-D tension test MPa
    hh = G #kinematic hardening parameter
    mm = 1.0 #kinematic hardening parameter
    beta = 0.5 #midpoint integration
    RR = np.sqrt(8/3)*Su

    #namedtuple used to organzie related variables, similar to a structure
    static = namedtuple('StaticParam',['E','v','G','K','Kmod','Su','hh','mm','beta','RR'])
    StaticParam = static(E,v,G,K,Kmod,Su,hh,mm,beta,RR)


    # ========================================================================
    # ---- INITIAL CONDITIONS

    # Initialize the state variables
    if StrainDriven:
        IniStress = -0.0*(np.array([1, 1, 1, 0, 0, 0]).reshape(6,1))
        IniStrain = np.linalg.solve(GetCe(StaticParam), IniStress) #Check if GetCe compacts to nxn
    elif StressDriven:
        IniStress =  0.0*(np.array([1, 1, 1, 0, 0, 0]).reshape(6,1))
        IniStrain =  0.0*(np.array([1, 1, 1, 0, 0, 0]).reshape(6,1))   

    #Structure for IniState (initial state parameters, static) and CurState (changing state parameters)
    state = namedtuple('state', ['eP','alphaISO','Stress0', 'Kappa', 'Psi'])
 
    eP = 0.0*(np.array([1, 1, 1, 0, 0, 0]).reshape(6,1))
    alphaISO = 0.0  
    Stress0 = 0.0*(np.array([1, 1, 1, 0, 0, 0]).reshape(6,1))
    Kappa = 0.0
    Psi = 0.0
    IniState = state(eP, alphaISO, Stress0, Kappa, Psi)

    # For first iteration
    CurStress = IniStress
    CurStrain = IniStrain
    CurState  = IniState

    # Variables used for plotting
    alphaISO_plot, j2_plot, j2e_plot, stress_var_plot, stress_var2_plot = [], [], [], [], [] #Initiliaze list format
    alphaISO_plot.append(0) #Python list allows for easy data addition
    strain[:,0] = CurStrain.T - IniStrain.T 
    stress[:,0] = CurStress.T
    j2_plot.append(0)
    j2e_plot.append(0)
    stress_var_plot.append(0)
    Stress0[:,0] = CurStress.T
    Iter = np.zeros(time.shape)


    # ========================================================================
    # ---- COMPUTATION CYCLES

    if StrainDriven:
        #StrainDriven
        for i in range(1, (len(strain[0]) )):

            NextStrain = strain[:,i] + IniStrain.T
            dStrain = strain[:,i] - strain[:, i-1] #Driving variable
            
            #Current BSRadialMap is a function, will be transformed into a class eventually
            NextStress, NextState, NextCep = BSRadialMap(dStrain, StaticParam, CurStress, CurState)

            # Update Stress, Strain, and State
            CurStress = NextStress
            CurState = NextState
            
            # Variables created for plotting purposes
            alphaISO_plot.append(CurState.alphaISO)
            stress = np.append(stress, CurStress, 1)
            j2_plot.append(GetJ2(CurStress))
            stress_var_plot.append(np.sqrt(2*j2_plot[i])*np.sqrt(3/2)*np.sign(stress[0,i] - stress[1,i]))
            stress_var2_plot.append((stress[0,i] - stress[1,i]))
            Stress0 = np.append(Stress0, CurState.Stress0, 1)
            
    elif StressDriven:
        # StressDriven driver
        # set tolerance value for iterative procedure(s)
        TOLERANCE = 1e-10 

        for i in range(0, len(stress[0])-1):

            # initialize strain epsilon_{n+1}^{(0)} = eps_{n} using the old state
            # (this is the initial approximation for eps_{n+1}
            if i == 0:
                # special settings for initial values at t_1
                NextStrain = np.array([0,0,0,0,0,0]).reshape(6,1)
                dStrain = np.array([0,0,0,0,0,0]).reshape(6,1)
                CurState = IniState
            else:
                NextStrain = CurStrain
                dStrain = np.array([0,0,0,0,0,0]).reshape(6,1)

            NextStress, NextState, Cep = BSRadialMap(dStrain, StaticParam, CurStress, CurState)

            RR = stress[:, i].reshape(6,1) - NextStress
            RR = RR.reshape(6,1)
            RR0 = normS(RR)

            # reset iteration counter
            kk = 0
            # iterate until convergence
            while normS(RR)/RR0 > TOLERANCE:
                
                # update strain from eps_{n+1}^{(k)} to eps_{n+1}^{(k+1)}
                dStrain = np.linalg.solve(Cep, RR)
                NextStrain = NextStrain + dStrain

                # compute material response for estimated strain state
                # NOTE: the state variables are taken at t_n
                NextStress, NextState, Cep = BSRadialMap(dStrain, StaticParam, CurStress, CurState)
                #print('NextStress:',NextStress)
                #print('Stress0:',NextState.Stress0)
                # check for equilibrium
                RR = stress[:,i].reshape(6,1) - NextStress
                RR = RR.reshape(6,1)
                kk = kk + 1
                # emergence exit if procedure does not converge            
                if kk > 3:
                    print('procedure slow to converge. Error : ', normS( RR )/RR0)
                
                if kk > 20:
                    print('procedure did not converge. Error : ', normS( RR )/RR0)
                    print('YOUR TANGENT Cep IS WRONG', normS( RR )/RR0)
                    break

                Iter[i] = kk
                CurStress = NextStress
                CurState = NextState


            # Update State variables for next step
            CurStress = NextStress
            CurStrain = NextStrain
            CurState  = NextState

            # Update variables for plotting purposes
            strain = np.append(strain, CurStrain, 1)
            alphaISO_plot.append(CurState.alphaISO)
            j2_plot.append(GetJ2(CurStress))
            stress_var_plot.append(np.sqrt(2*j2_plot[i])*np.sqrt(3/2)*np.sign(stress[3,i]))
            Stress0 = np.append(Stress0, CurState.Stress0, 1)
                                   
    DriverOutput = namedtuple('DriverOutput',['StaticParam','time','strain','stress','alphaISO','j2','stress_var','stress_var2', 'Stress0','Iter'])
    DriverOutput = DriverOutput(StaticParam, time, strain, stress, alphaISO_plot, j2_plot, stress_var_plot, stress_var2_plot, Stress0, Iter)
    
    return DriverOutput
    
    # =========================================================================
