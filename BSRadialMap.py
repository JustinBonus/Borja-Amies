def BSRadialMap(dStrain, StaticParam, CurStress, CurState):
    # Radial Return Algorithm for J2 Bounding Surface Plasticity Model 
    # Borja and Amies 1994
    # Written by Pedro Arduino
    # Copyright - Arduino Computational Geomechanics Group
    # March, 2019
    # Ported into Python/Jupyter Notebook by Justin Bonus
    # July, 2019
    #
    # Input:
    #    dStrain               ... Strain differential, tn to tn+1
    #    StaticParam.E         ... Young's modulus
    #    StaticParam.v         ... poissons ratio
    #    StaticParam.G         ... Shear modulus
    #    StaticParam.K         ... Bulk modulus
    #    StaticParam.hh        ... Kinematic Hardening parameter
    #    StaticParam.mm        ... Kinematic hardening parameter
    #    StaticParam.beta      ... Integration parameter (0=expl, 1=impl)
    #    StaticParam. RR       ... Bounding surface radius
    #    CurStress             ... Stress at tn
    #    CurState.eP           ... Plastic strain at tn
    #    CurState.alphaISO     ... Isotropic internal variable at tn
    ##   CurState.alphaKIN     ... Kinematic internal variablen at tn
    #    CurState.Stress0      ... Stress point at unloading
    #
    # Output:
    #    NextStress            ... Stress at tn+1
    #    CurState.eP           ... Plastic strain at tn
    #    CurState.alphaISO     ... Isotropic internal variable at tn+1
    #    CurState.alphaKIN     ... Kinematic internal variablen at tn+1
    #    Cep                   ... Consistent tangent modulus

    #Static Parameters
    G = StaticParam.G
    K = StaticParam.K
    mm = StaticParam.mm
    hh = StaticParam.hh
    beta = StaticParam.beta
    RR = StaticParam.RR

    tol_rel = 1.0e-10
    meye = np.eye(6,1); meye[1]=1; meye[2]=1 
    small = 1.0e-10
    debugFlag = 1 #Set 0 for true, 1 for false

    Ce = GetCe(StaticParam)
    NextCep = Ce

    dev_dStrain = dev(dStrain)
    vol_dStrain = trace(dStrain)
    res = np.zeros((2,1))
    
    
    NextStress0 = CurState.Stress0
    dev_NextStress0 = dev(NextStress0)
    dev_CurStress = dev(CurStress)

    norm_dev_CurStress = normS(dev_CurStress)
    norm_dev_NextStress0 = normS(dev_NextStress0)

    CurKappa = CurState.Kappa
    CurPsi = CurState.Psi
    #Assume next Kappa and Psi are current (will overwrite for loading)
    NextKappa = CurKappa
    NextPsi = CurPsi
    NextalphaISO = CurState.alphaISO

    numerator = innerProduct((-(1 + CurKappa)*dev_CurStress - CurKappa*(1 + CurKappa)*(dev_CurStress - dev_NextStress0)), dev_dStrain, 3)
    denominator = innerProduct(( (1 + CurKappa)*dev_CurStress - CurKappa*dev_NextStress0), (dev_CurStress - dev_NextStress0), 1)

    if np.absolute(denominator) < small:
        loadingCond = 0
    else:
        loadingCond = numerator/denominator


    if loadingCond > 0.0:
        if debugFlag == 1:
            print('Unloading Happened')

        NextStress0 = CurStress
        dev_NextStress0 = dev(NextStress0)
        loadingCond = 0

    if loadingCond == 0:
        # ====================================================================================================
        # unloading (or beginning of loading) just happened
        if debugFlag == 0:
            print('Initial Loading') 
        
        #Initial pseudo-elastic Psi and Kappa at unloading point
        NextPsi = 2.0*G 
        NextKappa = 1.0e6
        dev_dStrainNorm = normE(dev_dStrain)
             
        if np.absolute(dev_dStrainNorm) < small:
            NextStress = CurStress
        else:
            #dStrain deviatoric unit vector
            dev_dStrainDir = dev_dStrain/dev_dStrainNorm
            NextH = evalH(StaticParam, NextKappa)
            res[0] = NextPsi*(1.0 + 3.0*G*beta/NextH)/(2.0*G)-1.0 
            res[1] = vectorNorm((dev_CurStress+(1.0+NextKappa)*NextPsi*convert2StressLike(dev_dStrain)), 1)/RR - 1.0
            res_norm = np.sqrt(res[0]**2 + res[1]**2)

            #Initialize Newton variables
            iter_counter = 0
            iter_max = 50
            tol_mat = tol_rel*res_norm
            incVar = np.zeros((2,1))

            for iter_counter in range(1, iter_max+1):
                if debugFlag == 0: 
                    print('Iteration = ', num2str(iter_counter), '  Norm = ', num2str(res_norm))
                if res_norm < (tol_mat + small):
                    NextStress = CurStress + K*vol_dStrain*meye + NextPsi*convert2StressLike(dev_dStrain)
                    break 

                temp = (dev_CurStress + (1 + NextKappa)*NextPsi*convert2StressLike(dev_dStrain))
                temp = temp / vectorNorm(temp,1)
                Ktan = np.zeros((2,2))
                Ktan[0,0] = (1.0 + 3.0*G*beta/NextH) / (2.0*G)
                Ktan[0,1] = (-3.0*G*NextPsi*beta*mm/hh/(NextKappa**(mm+1.0))/(2.0*G))
                Ktan[1,0] = ((1.0+NextKappa) * innerProduct(temp, dev_dStrain, 3)) / RR
                Ktan[1,1] = (innerProduct(temp, NextPsi*convert2StressLike(dev_dStrain), 1)) / RR

                incVar = np.linalg.solve(Ktan, res)  #\ operator is left matrix division, minimizes AX - B

                NextPsi = NextPsi - incVar[0]
                NextKappa = NextKappa - incVar[1]

                NextH = evalH(StaticParam, NextKappa)

                #Calculate New Residual
                res[0] = NextPsi*(1.0 + 3.0*G*beta/NextH)/(2.0*G)-1.0 
                res[1] = vectorNorm(dev_CurStress + (1.0 + NextKappa)*NextPsi*convert2StressLike(dev_dStrain), 1)/RR - 1.0
                res_norm = np.sqrt(res[0]**2 + res[1]**2)

    else:

        # ====================================================================================================
        # Continuing Loading
        if debugFlag == 0: 
            print('Loading Continues')  
        
        
        CurH = evalH(StaticParam, CurKappa)
        #Pseudo-elastic assumption
        NextH = evalH(StaticParam, NextKappa)
        res[0] = NextPsi*(1.0 + 3.0*G*((1-beta)/CurH + beta/NextH))/(2.0*G)-1.0
        res[1] = vectorNorm(dev_CurStress + (1.0 + NextKappa)*NextPsi*convert2StressLike(dev_dStrain) + NextKappa*(dev_CurStress - dev_NextStress0), 1)/RR - 1.0
        res_norm = np.sqrt(res[0]**2+res[1]**2)

        #Initialize Newton variables
        iter_counter = 0
        iter_max = 50
        tol_mat = tol_rel*res_norm
        incVar = np.zeros((2,1))

        for iter_counter in range(1, iter_max+1): 
            if debugFlag == 0: 
                print('Iteration = ', num2str(iter_counter), '  Norm = ', num2str(res_norm))
            if res_norm < (tol_mat+small):
                NextStress = CurStress + K*vol_dStrain*meye + NextPsi*convert2StressLike(dev_dStrain)
                break 

            temp = (dev_CurStress + (1+NextKappa)*NextPsi*convert2StressLike(dev_dStrain) + NextKappa*(dev_CurStress - dev_NextStress0)); 
            temp = temp / vectorNorm(temp,1)
            Ktan = np.zeros((2,2))
            Ktan[0,0] = (1.0+3.0*G*((1-beta)/CurH + beta/NextH)) / (2.0*G)
            Ktan[0,1] = (-3.0*G*NextPsi*beta*mm/hh/(NextKappa**(mm+1.0))/(2.0*G))
            Ktan[1,0] = ((1.0+NextKappa) * innerProduct(temp, dev_dStrain, 3)) / RR
            Ktan[1,1] = (innerProduct(temp, dev_CurStress + NextPsi*convert2StressLike(dev_dStrain)-dev_NextStress0, 1)) / RR

            incVar = np.linalg.solve(Ktan,res) #Since Ktan is nxn: X = A\B | AX = B | X minimizes norm(AX - B)
            
            #Real next Psi, Kappa, and H
            NextPsi = NextPsi - incVar[0]
            NextKappa = NextKappa - incVar[1]
            NextH = evalH(StaticParam, NextKappa)

            #Calculate New Residual            
            res[0] = NextPsi*(1.0+3.0*G*((1-beta)/CurH + beta/NextH))/(2.0*G) - 1.0
            res[1] = vectorNorm(dev_CurStress+(1.0+NextKappa)*NextPsi*convert2StressLike(dev_dStrain)+ NextKappa*(dev_CurStress - dev_NextStress0), 1)/RR - 1.0
            res_norm = np.sqrt(res[0]**2 + res[1]**2)


    NextStress = CurStress + K*vol_dStrain*meye + NextPsi*convert2StressLike(dev_dStrain)

    # Update State
    state = namedtuple('state', ['eP','alphaISO','Stress0', 'Kappa', 'Psi'])
    NextState = state(0, NextalphaISO, NextStress0, NextKappa, NextPsi)
    
    #Get 6x6 representations of Ivol and Idev
    Ivol = GetIvol()
    Idev = GetIdev()
    #Ivol = np.array([[np.ones((3,3)), np.zeros((3,3))], [np.zeros((3,3)), np.zeros((3,3))]]).reshape(6,6) # 3Ivol
    #NOTE: NextCep is equivalent to 3KIvol + 2GIdev when NextPsi is at its elastic state: 2G
    NextCep = 3*K*Ivol + NextPsi*Idev
    
    return NextStress, NextState, NextCep
