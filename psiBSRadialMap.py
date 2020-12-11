def psiBSRadialMap(dStr, StaticParam, CurStrain, CurStress, CurState, StressDriven):
    # Radial Return Algorithm for J2 Bounding Surface Plasticity Model
    # Original model by Borjas & Amies 1994
    # Base radialmap by Pedro Arduino, Mar. 22 2019
    # Copyright Arduino Computational Geomechanics Group
    # Precomputed matrix enabled version by Justin Bonus, Jul. 2019
    #
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
    R = StaticParam.RR

    meye = np.eye(3,1); meye[1]=1; meye[2]=1 
    small = 1.0e-10 # Tolerance for small dStrain steps, AKA no stress change
    debugFlag = 1 #Set 0 for true, 1 for false

    #=================================
    # -----INDEX MAPPING AND ROTATIONS

    # Change these to passed values
    #res_x = 101
    #res_y = 101
    cen_x = int(np.floor(res_x/2))
    cen_y = int(np.floor(res_y/2))
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if StressDriven == 0:
        # Strain-driven loading
        #print('Strain-driven')
        dStrain = dStr

        dev_dStrain = dev(dStrain)
        vol_dStrain = trace(dStrain)

        NextStress0 = CurState.Stress0 
        dev_NextStress0 = dev(NextStress0)
        dev_CurStress = dev(CurStress)

        CurKappa = CurState.Kappa
        CurPsi = CurState.Psi
        NextKappa = CurKappa # Temporary
        NextPsi = CurPsi # Temporary

        NextalphaISO = 0

        # Unit-vectors
        hydro = np.array([[np.sqrt(1/3)],[np.sqrt(1/3)],[np.sqrt(1/3)]]) # Hydrostatic axis unit-vector
        north = np.array([[np.sqrt(2/3)],[-np.sqrt(1/6)],[-np.sqrt(1/6)]]) #\Pi-plane north unit-vector
        east = np.array([[0],[np.sqrt(1/2)],[-np.sqrt(1/2)]]) #\Pi-plane east unit-vector
        south = -north

        #=================================
        # -----UNLOADING CHECK
        numerator = innerProduct((-(1 + CurKappa)*dev_CurStress - CurKappa*(1 + CurKappa)*(dev_CurStress - dev_NextStress0)), dev_dStrain, 3)
        denominator = innerProduct(( (1 + CurKappa)*dev_CurStress - CurKappa*dev_NextStress0), (dev_CurStress - dev_NextStress0), 1)

        if np.absolute(denominator) < small:
            loadingCond = 0
        else:
            loadingCond = numerator/denominator

        if loadingCond > 0.0:
            if debugFlag == 0:
                print('Unloading Happened')
            NextStress0 = CurStress
            dev_NextStress0 = dev(NextStress0)
            loadingCond = 0
        else:
            loadingCond = 0

        unloading = True # Set to True in order to enter while loop
        #=================================
        
        if loadingCond == 0:
            # Check if zero vector (breaks angle-finders), set angles accordingly
            if np.all(dev_NextStress0 == 0):
                northUnloadingAngle = 0.
                eastUnloadingAngle = 0.
            else:
                northUnloadingAngle = angleNormal(north, dev_NextStress0)
                eastUnloadingAngle = angleNormal(east, dev_NextStress0)

            # Check if the unloading stress is west or east of centerline
            if eastUnloadingAngle > np.pi/2:
                northUnloadingAngle = -northUnloadingAngle

            #=========================
            # -----Rotate North & Map to Index Space
            #northNextStress0 = quatHydroRotator(dev_NextStress0, -northUnloadingAngle)
            #northCurStress = quatHydroRotator(dev_CurStress, -northUnloadingAngle)
            #indexNorthNextStress0 = stressToIndex(northNextStress0, res_x, res_y, R, sliced=False) # UL Index float space
            #indexNorthCurStress = stressToIndex(northCurStress, res_x, res_y, R, sliced=False) # UL Index float space
            
            indexNextStress0 = stressToIndex(NextStress0, res_x, res_y, R, sliced=False) - cen_y # O\Pi Index float space
            indexCurStress = stressToIndex(CurStress, res_x, res_y, R, sliced=False) - cen_y # O\Pi Index float space
            indexNorthNextStress0 = indexRotator(indexNextStress0, -northUnloadingAngle) + cen_y # A\Pi Index float space
            indexNorthCurStress = indexRotator(indexCurStress, -northUnloadingAngle) + cen_y # A\Pi Index float space         
            
            #=========================
            # -----Access H' Matrix
            # Find layer containing appropiate 1D North NextStress0 position

            if indexNorthNextStress0[1] == 0:
                indexNorthNextStress0[1] = 1
            layerH = cen_y - int(np.round(indexNorthNextStress0[1]))
            jH = int(np.round(indexNorthCurStress[1]))
            iH = int(np.round(indexNorthCurStress[0]))
            cartH = hardeningM[:, :, layerH] # Holds relevant H' matrix layer
            CurH = cartH[jH, iH] # Determine current H', could pass between time-step


            # ========================        
            # -----Rotate South and Access \psi Matrix
            # Will need to find the \psi layer containing
            # the appropiate 1D SouthNorth CurH position with respect to 1D North NextStress0.
            # Then search, for paired psi and NextH values to satisfy formula 

            Iter_H = 0 # Redundant, currently used for tracking random  values of interest
            CurKappa = float((CurH/hh)**(1/mm)) # Backcalculated \kappa value, depends on hardening function
            midpoint = float(cen_y - (CurKappa * (cen_y - int(np.round(indexNorthNextStress0[1]))))/(CurKappa + 1)) # 1D Midpoint of \kappa contour
            midpoint = np.array([[cen_x],[midpoint]]) # Float index space
            radius = (float(midpoint[0] - indexNorthCurStress[0])**2 + float(midpoint[1] - indexNorthCurStress[1])**2)**0.5 # Distance from contour midpoint to CurStress
            
            subG = int(cen_y - indexNorthNextStress0[1]) # Subgroup holding appropiate NextStress0; range(0, cen_y)
            curS = int(np.round(midpoint[1])) # Layer of subgroup holding CurStress position; range(0, res_y)
            Iter_H = curS # Redundant
            layerP = subG*res_y + curS # (Subgroup relative to NextStress0) + (Layer of subgroup holding CurStress)
            cartP = psiM[:,:,layerP] # Pull appropiate \psi layer into RAM
            
            adjMidpoint = midpoint - cen_y # Recenter contour midpoint
            adjIndexNorthCurStress = np.array([[0],[0]]) # Initialize 
            adjIndexNorthCurStress[0] = int(np.round(indexNorthCurStress[0])) - cen_x # O\Pi, integer
            adjIndexNorthCurStress[1] = int(np.round(indexNorthCurStress[1])) - cen_y # O\Pi, integer
            adjIndexNorthCurStress[1] = adjIndexNorthCurStress[1] - adjMidpoint[1] # Offest for relative distance
            south = np.array([[0],[1]]) # Upper-left centered, (+x = right, +y = down)
            theta = angleIndex(south, adjIndexNorthCurStress) # Angle between south and float-index CurStress
            if np.isnan(theta):
                angle = 0

            indexNorthNorthCurStress = np.array([[cen_x],[int(np.round(midpoint[1] - radius))]]) # North aligned
            indexSouthNorthCurStress = np.array([[cen_x],[int(np.round(midpoint[1] + radius))]]) # South aligned
            southNorthDStrain = quatHydroRotator(dev_dStrain, -theta)  # dev_dStrain for South aligned CurStress, North aligned Stress0
            #northNorthDStrain = -southNorthDStrain # dev_dStrain for North aligned CurStress, North aligned Stress0

            # Search line for \psi and H' pairs
            rr, cc = psiSearch(southNorthDStrain, indexSouthNorthCurStress, res_x, res_y, R)    
            
            indexNextH = np.array([[0],[0]]) # Initialize array
            indexPsi = np.array([[0],[0]]) # Initialize array
            psi = 2*G # Initialize psi
            tol_R = 1e-1 # Tolerance for equation 37, major speed increase when looser
            hold_R = 1e6 # Initialize tolerance comparison eq. 37
            Iter_Psi = 0 # Reset psi-search iterations
            
            # 'stop' is a very important parameter
            # High values will slow the program/increase chance of error
            # Low values will limit the stress movement/increase error
            # Will eventually become 'predictively' set
            stop = 30 # How many points to evaluate on search line (more needed for larger steps/finer resolution)
            stop = np.linspace(0, stop-1, stop) 
            for j, i, s in zip(rr, cc, stop):
                NextH = float(cartH[j,i]) # H' value of NextStress cell
                # Skip if cell is on/outside of bounding surface
                if NextH == 0:
                    continue
                    
                NextKappa = float((NextH/hh)**(1/mm)) # Backcompute kappa from hardening function
                psi = float(cartP[j,i]) # \psi value of specific \psi-H' pair
                
                # Solve base equation, zero res_search and res_R means pair matches true function (i.e. no error)
                #res_search = (2*G)/(1 + 3*G*(((1-beta)/CurH) + (beta)/NextH)) - psi
                res_R = np.abs(np.linalg.norm(dev_CurStress + psi*dev_dStrain + NextKappa*(dev_CurStress + psi*dev_dStrain - dev_NextStress0))-R)
                Iter_Psi = Iter_Psi + 1

                # If in an acceptable range of true solution, break loop to save time
                if res_R <= tol_R:
                    hold = False
                    break

                #Hold best solution if one below tolerance isn't found
                if res_R < hold_R:
                    hold_NextKappa = NextKappa
                    hold_NextH = NextH
                    hold_psi = psi
                    hold_Iter_Psi = Iter_Psi
                    hold_R = res_R
                    hold = True

            # If tolerances were not met, give the best pair found
            if hold == True:
                NextKappa = hold_NextKappa 
                NextH = hold_NextH
                psi = hold_psi
                Iter_Psi = hold_Iter_Psi
                res_R = hold_R

        #==================================
        # -----Solve Constitutive Equation
        #NextStress = CurStress + K*vol_dStrain*meye + psi*dev_dStrain
        NextStress = CurStress + K*vol_dStrain*meye + (1/3)*psi*dev_dStrain
        NextStrain = CurStrain + vol_dStrain + dev_dStrain
        
        # Update State
        Iter_H = res_R # Redundant
        state = namedtuple('state', ['eP','alphaISO', 'Stress0', 'Iter_H', 'Iter_Psi', 'Kappa', 'H', 'Psi'])
        NextState = state(0, NextalphaISO, NextStress0, Iter_H, Iter_Psi, NextKappa, NextH, psi)
    
        return NextStrain, NextStress, NextState
    
    # -----------------------------------------------------------------------------------------------------------
    # STRESS - DRIVEN    
    elif StressDriven == 1:
        
        # Initialize state
        dev_CurStrain = dev(CurStrain)
        dev_CurStress = dev(CurStress)      
        
        dStress = dStr        
        dev_dStress = dev(dStress)
        vol_dStress = trace(dStress)
              
        NextStress = CurStress + dStress
        dev_NextStress = dev(NextStress)
        
        NextStress0 = CurState.Stress0 
        dev_NextStress0 = dev(NextStress0)
        norm_dev_NextStress0 = normS(dev_NextStress0)        
        
        # Initialize state parameters
        CurKappa = CurState.Kappa
        CurPsi = CurState.Psi
        NextKappa = CurKappa # Temporary
        NextPsi = CurPsi # Temporary
        NextalphaISO = 0
        
        north = np.array([[np.sqrt(2/3)],[-np.sqrt(1/6)],[-np.sqrt(1/6)]]) #\Pi-plane north unit-vector
        east = np.array([[0],[np.sqrt(1/2)],[-np.sqrt(1/2)]]) #\Pi-plane east unit-vector
        
        #=================================
        # -----UNLOADING CHECK
        unloaded = False
        numerator = innerProduct((-(1 + CurKappa)*dev_CurStress - CurKappa*(1 + CurKappa)*(dev_CurStress - dev_NextStress0)), dev_dStress, 1)
        denominator = innerProduct(( (1 + CurKappa)*dev_CurStress - CurKappa*dev_NextStress0), (dev_CurStress - dev_NextStress0), 1)
        small = 1e-10
        if np.absolute(denominator) < small:
            loadingCond = 0
        else:
            loadingCond = numerator/denominator

        if loadingCond > 0.0:
            # Reset unloading stress to current position
            if debugFlag == 0:
                print('Unloading Happened')
            NextStress0 = CurStress
            dev_NextStress0 = dev(NextStress0)
            norm_dev_NextStress0 = normS(dev_NextStress0) 
            loadingCond = 0
            unloaded = True
        else:
            loadingCond = 0
            unloaded = False
        
        # Check if zero vector unloading point (breaks angle-finders), set angles accordingly
        if np.all(NextStress0 == 0):
            northUnloadingAngle = 0.
            eastUnloadingAngle = 0.
        else:
            northUnloadingAngle = angleNormal(north, dev_NextStress0)
            eastUnloadingAngle = angleNormal(east, dev_NextStress0)
        # Check if the unloading stress is west or east of centerline
        if eastUnloadingAngle > np.pi/2:
            northUnloadingAngle = -northUnloadingAngle
            
        #============================================
        # ----- Map to Index Space and Rotate North
        # Rotate and map stress state to index notation
        northCurStress = quatHydroRotator(dev_CurStress, -northUnloadingAngle)
        northNextStress = quatHydroRotator(dev_NextStress, -northUnloadingAngle)
        northNextStress0 = quatHydroRotator(dev_NextStress0, -northUnloadingAngle)
        indexNorthCurStress = stressToIndex(northCurStress, res_x, res_y, R)
        indexNorthNextStress = stressToIndex(northNextStress, res_x, res_y, R)
        indexNorthNextStress0 = stressToIndex(northNextStress0, res_x, res_y, R)
        
        #indexNorthNextStress0 = np.array([[cen_x],[int(np.round(cen_y * (1 - norm_dev_NextStress0/R)))]])
        
        if indexNorthNextStress0[1] <= 0:
            indexNorthNextStress0[1] = 1

        #===========================================
        # ----- Access H' Matrix
        layer = cen_y - int(indexNorthNextStress0[1])  # Relevant H' matrix layer for NextStress0
        if unloaded == True:
            maxKappa = 1e5
            CurH = float((maxKappa*hh)**(mm))
        elif unloaded == False:
            CurH = float(hardeningM[int(indexNorthCurStress[1]),int(indexNorthCurStress[0]), layer])
        NextH = float(hardeningM[int(indexNorthNextStress[1]),int(indexNorthNextStress[0]), layer])
        NextKappa = float((NextH / hh)**(1/mm))
        trap = (1 + 3*G*(((1-beta)/CurH) + (beta/float(NextH))))
        psi = 2*G / trap
        
        #============================================
        # ----- Solve Constitutive Equation
        dev_dStrain = (dev_dStress * trap) / (2*G)
        vol_dStrain = (dStress - psi*dev_dStrain) / K
        #NextStrain = CurStrain + vol_dStrain + dev_dStrain
        NextStrain = CurStrain + vol_dStrain + 2*dev_dStrain

        Iter_H = 0 # No iterations required in stress-driving
        Iter_Psi = 0 # No iterations required in stress-driving
        
        # Update State
        state = namedtuple('state', ['eP','alphaISO', 'Stress0', 'Iter_H', 'Iter_Psi', 'Kappa', 'H', 'Psi'])
        NextState = state(0, NextalphaISO, NextStress0, Iter_H, Iter_Psi, NextKappa, NextH, psi)
    
        return NextStrain, NextStress, NextState
