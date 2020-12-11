def evalH(param, kappa):
    #s = evalH(param, kappa) returns the value of H=h*kappa^m (AKA exponential hardening moduli)
    if kappa < 0:
        s = 1.0e-10
        return s
    s = param.hh*(kappa)**param.mm
    return s
