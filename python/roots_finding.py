import numpy as np
import math
from scipy import optimize as opt
import dispersion as disp

def free_solid_free_fk_real(frequency,wavenumber,nlayers,vp,vs,density,thickness):
    val=disp.free_solid_free(frequency,wavenumber,nlayers,vp,vs,density,thickness)
    val=val.real
    return val

def free_solid_free_cf_real(phase_velocity,frequency,nlayers,vp,vs,density,thickness):
    val=disp.free_solid_free(frequency,2*math.pi*frequency/phase_velocity,nlayers,vp,vs,density,thickness)
    val=val.real
    return val

def rigid_solid_rigid_fk_real(frequency,wavenumber,nlayers,vp,vs,density,thickness):
    val=disp.rigid_solid_rigid(frequency,wavenumber,nlayers,vp,vs,density,thickness)
    val=val.real
    return val

def rigid_solid_rigid_cf_real(phase_velocity,frequency,nlayers,vp,vs,density,thickness):
    val=disp.rigid_solid_rigid(frequency,2*math.pi*frequency/phase_velocity,nlayers,vp,vs,density,thickness)
    val=val.real
    return val

def free_solid_halfspace_fk_real(frequency,wavenumber,nlayers,vp,vs,density,thickness):
    val=disp.free_solid_halfspace(frequency,wavenumber,nlayers,vp,vs,density,thickness)
    val=val.real
    return val

def free_solid_halfspace_cf_real(phase_velocity,frequency,nlayers,vp,vs,density,thickness):
    val=disp.free_solid_halfspace(frequency,2*math.pi*frequency/phase_velocity,nlayers,vp,vs,density,thickness)
    val=val.real
    return val

def halfspace_solid_halfspace_fk_real(frequency,wavenumber,nlayers,vp,vs,density,thickness):
    val=disp.halfspace_solid_halfspace(frequency,wavenumber,nlayers,vp,vs,density,thickness)
    val=val.real
    return val

def halfspace_solid_halfspace_cf_real(phase_velocity,frequency,nlayers,vp,vs,density,thickness):
    val=disp.halfspace_solid_halfspace(frequency,2*math.pi*frequency/phase_velocity,nlayers,vp,vs,density,thickness)
    val=val.real
    return val

def fkroots(func, fmin, fmax, df, wavenumber, nlayers, vp, vs, density, thickness):
    """
    Compute all roots of a characteristic function by frequency sweep between fmin and fmax for a given wavenumber
    
    Input:
        func: function - characteristic function
        fmin,fmax,df: float - minimum, maximum and step for frequency sweep
        wavenumber: float
        nlayers: integer - number of layers in the system (including half-spaces)
        vp, vs, density, thickness: 1D numpy arrays - elastic properties and thicknesses for the nlayers (any random thickness for half-spaces) 

    Output:
        1D numpy array containing all the frequency roots
    """
    roots = np.zeros(100)
    nf=round((fmax-fmin)/df)
    nr=0
    r=0
    for i in range(nf):
        f1=i*df+fmin
        f2=f1+df
        if (func(f1,wavenumber,nlayers,vp,vs,density,thickness)*func(f2,wavenumber,nlayers,vp,vs,density,thickness)<=0):
            r = opt.root_scalar(func,args=(wavenumber,nlayers,vp,vs,density,thickness),bracket=[f1,f2])
            roots[nr]=r.root
            nr+=1
    roots = np.resize(roots,nr)
    return roots

def cfunc(func):
    def newfunc(frequency,phase_velocity,nlayers,vp,vs,density,thickness):
        return func(phase_velocity,frequency,nlayers,vp,vs,density,thickness)
    return newfunc

def cfroots(func, fmin, fmax, df, phase_velocity, nlayers, vp, vs, density, thickness):
    """
    Compute all roots of a characteristic function by frequency sweep between fmin and fmax for a given phase velocity
    
    Input:
        func: function - characteristic function
        fmin,fmax,df: float - minimum, maximum and step for frequency sweep
        c: float
        nlayers: integer - number of layers in the system (including half-spaces)
        vp, vs, density, thickness: 1D numpy arrays - elastic properties and thicknesses for the nlayers (any random thickness for half-spaces) 

    Output:
        1D numpy array containing all the frequency roots
    """
    roots = np.zeros(100)
    nf=round((fmax-fmin)/df)
    nr=0
    r=0
    c=phase_velocity
    func2=cfunc(func)
    for i in range(nf):
        f1=i*df+fmin
        f2=f1+df
        if (func(c,f1,nlayers,vp,vs,density,thickness)*func(c,f2,nlayers,vp,vs,density,thickness)<=0):
            r = opt.root_scalar(func2,args=(phase_velocity,nlayers,vp,vs,density,thickness),bracket=[f1,f2])
            roots[nr]=r.root
            nr+=1
    roots = np.resize(roots,nr)
    return roots

def fk_dispersion_curve(func, frequency, wavenumber, kmax, dk, nlayers, vp, vs, density, thickness):
    """
    Trace a dispersion curve in the f-k domain starting from a given f-k root

    Input:
        func: function - characteristic function
        frequency, wavenumber: float - initial f-k root from which to trace the dispersion curve
        kmax, dk: float - maximum wavenumber and step for the dispersion curve
        nlayers: integer - number of layers in the system (including half-spaces)
        vp, vs, density, thickness: 1D numpy arrays - elastic properties and thicknesses for the nlayers (any random thickness for half-spaces) 

    Output:
        2D numpy array containing all the f-k roots (the dispersion curve)
    """
    nk = round((kmax-wavenumber)/dk)
    nk=int(nk)
    fk = np.zeros((2,nk))
    k = wavenumber
    f = frequency
    fk[0,0] = k
    fk[1,0] = f

    # compute the first 5 roots
    k += dk
    r,infodict,ier,mesg = opt.fsolve(func,args=(k,nlayers,vp,vs,density,thickness),x0=fk[1,0],full_output=True)
    if ier==1:
        fk[:,1] = [k,r]
    else:
        fk = fk[0:2,0:0]
        return fk

    k += dk
    r,infodict,ier,mesg = opt.fsolve(func,args=(k,nlayers,vp,vs,density,thickness),x0=2*fk[1,1]-fk[1,0],full_output=True)
    if ier==1:
        fk[:,2] = [k,r]
    else:
        fk = fk[0:2,0:1]
        return fk

    k += dk
    r,infodict,ier,mesg = opt.fsolve(func,args=(k,nlayers,vp,vs,density,thickness),x0=2*fk[1,2]-fk[1,1],full_output=True)
    if ier==1:
        fk[:,3] = [k,r]
    else:
        fk = fk[0:2,0:2]
        return fk

    k += dk
    r,infodict,ier,mesg = opt.fsolve(func,args=(k,nlayers,vp,vs,density,thickness),x0=2*fk[1,3]-fk[1,2],full_output=True)
    if ier==1:
        fk[:,4] = [k,r]
    else:
        fk = fk[0:2,0:3]
        return fk

    k += dk
    r,infodict,ier,mesg = opt.fsolve(func,args=(k,nlayers,vp,vs,density,thickness),x0=2*fk[1,4]-fk[1,3],full_output=True)
    if ier==1:
        fk[:,5] = [k,r]
    else:
        fk = fk[0:2,0:4]
        return fk

    # compute the rest of the roots using quadratic extrapolation with alternate points
    ik = 5
    for ik in range(6,nk,1):
        k += dk
        f0=fk[1,ik-6]-3*fk[1,ik-4]+3*fk[1,ik-2]
        r,infodict,ier,mesg = opt.fsolve(func,args=(k,nlayers,vp,vs,density,thickness),x0=f0,full_output=True)
        if ier==1:
            fk[:,ik] = [k,r]
        else:
            print(mesg)
            break

    fk = fk[0:2,0:ik]

    return fk

def cf_dispersion_curve(func, phase_velocity, frequency, fmax, df, nlayers, vp, vs, density, thickness):
    """
    Trace a dispersion curve in the c-f domain starting from a given c-f root

    Input:
        func: function - characteristic function
        phase_velocity, frequency: float - initial c-f root from which to trace the dispersion curve
        fmax, df: float - maximum frequency and step for the dispersion curve
        nlayers: integer - number of layers in the system (including half-spaces)
        vp, vs, density, thickness: 1D numpy arrays - elastic properties and thicknesses for the nlayers (any random thickness for half-spaces) 

    Output:
        2D numpy array containing all the c-f roots (the dispersion curve)
    """
    nf = round((fmax-frequency)/df)
    nf=int(nf)
    cf = np.zeros((2,nf))
    f = frequency
    c = phase_velocity
    cf[:,0] = [f,c]

    # compute the first 5 roots
    f += df
    r,infodict,ier,mesg = opt.fsolve(func,args=(f,nlayers,vp,vs,density,thickness),x0=cf[1,0],full_output=True)
    if ier==1:
        cf[:,1] = [f,r]
    else:
        cf = cf[0:2,0:0]
        return cf

    f += df
    r,infodict,ier,mesg = opt.fsolve(func,args=(f,nlayers,vp,vs,density,thickness),x0=2*cf[1,1]-cf[1,0],full_output=True)
    if ier==1:
        cf[:,2] = [f,r]
    else:
        cf = cf[0:2,0:1]
        return cf

    f += df
    r,infodict,ier,mesg = opt.fsolve(func,args=(f,nlayers,vp,vs,density,thickness),x0=2*cf[1,2]-cf[1,1],full_output=True)
    if ier==1:
        cf[:,3] = [f,r]
    else:
        cf = cf[0:2,0:2]
        return cf

    f += df
    r,infodict,ier,mesg = opt.fsolve(func,args=(f,nlayers,vp,vs,density,thickness),x0=2*cf[1,3]-cf[1,2],full_output=True)
    if ier==1:
        cf[:,4] = [f,r]
    else:
        cf = cf[0:2,0:3]
        return cf

    f += df
    r,infodict,ier,mesg = opt.fsolve(func,args=(f,nlayers,vp,vs,density,thickness),x0=2*cf[1,4]-cf[1,3],full_output=True)
    if ier==1:
        cf[:,5] = [f,r]
    else:
        cf = cf[0:2,0:4]
        return cf

    # compute the rest of the roots using quadratic extrapolation with alternate points
    jf = 5
    for jf in range(6,nf,1):
        f += df
        c0=cf[1,jf-6]-3*cf[1,jf-4]+3*cf[1,jf-2]
        r,infodict,ier,mesg = opt.fsolve(func,args=(f,nlayers,vp,vs,density,thickness),x0=c0,full_output=True)
        if ier==1:
            cf[:,jf] = [f,r]
        else:
            print(mesg)
            break

    cf = cf[0:2,0:jf]

    return cf