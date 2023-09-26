import numpy as np
import math
from scipy import optimize as opt
import dispersion as disp

def leaky(phase_velocity,frequency, attenuation, nlayers,vp,vs,density,thickness):
    val=disp.halfspace_solid_halfspace_leaky(frequency,complex(2*math.pi*frequency/phase_velocity,attenuation),nlayers,vp,vs,density,thickness)
    return val

def leaky_log_abs(ca,frequency,nlayers,vp,vs,density,thickness):
    val=disp.halfspace_solid_halfspace_leaky(frequency,complex(2*math.pi*frequency/ca[0],ca[1]),nlayers,vp,vs,density,thickness)
    return math.log10(abs(val))

def leaky_abs(ca,frequency,nlayers,vp,vs,density,thickness):
    val=disp.halfspace_solid_halfspace_leaky(frequency,complex(2*math.pi*frequency/ca[0],ca[1]),nlayers,vp,vs,density,thickness)
    return abs(val)

def free_leaky(phase_velocity,frequency, attenuation, nlayers,vp,vs,density,thickness):
    val=disp.free_solid_halfspace_leaky(frequency,complex(2*math.pi*frequency/phase_velocity,attenuation),nlayers,vp,vs,density,thickness)
    return val

def free_leaky_log_abs(ca,frequency,nlayers,vp,vs,density,thickness):
    val=disp.free_solid_halfspace_leaky(frequency,complex(2*math.pi*frequency/ca[0],ca[1]),nlayers,vp,vs,density,thickness)
    return math.log10(abs(val))

def free_leaky_abs(ca,frequency,nlayers,vp,vs,density,thickness):
    val=disp.free_solid_halfspace_leaky(frequency,complex(2*math.pi*frequency/ca[0],ca[1]),nlayers,vp,vs,density,thickness)
    return abs(val)

def leaky_f_real(func):
    def newfunc(frequency,attenuation,phase_velocity,nlayers,vp,vs,density,thickness):
        val=func(phase_velocity,frequency,attenuation,nlayers,vp,vs,density,thickness)
        return val.real
    return newfunc

def leaky_f_imag(func):
    def newfunc(frequency,attenuation,phase_velocity,nlayers,vp,vs,density,thickness):
        val=func(phase_velocity,frequency,attenuation,nlayers,vp,vs,density,thickness)
        return val.imag
    return newfunc

def leaky_a_real(func):
    def newfunc(attenuation,frequency,phase_velocity,nlayers,vp,vs,density,thickness):
        val=func(phase_velocity,frequency,attenuation,nlayers,vp,vs,density,thickness)
        return val.real
    return newfunc

def leaky_a_imag(func):
    def newfunc(attenuation,frequency,phase_velocity,nlayers,vp,vs,density,thickness):
        val=func(phase_velocity,frequency,attenuation,nlayers,vp,vs,density,thickness)
        return val.imag
    return newfunc

def leaky_roots_real(func,f1,f2,a1,a2,c,n,vp,vs,rho,d):
    func_f_real=leaky_f_real(func)
    func_a_real=leaky_a_real(func)
    real_roots=np.zeros((2,4))

    if func(c,f2,a1,n,vp,vs,rho,d).real * func(c,f1,a1,n,vp,vs,rho,d).real <= 0:
        r = opt.root_scalar(func_f_real,args=(a1,c,n,vp,vs,rho,d),bracket=[f1,f2])
        real_roots[:,3]=[r.root,a1]
    else:
        real_roots = np.delete(real_roots,3,1)

    if func(c,f2,a2,n,vp,vs,rho,d).real * func(c,f2,a1,n,vp,vs,rho,d).real <= 0:
        r = opt.root_scalar(func_a_real,args=(f2,c,n,vp,vs,rho,d),bracket=[a1,a2])
        real_roots[:,2]=[f2,r.root]
    else:
        real_roots = np.delete(real_roots,2,1)

    if func(c,f1,a2,n,vp,vs,rho,d).real * func(c,f2,a2,n,vp,vs,rho,d).real <= 0:
        r = opt.root_scalar(func_f_real,args=(a2,c,n,vp,vs,rho,d),bracket=[f1,f2])
        real_roots[:,1]=[r.root,a2]
    else:
        real_roots = np.delete(real_roots,1,1)

    if func(c,f1,a1,n,vp,vs,rho,d).real * func(c,f1,a2,n,vp,vs,rho,d).real <= 0:
        r = opt.root_scalar(func_a_real,args=(f1,c,n,vp,vs,rho,d),bracket=[a1,a2])
        real_roots[:,0]=[f1,r.root]
    else:
        real_roots = np.delete(real_roots,0,1)

    return real_roots

def leaky_roots_imag(func,f1,f2,a1,a2,c,n,vp,vs,rho,d):
    func_f_imag=leaky_f_imag(func)
    func_a_imag=leaky_a_imag(func)
    imag_roots=np.zeros((2,4))

    if func(c,f2,a1,n,vp,vs,rho,d).imag * func(c,f1,a1,n,vp,vs,rho,d).imag <= 0:
        r = opt.root_scalar(func_f_imag,args=(a1,c,n,vp,vs,rho,d),bracket=[f1,f2])
        imag_roots[:,3]=[r.root,a1]
    else:
        imag_roots = np.delete(imag_roots,3,1)

    if func(c,f2,a2,n,vp,vs,rho,d).imag * func(c,f2,a1,n,vp,vs,rho,d).imag <= 0:
        r = opt.root_scalar(func_a_imag,args=(f2,c,n,vp,vs,rho,d),bracket=[a1,a2])
        imag_roots[:,2]=[f2,r.root]
    else:
        imag_roots = np.delete(imag_roots,2,1)

    if func(c,f1,a2,n,vp,vs,rho,d).imag * func(c,f2,a2,n,vp,vs,rho,d).imag <= 0:
        r = opt.root_scalar(func_f_imag,args=(a2,c,n,vp,vs,rho,d),bracket=[f1,f2])
        imag_roots[:,1]=[r.root,a2]
    else:
        imag_roots = np.delete(imag_roots,1,1)

    if func(c,f1,a1,n,vp,vs,rho,d).imag * func(c,f1,a2,n,vp,vs,rho,d).imag <= 0:
        r = opt.root_scalar(func_a_imag,args=(f1,c,n,vp,vs,rho,d),bracket=[a1,a2])
        imag_roots[:,0]=[f1,r.root]
    else:
        imag_roots = np.delete(imag_roots,0,1)

    return imag_roots

def intersection(real_roots,imag_roots):
    ar = (real_roots[1,1]-real_roots[1,0])/(real_roots[0,1]-real_roots[0,0])
    br = -ar*real_roots[0,0]+real_roots[1,0]
    ai = (imag_roots[1,1]-imag_roots[1,0])/(imag_roots[0,1]-imag_roots[0,0])
    bi = -ai*imag_roots[0,0]+imag_roots[1,0]

    f = - (bi-br)/(ai-ar)
    a = ar * f + br

    return f,a


def leaky_roots(func, phase_velocity, fmin, fmax, df, amin, amax, da, nlayers, vp, vs, density, thickness):
    """
    Compute all complex roots of a characteristic function by frequency/attenuation sweep between fmin and fmax, amin and amax for a given phase velocity
    
    Input:
        func: function - characteristic function
        phase_velocity: float
        fmin,fmax,df: float - minimum, maximum and step for frequency sweep
        amin,amax,da: float - minimum, maximum and step for attenuation sweep (imaginary wavenumber)
        nlayers: integer - number of layers in the system (including half-spaces)
        vp, vs, density, thickness: 1D numpy arrays - elastic properties and thicknesses for the nlayers (any random thickness for half-spaces) 

    Output:
        2D numpy array containing all the frequency/attenuation roots
    """

    roots = np.zeros((2,20))
    nf=round((fmax-fmin)/df)
    na=round((amax-amin)/da)
    nr=0
    c=phase_velocity

    func_f_real=leaky_f_real(func)
    func_f_imag=leaky_f_imag(func)
    func_a_real=leaky_a_real(func)
    func_a_imag=leaky_a_imag(func)

    for i in range(nf-1):
        f1=fmin+i*df
        f2=f1+df
        for j in range(na-1):
            a1=amin+j*da
            a2=a1+da
            real_roots = leaky_roots_real(func,f1,f2,a1,a2,c,nlayers,vp,vs,density,thickness)
            imag_roots = leaky_roots_imag(func,f1,f2,a1,a2,c,nlayers,vp,vs,density,thickness)
            if real_roots.size==4 and imag_roots.size==4:
                f_root,a_root=intersection(real_roots, imag_roots)
                if f_root>=f1 and f_root<=f2 and a_root>=a1 and a_root<=a2:
                    roots[:,nr] = [f_root,a_root]
                    nr += 1

    roots = roots[0:2,0:nr]
    return roots

def leaky_dispersion_curve(func, fca_roots, fmax, df, threshold, nlayers, vp, vs, density, thickness):
    """
    Trace a dispersion curve in the f-c-a domain starting from two given f-c-a roots

    Input:
        func: function - characteristic function
        fca_roots: 3x2 numpy array with 2 f-c-a roots
        fmax, df: float - maximum frequency and step for the dispersion curve
        nlayers: integer - number of layers in the system (including half-spaces)
        vp, vs, density, thickness: 1D numpy arrays - elastic properties and thicknesses for the nlayers (any random thickness for half-spaces) 

    Output:
        2D numpy array containing all the f-c-a roots (the dispersion curve)
    """
    nf = round((fmax-fca_roots[0])/df)
    nf=int(nf)
    fca = np.zeros((3,nf+1))
    fca[:,0] = fca_roots

    bounds=opt.Bounds([0,0],[float('inf'),float('inf')],False)
    options={}
    options['maxiter']=2000
    options['disp']=False
    method='Nelder-Mead'
    options['maxfev']=5000
    #options['xatol']=5000
    options['fatol']=10
    #method='TNC'

    # compute the first 5 roots
    f = fca_roots[0] + df
    res = opt.minimize(func,args=(f,nlayers,vp,vs,density,thickness),x0=fca[1:3,0],method=method,options=options,bounds=bounds)
    if res.success==True and res.x[0]>0 and res.x[1]>0 and (res.fun<threshold or threshold<=0):
        fca[:,1] = [f,res.x[0],res.x[1]]
    else:
        fca = fca[0:3,0:1]
        print("f = %f" %f)
        print(res)
        return fca

    f += df
    res = opt.minimize(func,args=(f,nlayers,vp,vs,density,thickness),x0=2*fca[1:3,1]-fca[1:3,0],method=method,options=options,bounds=bounds)
    if res.success==True and res.x[0]>0 and res.x[1]>0 and (res.fun<threshold or threshold<=0):
        fca[:,2] = [f,res.x[0],res.x[1]]
    else:
        fca = fca[0:3,0:2]
        print("f = %f" %f)
        print(res)
        return fca

    f += df
    res = opt.minimize(func,args=(f,nlayers,vp,vs,density,thickness),x0=2*fca[1:3,2]-fca[1:3,1],method=method,options=options,bounds=bounds)
    if res.success==True and res.x[0]>0 and res.x[1]>0 and (res.fun<threshold or threshold<=0):
        fca[:,3] = [f,res.x[0],res.x[1]]
    else:
        fca = fca[0:3,0:3]
        print("f = %f" %f)
        print(res)
        return fca
    
    f += df
    res = opt.minimize(func,args=(f,nlayers,vp,vs,density,thickness),x0=2*fca[1:3,3]-fca[1:3,2],method=method,options=options,bounds=bounds)
    if res.success==True and res.x[0]>0 and res.x[1]>0 and (res.fun<threshold or threshold<=0):
        fca[:,4] = [f,res.x[0],res.x[1]]
    else:
        fca = fca[0:3,0:4]
        print("f = %f" %f)
        print(res)
        return fca

    f += df
    res = opt.minimize(func,args=(f,nlayers,vp,vs,density,thickness),x0=2*fca[1:3,4]-fca[1:3,3],method=method,options=options,bounds=bounds)
    if res.success==True and res.x[0]>0 and res.x[1]>0 and (res.fun<threshold or threshold<=0):
        fca[:,5] = [f,res.x[0],res.x[1]]
    else:
        fca = fca[0:3,0:5]
        print("f = %f" %f)
        print(res)
        return fca

    # compute the rest of the roots using quadratic extrapolation with alternate points
    jf = 5
    for jf in range(6,nf,1):
        f += df
        ca=fca[1:3,jf-6]-3*fca[1:3,jf-4]+3*fca[1:3,jf-2]
        res = opt.minimize(func,args=(f,nlayers,vp,vs,density,thickness),x0=ca,method=method,options=options,bounds=bounds)
        if res.success==True and res.x[0]>0 and res.x[1]>0 and (res.fun<threshold or threshold<=0):
            fca[:,jf] = [f,res.x[0],res.x[1]]
        else:
            print("f = %f" %f)
            print(res)
            break

    fca = fca[0:3,0:jf]

    return fca
