# Numpy
# Typing
from numpy import ndarray
# Numpy array
from numpy import array, append, zeros, ones, where, linspace, arange, searchsorted
# Numpy common math function
from numpy import exp, sqrt, arctan2, cos, sin, angle, radians, sign, log, ceil
# Numpy constant
from numpy import pi, logical_and



def gaussianFunc (x, *p)->ndarray:
    """
    x: array like, shape (n,)\n
    p: parameters\n
        p[0]: amp\n
        p[1]: sigma\n
        p[2]: peak position\n
    """

    return p[0] *exp( -( (x-p[2]) /p[1] )**2 /2)
def derivativeGaussianFunc (x, *p)->ndarray:
    """
    return derivative Gaussian
    x: array like, shape (n,) \n
    p: parameters \n
        p[0]: amp \n
        p[1]: sigma \n
        p[2]: peak position \n 
    """

    if p[1] != 0. :
        return -p[0] / p[1]**2 *(x-p[2]) *exp( -( (x-p[2]) /p[1] )**2 /2)
    else :
        return zeros(len(x))

def constFunc (x, *p)->ndarray:
    """
    return constant array
    x: array like, shape (n,) \n
    p[0]: value \n
    """
    return p[0]*ones(len(x))

def rectPulseFunc (x, *p)->ndarray:
    """
    return constant array
    x: array like, shape (n,) \n
    p[0]: amp \n
    p[1]: width \n
    p[2]: start \n
    """
    condition = logical_and(abs(x)>=p[2], abs(x)<=(p[1]+p[2]))
    #print(condition)
    return where(condition, p[0], 0)

def GERP (x, *p)->ndarray:
    """
    return Gaussian Edge Rectangular Pulse array
    x: array like, shape (n,) \n
    p[0]: amp \n
    p[1]: width \n
    p[2]: start \n
    p[3]: edge width \n
    p[4]: edge sigma \n
    """
    amp = p[0]
    width = p[1]
    edgeLength = p[3]
    flat_start = p[2]
    peakLength = edgeLength*2
    flatLength = x -peakLength
    peakSigma = peakLength

    startPos = edgeLength
    endPos = startPos +flatLength

    startEdge = [ p[0], p[4], startPos ]
    gaussUp = where( relativeTime<startPos, gaussianFunc(x, startEdge),0. )
    endEdge = [ flatHieght, peakSigma, endPos ]
    gaussDn = where( relativeTime>endPos, gaussianFunc(relativeTime, endEdge),0. )
    step = where( (relativeTime>=startPos) & (relativeTime<=endPos), constFunc(relativeTime, [flatHieght]),0. )
    return gaussUp +step +gaussDn
def linearFunc (t, *p)->ndarray:
    """
    return constant array
    x: array like, shape (n,) \n
    p[0]: slope \n
    p[1]: intersection \n
    """
    return p[0]*t+p[1]

def DRAGFunc ( t, *p )->ndarray:
    """
    return gaussian +1j*derivative Gaussian\n
    x: array like, shape (n,), the element is complex number \n
    p[0]: amp \n
    p[1]: sigma \n
    p[2]: peak position \n
    p[3]: derivative Gaussian amplitude ratio \n
    """
    gaussParas = (p[0],p[1],p[2])
    return gaussianFunc( t, *gaussParas )+ 1j*p[3]*derivativeGaussianFunc( t, *gaussParas )

if __name__ == '__main__':
    from numpy import linspace
    import matplotlib.pyplot as plt

    x = linspace(0,9,100)
    p = (1,3,0)
    plt.plot(x,derivativeGaussianFunc(x,*p))
    plt.show()