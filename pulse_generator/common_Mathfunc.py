# Numpy
# Typing
from numpy import ndarray
# Numpy array
from numpy import array, append, zeros, ones, where, linspace, arange
# Numpy common math function
from numpy import exp, sqrt, arctan2, cos, sin, angle, radians, sign, log, ceil
# Numpy constant
from numpy import pi



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

def constFunc (t, *p)->ndarray:
    """
    return constant array
    x: array like, shape (n,) \n
    p[0]: value \n
    """
    return p[0]*ones(len(t))

def linearFunc (t, *p)->ndarray:
    """
    return constant array
    x: array like, shape (n,) \n
    p[0]: slope \n
    p[1]: intersection \n
    """
    return p[0]*t+p[1]


if __name__ == '__main__':
    from numpy import linspace
    x = linspace(0,10,10)
    p = (0,1,2)
    print(gaussianFunc(x,*p))