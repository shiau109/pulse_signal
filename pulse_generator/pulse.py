# import sys
# sys.path.insert(0, r'../')
# Numpy Series
# Typing
from numpy import ndarray

from numpy import array, linspace, empty

from common_Mathfunc import gaussianFunc
class Waveform():
    """ Data format for the f(x) with linear x input"""
    def __init__( self, x0:float=0, dx:float=1, Y:ndarray=array([]) ):
        self.x0 = x0
        self.dx = dx
        self.Y = Y

    @property
    def Y ( self )->ndarray:
        return self._Y
    @Y.setter
    def Y ( self, value:ndarray ):
        self._Y = value
        self._points = value.shape[-1]

    @property
    def x0 ( self )->ndarray:
        return self._x0
    @x0.setter
    def x0 ( self, value:float ):
        self._x0 = value

    @property
    def dx ( self )->ndarray:
        return self._dx
    @dx.setter
    def dx ( self, value:ndarray ):
        self._dx = value

    @property
    def points ( self )->int:
        """ Array length of Y """
        return self._points

    def get_xAxis ( self ):
        return linspace(self.x0, self.x0+self.dx*self.points,self.points, endpoint=False)


class PulseBuilder():
    """ Store the necessary information for waveform """
    modeList = ["Direct,Mixer"]
    def __init__ ( self ):

        self._duration = None
        self._carrierFrequency = None
        self._carrierPhase = None
        self._envelopeFunc = None
        self._parameters = None

    @property
    def signal ( self )->Waveform:
        """ The output waveform """
        return self._signal
    @signal.setter
    def signal ( self, value ):
        if isinstance(value, Waveform):
            self._signal = value
        elif isinstance(value, ndarray):
            self._signal.Y = value

    @property
    def envelope ( self )->Waveform:
        """ The envelope of the signal."""
        return self._envelope
    @envelope.setter
    def envelope ( self, value ):
        if isinstance(value, Waveform):
            self._envelope = value
        elif isinstance(value, ndarray):
            self._envelope.Y = value


    @property
    def carrierFrequency ( self )->float:
        """ The carrier frequency of the signal, unit depended on dt."""
        return self._carrierFrequency
    @carrierFrequency.setter
    def carrierFrequency ( self, value:float ):
        self._carrierFrequency = value

    @property
    def carrierPhase ( self )->float:
        """ The carrier phase of the signal, unit depended on dt."""
        return self._carrierFrequency
    @carrierPhase.setter
    def carrierPhase ( self, value:float ):
        self._carrierPhase = value

    @property
    def duration ( self )->float:
        """ The duration time of the pulse, unit depended on dt."""
        return self._duration
    @duration.setter
    def duration ( self, value:float ):
        self._duration = value

    @property
    def envelopeFunc ( self ):
        """ The function to form the envelope."""
        return self._envelopeFunc
    @envelopeFunc.setter
    def envelopeFunc ( self, value ):
        self._envelopeFunc = value

    @property
    def parameters ( self )->tuple:
        """ The parameters for the function to form the envelope."""
        return self._parameters
    @parameters.setter
    def parameters ( self, value:tuple ):
        self._parameters = value

    def generate_envelope( self, t0:float, dt:float ):
        """ For a given dt and t0, calculate the waveform"""
        points = int( -(self.duration //-dt) )
        self.envelope = Waveform(t0, dt, empty(points))

        time = self.envelope.get_xAxis()
        self.envelope.Y = self.envelopeFunc( time, *self.parameters )
        


# API
def get_Pulse_gauss ( duration:float, parameters:tuple, carrierFrequency:float=0, carrierPhase:float=0  ):
    """ 
    If the carrier frequency of the pulse is 0 
    p0: Amplitude 
    p1: sigma
    p2: Peak Position
    """
    newPulse = PulseBuilder()
    newPulse.carrierFrequency = carrierFrequency
    newPulse.carrierPhase = carrierPhase
    newPulse.duration = duration
    newPulse.envelopeFunc = gaussianFunc
    newPulse.parameters = parameters
    return newPulse





if __name__ == '__main__':
    import matplotlib.pyplot as plt
    a = get_Pulse_gauss( 10, ( 1, 4, 15 ), 10 )
    a.generate_envelope( 10, 0.5 )
    plt.plot( a.envelope.get_xAxis(), a.envelope.Y, "o")
    plt.show()