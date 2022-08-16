# import sys
# sys.path.insert(0, r'../')
# Numpy Series
# Typing
from numpy import ndarray, complex128, issubdtype
# Array
from numpy import array, linspace, empty, append
# Math
from numpy import cos, sin, exp, arctan2, radians, sign
# const
from numpy import pi

from pulse_generator.common_Mathfunc import gaussianFunc, DRAGFunc
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
    def points( self )->int:
        """ Array length of Y """
        return self._points

    def append( self, appended ):
        if appended.dx == self.dx:
            self.Y = append(self.Y, appended.Y)
        else:
            raise ValueError("dx are different")

    def get_xAxis ( self ):
        return linspace(self.x0, self.x0+self.dx*self.points,self.points, endpoint=False)


class Pulse():
    """ Store the necessary information for waveform """
    def __init__ ( self ):

        self._duration = None
        self._carrierFrequency = None
        self._carrierPhase = None
        self._envelopeFunc = None
        self._parameters = None



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
        return self._carrierPhase
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

    def generate_envelope( self, t0:float, dt:float )->Waveform:
        """ For a given dt and t0, calculate the envelop waveform"""
        points = int( -(self.duration //-dt) )
        envelope = Waveform(t0, dt, empty(points))

        time = envelope.get_xAxis()
        envelope.Y = self.envelopeFunc( time, *self.parameters )
        
        return envelope

    def generate_signal( self, t0:float, dt:float )->Waveform:
        """ For a given dt and t0, calculate the signal waveform"""

        envelope = self.generate_envelope( t0, dt )
        signal = Waveform(envelope.x0, envelope.dx, empty(envelope.Y.shape[-1]))
        time = envelope.get_xAxis()
        if issubdtype(envelope.Y.dtype,complex):
            phase_I = 2.*pi*self.carrierFrequency*time +self.carrierPhase
            print(time[0])
            phase_Q = phase_I +pi/2
            LO_I = cos( phase_I )
            LO_Q = cos( phase_Q )
            signal.Y = envelope.Y.real*LO_I +envelope.Y.imag*LO_Q
        else:
            signal.Y = envelope.Y*cos( 2.*pi*self.carrierFrequency*time +self.carrierPhase)

        return signal

    def generate_IQSignal( self, t0:float, dt:float, IFFreq:float, IQMixer:tuple=(1,90,0,0) ):
        """
        For a given dt and t0, calculate the output for IQMixer 
        to form the target signal waveform.\n

        IFFreq: The Intermediate frequency of I/Q ( Unit in dt ) \n
        IQMixer: The parametrs for calibrate IQmixer\n
            p1: I/Q Amplitude Balance ( dimensionless ratio )\n
            p2: Phase Balance ( unit in angle )\n
            p3: I offset\n
            p4: Q offset\n
        The LO frequency should be RF+IF (RF is carrier frequency)
        """
        envelope = self.generate_envelope( t0, dt )
        time = envelope.get_xAxis()

        Signal_I = Waveform(envelope.x0, envelope.dx, empty(envelope.Y.shape[-1]))
        Signal_Q = Waveform(envelope.x0, envelope.dx, empty(envelope.Y.shape[-1]))
        
        ampBalance = IQMixer[0]
        phaseBalance = IQMixer[1]
        offsetI = IQMixer[2]
        offsetQ = IQMixer[3]
        LOShiftSign = sign(sin(radians(phaseBalance)))

        envelopeIQ = abs( envelope.Y )
        envelopeI = envelopeIQ /cos(radians(abs(phaseBalance)-90))
        envelopeQ = envelopeI /ampBalance
        phi_envelope = arctan2( envelope.Y.imag, envelope.Y.real )
        phi_Q = phi_envelope -LOShiftSign*pi/2 +self.carrierPhase
        phi_I = phi_Q -radians(phaseBalance) +pi
        Signal_I.Y = envelopeI *cos( 2. *pi *IFFreq *time +phi_I) -offsetI
        Signal_Q.Y = envelopeQ *cos( 2. *pi *IFFreq *time +phi_Q) -offsetQ

        return Signal_I, Signal_Q

def convert_envtoIQ( envelope:ndarray, IFFreq:float, IQMixer:tuple=(1,90,0,0) ):
    """
    There is no carrier phase parameter 
    """
    ampBalance = IQMixer[0]
    phaseBalance = IQMixer[1]
    offsetI = IQMixer[2]
    offsetQ = IQMixer[3]
    time = array(range(len(envelope)))

    LOShiftSign = sign(sin(radians(phaseBalance)))
    
    envelopeIQ = abs( envelope )
    envelopeI = envelopeIQ /cos(radians(abs(phaseBalance)-90))
    envelopeQ = envelopeI /ampBalance
    phi_envelope = arctan2( envelope.imag, envelope.real )
    phi_Q = phi_envelope -LOShiftSign*pi/2
    phi_I = phi_Q -radians(phaseBalance) +pi
    Signal_I = envelopeI *cos( 2. *pi *IFFreq *time +phi_I) -offsetI
    Signal_Q = envelopeQ *cos( 2. *pi *IFFreq *time +phi_Q) -offsetQ

    return Signal_I, Signal_Q


def simulate_IQMixer ( I:Waveform, Q:Waveform, LOFreq:float, IQMixer:tuple=(1,90,0,0)):
    time = I.get_xAxis()
    Signal_RF = Waveform(I.x0,I.dx,empty(time.shape[-1]))
    ampBalance = IQMixer[0]
    phaseBalance = IQMixer[1]
    offsetI = IQMixer[2]
    offsetQ = IQMixer[3]
    mixed_I = (I.Y+offsetI)*cos( 2.*pi*LOFreq*time )
    mixed_Q = (Q.Y+offsetQ)*ampBalance*cos( 2.*pi*LOFreq*time +radians(phaseBalance)) 
    Signal_RF.Y = mixed_I+mixed_Q
    return Signal_RF

# API
def get_Pulse_gauss ( duration:float, parameters:tuple, carrierFrequency:float=0, carrierPhase:float=0  ):
    """ 
    Get a Pulse Object
    p0: Amplitude 
    p1: sigma
    p2: Peak Position
    """
    newPulse = Pulse()
    newPulse.carrierFrequency = carrierFrequency
    newPulse.carrierPhase = carrierPhase
    newPulse.duration = duration
    newPulse.envelopeFunc = gaussianFunc
    newPulse.parameters = parameters

    return newPulse


def get_Pulse_DRAG ( duration:float, parameters:tuple, carrierFrequency:float=0, carrierPhase:float=0  ):
    """ 
    Get a Pulse Object \n
    p0: Amplitude \n
    p1: sigma \n
    p2: Peak Position \n
    p3: derivative Gaussian amplitude ratio \n
    """
    newPulse = Pulse()
    newPulse.carrierFrequency = carrierFrequency
    newPulse.carrierPhase = carrierPhase
    newPulse.duration = duration
    newPulse.envelopeFunc = DRAGFunc
    newPulse.parameters = parameters

    return newPulse


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    IFFreq = -0.089
    RFFreq = 5
    LOFreq = RFFreq-IFFreq
    MixerSetting = (1, -90, 0, 0) #IQMixer: tuple = (1, 90, 0, 0)
    testPunse = get_Pulse_DRAG( 30, ( 1, 7.5, 0, 10 ), RFFreq, pi/2 )
    WF_envelope = testPunse.generate_envelope( -15, 1 )

    WF_signal = testPunse.generate_signal( -15, 0.001 )

    WF_signal_I, WF_signal_Q = testPunse.generate_IQSignal( -15, 0.001, IFFreq, IQMixer=MixerSetting)

    WF_signal_RF = simulate_IQMixer( WF_signal_I, WF_signal_Q, LOFreq, IQMixer=MixerSetting)

    # Plot setting
    fig, ax = plt.subplots(3,1,sharex=True)

    # Compare signal and envelope
    ax[0].plot( WF_envelope.get_xAxis(), abs(WF_envelope.Y), label="ABS envelope" )

    ax[0].plot( WF_envelope.get_xAxis(), WF_envelope.Y.real, label="I envelope" )
    ax[0].plot( WF_envelope.get_xAxis(), WF_envelope.Y.imag, label="Q envelope" )

    ax[0].plot( WF_signal.get_xAxis(), WF_signal.Y, label="Target Signal" )
    ax[0].legend()

    # Compare IQ signal and envelope
    ax[1].plot( WF_envelope.get_xAxis(), abs(WF_envelope.Y), label="envelope" )
    ax[1].plot( WF_signal_I.get_xAxis(), WF_signal_I.Y, label="I" )
    ax[1].plot( WF_signal_Q.get_xAxis(), WF_signal_Q.Y, label="Q" )
    ax[1].legend()

    # Compare Desire sigmal and remixed signal by IQ signal
    ax[2].plot( WF_signal_RF.get_xAxis(), WF_signal_RF.Y, label="RF Signal" )
    ax[2].plot( WF_signal.get_xAxis(), WF_signal.Y, label="Target Signal" )
    ax[2].legend()

    plt.show()

