import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

sx = qt.sigmax()
sy = qt.sigmay()
sz = qt.sigmaz()

pulse_dict = {
        'I' : qt.identity(2),
        'X' : -1j * sx,
        'Y' : -1j * sy,
        'Z' : -1j * sz
}

# -- Functions --

# Get bloch vector of state
def get_bloch(state):
    """Returns the Bloch sphere coordinates (theta, phi) of a quantum state.
    
    :param state: The quantum state
    :type state: 2x1 matrix
    """
    theta = 2 * np.arccos(min(abs(state[0,0]), 1))
    phi = np.angle(state[1,0]) - np.angle(state[0,0])

    return theta, phi

# Converts spherical coords to cartesian coords
def spherical_to_cartesian(theta, phi):
    """Returns the cartesian coordinates (x, y, z) corresponding to the given spherical coordinates.
    
    :param theta: The polar angle
    :type theta: float
    :param phi: The azimuthal angle
    :type phi: float
    """
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    return x, y, z

def xy4(time_interval=1):
    """Returns a :class:`qubitvis.Sequence` object with the pulses (gates) X, Y, X, Y.
    
    :param time_interval: The time that elapses between each pulse, defaults to 1
    :type time_interval: int or float, optional
    """
    return Sequence('XYXY', time_interval)

# -- Classes --

class Sequence:
    """This is a class that represents a pulse sequence.
    
    :param pulses: Represents the sequence of pulses - Can be a string consisting of X, Y, Z, and/or I (representing the X, Y, Z, and identity gates correspondingly), or an array of 2x2 matrices representing each gate, defaults to []
    :type pulses: list, optional
    :param time_interval: The time that elapses between each pulse, defaults to 1
    :type time_interval: int or float, optional
    """

    def __init__(self, pulses=[], time_interval=1):
        self.t_int = time_interval
        self.set_pulses(pulses)
    
    def add_pulse(self, pulse):
        """Appends a pulse (gate) to the sequence.
        
        :param pulse: A 2x2 matrix representing the gate, or one of the strings `'X'`, `'Y'`, `'Z'`, or `'I'` (representing the X, Y, Z, and identity gates correspondingly)
        :type pulse: str or 2x2 matrix
        """
        if isinstance(pulse, str):
            self.pulses.append( pulse_dict[pulse] )
        else:
            self.pulses.append( qt.Qobj(pulse) )
    
    def set_pulses(self, pulses):
        """Sets the sequence of pulses (gates).
        
        :param pulses: Represents the sequence of pulses - Can be a string consisting of X, Y, Z, and/or I (representing the X, Y, Z, and identity gates correspondingly), or an array of 2x2 matrices representing each gate
        :type pulses: list
        """
        self.pulses = []
        for pulse in pulses:
            self.add_pulse(pulse)
    
    def get_pulses(self):
        """Returns the sequence of pulses (gates) as an array of :class:`qutip.Qobj`.
        """
        return self.pulses
    
    def U(self, H=[[0,0], [0,0]]):
        """Returns the unitary operator for one iteration of the entire sequence of pulses.
        
        :param H: The Hamiltonian for free evolution, defaults to [[0,0], [0,0]]
        :type H: 2x2 matrix, optional
        """
        H = qt.Qobj(H)
        free_ev = (-1j * self.t_int * H).expm()
        U = qt.identity(2)
        for pulse in reversed(self.pulses):
            U *= pulse * free_ev
        
        return U
    
    def U_t(self, time, H=[[0,0], [0,0]]):
        """Returns the unitary operator for the sequence of pulses at the given time.
        
        :param time: Time (starting from the first pulse)
        :type time: float
        :param H: The Hamiltonian for free evolution, defaults to [[0,0], [0,0]]
        :type H: 2x2 matrix, optional
        """
        H = qt.Qobj(H)
        
        if len(self.pulses) == 0:
            return (-1j * time * H).expm()  
        
        t_sequence = len(self.pulses) * self.t_int
        n = time//t_sequence
        # Correct for floating point error
        if np.allclose(time/t_sequence, np.round(time/t_sequence)):
            n = np.round(time/t_sequence)
        U = self.U(H) ** n
        
        time -= t_sequence * n
        num_pulses = int(time//self.t_int)
        # Correct for floating point error
        if np.allclose(time/self.t_int, np.round(time/self.t_int)):
            num_pulses = int(np.round(time/self.t_int))
        time -= num_pulses * self.t_int
        
        free_ev = (-1j * self.t_int * H).expm()
        U2 = (-1j * time * H).expm()
        for pulse in reversed(self.pulses[:num_pulses]):
            U2 *= pulse * free_ev
        
        return U2 * U

class Qubit:
    """This is a class that represents a qubit.
    
    :param state: The initial state, defaults to [[1],[0]]
    :type state: 2x1 matrix, optional
    :param H: The Hamiltonian, defaults to [[0,0],[0,0]]
    :type H: 2x2 matrix, optional
    :param steps: The default number of time steps calculated, defaults to 500
    :type steps: int, optional
    :param dt: The time that elapses between time steps, defaults to 0.1
    :type dt: float, optional
    :param sequence: A sequence of pulses that the qubit repeatedly executes, defaults to a :class:`qubitvis.Sequence` with no pulses
    :type sequence: :class:`qubitvis.Sequence`, optional
    """
    def __init__(self, state=[[1],[0]], H=[[0,0],[0,0]], steps=500, dt=0.1, sequence=Sequence()):
        self.psi = qt.Qobj(state)
        self.H = qt.Qobj(H)
        self.steps = steps
        self.dt = dt
        self.seq = sequence
        
        self.B = None
    
    def magnetic_field(self, B):
        """Sets the magnetic field acting on the qubit (and updates the Hamiltonian correspondingly).
        
        :param B: Represents the magnetic field - The first element is the X component, the second is the y component, and the third is the z component
        :type B: list
        """
        self.B = B
        self.H = B[0]*sx + B[1]*sy + B[2]*sz
    
    def hamiltonian(self, H=None):
        """Returns and/or sets the Hamiltonian.
        
        :param H: The Hamiltonian - Does not set the Hamiltonian if the value of this parameter is None, defaults to None
        :type H: 2x2 matrix, optional
        """
        if H != None:
            self.B = None
            self.H = qt.Qobj(H)
        return self.H
    
    def set_state(self, state):
        """Sets the qubit's initial state.
        
        :param state: The new initial state
        :type state: 2x1 matrix
        """
        self.psi = qt.Qobj(state)
    
    def state(self, time=0):
        """Returns the qubit's state.
        
        :param time: Returns the state at the time given by this parameter
        :type time: float, optional
        """
        return self.seq.U_t(time, self.H) * self.psi
        
    def sequence(self, sequence=None):
        """Returns and/or sets the sequence of pulses (gates).
        
        :param sequence: The sequence of pulses (gates) - Does not set the sequence if the value of this parameter is None, defaults to None
        :type sequence: :class:`qubitvis.Sequence`, optional
        """
        if isinstance(sequence, Sequence):
            self.seq = sequence
        return self.seq
    
    def show(self, time=0):
        """Shows the qubit as a Bloch sphere.
        
        :param time: Shows the qubit state at the time given by this parameter
        :type time: float, optional
        """
        b = qt.Bloch()
        b.add_states(self.state(time))
        b.show()
    
    # -- Animation --
    
    def _getxyz(self, steps, dt):
        x, y, z = [], [], []
        t = 0

        while t < steps*dt:
            theta, phi = get_bloch(self.state(t))
            x_i, y_i, z_i = spherical_to_cartesian(theta, phi)

            x.append(x_i)
            y.append(y_i)
            z.append(z_i)
            
            t += dt
            
        return x, y, z
    
    def animate(self, steps=None, dt=None):
        """Shows an animation of the qubit as a Bloch sphere. Note that if you are using a Jupyter Notebook, you must include :code:`%matplotlib notebook` at the start of the notebook.
        
        :param steps: The number of time steps, defaults to the value set during the object's initialization
        :type steps: int, optional
        :param dt: The time elapsed between each time step, defaults to the value set during the object's initialization
        :type dt: float, optional
        """
        if steps is None:
            steps = self.steps
        if dt is None:
            dt = self.dt
        
        fig = plt.figure()
        ax = Axes3D(fig,azim=-40,elev=30)
        b = qt.Bloch(axes=ax)
        
        x, y, z = self._getxyz(steps, dt)
        
        def a(i):
            b.clear()
            b.add_vectors([x[i], y[i], z[i]])
            b.add_points([x[:i+1], y[:i+1], z[:i+1]])
            if self.B != None:
                b.add_vectors(self.B)
            b.make_sphere()
        
        return animation.FuncAnimation(fig, a, np.arange(self.steps), init_func=lambda:ax, repeat=False)
    
    def path(self, steps=None, dt=None):
        """Shows the qubit's path and initial state on the Bloch sphere.
        
        :param steps: The number of time steps, defaults to the value set during the object's initialization
        :type steps: int, optional
        :param dt: The time elapsed between each time step, defaults to the value set during the object's initialization
        :type dt: float, optional
        """
        if steps is None:
            steps = self.steps
        if dt is None:
            dt = self.dt
    
        x, y, z = self._getxyz(steps, dt)
        
        b = qt.Bloch()
        b.add_states(self.psi)
        if self.B != None:
            b.add_vectors(self.B)
        b.add_points([x, y, z])
        b.show()