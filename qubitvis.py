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
    theta = 2 * np.arccos(min(abs(state[0,0]), 1))
    phi = np.angle(state[1,0]) - np.angle(state[0,0])

    return theta, phi

# Converts spherical coords to cartesian coords
def spherical_to_cartesian(theta, phi):
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    return x, y, z

def xy4(time_interval=10):
    return Sequence('XYXY', time_interval)

# -- Classes --

class Sequence:
    def __init__(self, pulses=[], time_interval=10):
        self.t_int = time_interval
        self.set_pulses(pulses)
    
    def add_pulse(self, pulse):
        if isinstance(pulse, str):
            self.pulses.append( pulse_dict[pulse] )
        else:
            self.pulses.append( qt.Qobj(pulse) )
    
    def set_pulses(self, pulses):
        self.pulses = []
        for pulse in pulses:
            self.add_pulse(pulse)
    
    def U(self, H=[[0,0], [0,0]]):
        H = qt.Qobj(H)
        free_ev = (-1j * self.t_int * H).expm()
        U = qt.identity(2)
        for pulse in reversed(self.pulses):
            U *= pulse * free_ev
        
        return U
    
    def U_t(self, time, H=[[0,0], [0,0]]):
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
    def __init__(self, state=qt.basis(2,0), H=qt.Qobj([[0,0],[0,0]]), steps=500, dt=1, sequence=Sequence()):
        self.psi = qt.Qobj(state)
        self.H = qt.Qobj(H)
        self.steps = steps
        self.dt = dt
        self.sequence = sequence
        
        self.B = None
    
    def magnetic_field(self, B):
        self.B = B
        self.H = B[0]*sx + B[1]*sy + B[2]*sz
    
    def hamiltonian(self, H=None):
        if H != None:
            self.B = None
            self.H = qt.Qobj(H)
        return self.H
    
    def set_state(self, state):
        self.psi = qt.Qobj(state)
    
    def state(self, time=0):
        return self.sequence.U_t(time, self.H) * self.psi
        #return (-1j * time * self.H).expm() * self.psi
    
    def show(self, time=0):
        b = qt.Bloch()
        b.add_states(self.state(time))
        b.show()
    
    # -- Animation --
    
    def getxyz(self, steps, dt):
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
        if steps is None:
            steps = self.steps
        if dt is None:
            dt = self.dt
        
        fig = plt.figure()
        ax = Axes3D(fig,azim=-40,elev=30)
        b = qt.Bloch(axes=ax)
        
        x, y, z = self.getxyz(steps, dt)
        
        def a(i):
            b.clear()
            b.add_vectors([x[i], y[i], z[i]])
            b.add_points([x[:i+1], y[:i+1], z[:i+1]])
            if self.B != None:
                b.add_vectors(self.B)
            b.make_sphere()
        
        return animation.FuncAnimation(fig, a, np.arange(self.steps), init_func=lambda:ax, repeat=False)
    
    def path(self, steps=None, dt=None):
        if steps is None:
            steps = self.steps
        if dt is None:
            dt = self.dt
    
        x, y, z = self.getxyz(steps, dt)
        
        b = qt.Bloch()
        b.add_states(self.psi)
        if self.B != None:
            b.add_vectors(self.B)
        b.add_points([x, y, z])
        b.show()