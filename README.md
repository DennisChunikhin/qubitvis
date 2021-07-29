# QubitVis
**QubitVis** is a python package for visualizing qubits and the time evolution of qubits.  
<br/>
Note: If you are using a Jupyter Notebook, include `%matplotlib notebook` at the start of your notebook for certain functions such as `animate()` to work properly.
## Methods
### Qubit Methods
| Method | Description |
| ------------ | ------------ |
| animate(steps=None, dt=None) | Shows an animation of the qubit as a Bloch sphere. |
| hamiltonian(H=None) | Returns and/or sets the Hamiltonian. |
| magnetic_field(B) | Sets the magnetic field acting on the qubit (and updates the Hamiltonian correspondingly). |
| path(steps=None, dt=None) | Shows the qubit's path and initial state on the Bloch sphere. |
| sequence(sequence=None) | Returns and/or sets the sequence of pulses (gates). |
| set_state(state) | Sets the qubit's initial state. |
| show(time=0) | Shows the qubit as a Bloch sphere. |
| state(time=0) | Returns the qubit's state. |
### Sequence Methods
| Method | Description |
| ------------ | ------------ |
| U(H=[[0, 0], [0, 0]]) | Returns the unitary operator for one iteration of the entire sequence of pulses. |
| U_t(time, H=[[0, 0], [0, 0]]) | Returns the unitary operator for the sequence of pulses at the given time. |
| add_pulse(pulse) | Appends a pulse (gate) to the sequence. |
| get_pulses() | Returns the sequence of pulses (gates) as an array of `qutip.Qobj` objects. |
| set_pulses(pulses) | Sets the sequence of pulses (gates). |
## Links
[Demo](https://github.com/DennisChunikhin/qubitvis/blob/master/Demo.ipynb)  
[Documentation](https://dennischunikhin.github.io/qubitvis/build/html/index.html)