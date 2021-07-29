from setuptools import setup

setup(
    name='QubitVis',
    url='https://github.com/DennisChunikhin/qubitvis',
    author='Dennis Chunikhin',
    author_email='dennis.chunikhin@gmail.com',
    packages=['qubitvis'],
    install_requires=['numpy', 'matplotlib', 'qutip'],
    version='0.0.1',
    license='MIT',
    description='A python package for visualizing qubits and the time evolution of qubits.',
    long_description='A python package for visualizing qubits and the time evolution of qubits. The package can also visualize the effects of pulse sequences on the time evolution of qubits.'
)