# QCSimulationPackage

This is just a toy project to get familiar with packages and GitHub.

In the QCircuit.py there are two classes: **QCircuit** and **QCChain**:
QCircuit updates the state of the system after every addition of a gate.
QChain keeps track of all gates, and one need to call .simulate to know its state.
*QChain is also more flexible*, one can add and delete gates at any place using .add and .delq.
Finally, QChain has a .run module which permits to launch many simulations, and a .plot option.

Demonstration, the 'Hello, world!' of Quantum Computing: Creating Entangled Pairs.
