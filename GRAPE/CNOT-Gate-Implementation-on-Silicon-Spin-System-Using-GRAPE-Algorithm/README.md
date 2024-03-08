# CNOT Gate Implementation on Silicon Spin System Using GRAPE Algorithm

To implement the CNOT gate on silicon spin systems, we utilized the GRAPE algorithm with QUTIP. 
The system consists of two phosphorus nuclei and one shared electron between them (2P1e)(https://www.nature.com/articles/s41586-021-04292-7). 
The electron's rotation around the nuclei creates hyperfine interaction, which can be used to generate entangling gates such as CNOT. 
However, since we have two nuclei in the system, there are two hyperfine interactions, which are permanent properties of the system 
and cannot be turned off, tuned or manipulated. Therefore, we attempted to implement the CNOT gate using one of these hyperfine interactions, 
while the other hyperfine interaction is a disturbance to the system. To achieve a perfect CNOT gate, we would need to wait long enough 
for the second hyperfine interaction to become negligible, and the GRAPE algorithm was used to examine this effect.

The algorithm was executed with various gate evolution times (1, 10, 100, and 1000 microseconds) while maintaining constant min/max amplitudes.
The theoretical fidelity was found to be optimal at a gate evolution time of 100 microseconds, which can be attributed to the 
interruption of the second hyperfine interaction in the system. In addition, the algorithm was executed for the same gate evolution time, 
but with varying min/max control amplitudes (ranging from 2mT to 40mT). These two scenarios were then tested for different frequencies in the system, 
using both modulated and non-modulated control Hamiltonians. It is important to note that a minimum of 10,000 time slots were required for the algorithm 
to successfully implement the CNOT gate. Any fewer time slots than this will prevent the algorithm from executing the gate properly.

You can find the codes and the results of the codes in this folder


P.S Another solution for implementing a CNOT gate without having to wait for a long time could be to polarize the state of the nuclei. 
For instance, if we are implementing a CNOT gate between the electron and the down state of the nuclei, we can bring the second nucleus to the 
down state as well, so that they cannot interfere with each other. However, this method goes against the universality of computation, 
where we should be able to implement gates perfectly and independently of the states of the system.

This folder contains CNOT gates created with different durations (1, 10, 100, and 1000 microseconds) and the outputs.
We utilized the GRAPE algorithm in two distinct ways, namely modulated and non-modulated GRAPE.

  - For non-modulated GRAPE, instead of taking into account the cosine(sine) term in the control Hamiltonian, 
  our list of controls are in both x and y direction with the proper coefficients 
  i.e., ```H_c = [-gamma_E*tensor(px, p0, p0)- gamma_n*tensor(p0, px, p0)- gamma_n*tensor(p0, p0, px),
       -gamma_E*tensor(py, p0, p0)- gamma_n*tensor(p0, py, p0)- gamma_n*tensor(p0, p0, py)] ```
  - For modulated GRAPE, we take into account the cosine term and modulate the control Hamiltonian for each time slot 
  with an appropriate radial frequency. In this method, we only use the control Hamiltonian for either the x or y direction, not both. 
  The important point is that the GRAPE algorithm receives a list of control Hamiltonians as a parameter, 
  and since we are modulating for each time slot, the control is now a list of lists.
  
  
The algorithm was executed for different frequencies in the system, while keeping the min/max control amplitudes constant throughout the execution. 
The minimum number of time slots is 10000, and the maximum number of time slots is 1000000. 

If you use the code, please cite the software and the paper.
For paper:

@misc{üstün2023singlestep,
      title={Single-Step Parity Check Gate Set for Quantum Error Correction}, 
      author={Gözde Üstün and Andrea Morello and Simon Devitt},
      year={2023},
      eprint={2306.08849},
      archivePrefix={arXiv},
      primaryClass={quant-ph}
}

For software:

Üstün, G., Morello, A., Devitt, S., & Pitchford, A. (2023). CNOT Gate Implementation on Silicon Spin System Using GRAPE Algorithm (Version 1.0.0) [Computer software]. https://github.com/apassenger/CNOT_gate_implementation_with_Grape
or 
@software{Ustun_CNOT_Gate_Implementation_2023,
author = {Üstün, Gözde and Morello, Andrea and Devitt, Simon and Pitchford,  Alexander},
month = apr,
title = {{CNOT Gate Implementation on Silicon Spin System Using GRAPE Algorithm}},
url = {https://github.com/apassenger/CNOT_gate_implementation_with_Grape},
version = {1.0.0},
year = {2023}
}

Citation will not take anything from your work !
