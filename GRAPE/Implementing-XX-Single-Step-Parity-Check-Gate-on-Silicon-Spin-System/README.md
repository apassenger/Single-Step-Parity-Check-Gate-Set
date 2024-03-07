# Implementing-XX-Single-Step-Parity-Check-Gate-on-Silicon-Spin-System

This repository contains the codes written for the research paper entitled "Single-Step Parity Check Gate for Quantum Error Correction." 
We selected the 2P1e system (two Phosphorus nuclei and one shared electron between them) to implement this gate, as silicon spin systems have the 
inherent capability to perform multi-body measurements due to their readout capabilities. To assess the feasibility of implementing this gate, 
we investigated the experimental limitations that might arise. Moreover, we explored whether this gate can be implemented by adhering to the experimental 
limitations. Most importantly, we aimed to determine whether a single pulse could implement this gate. To achieve these objectives, we utilized the 
Gradient Ascent Pulse Engineering (GRAPE) algorithm in QUTIP. In particular, the QUTIP GRAPE algorithm is wrapped by the L-BFGS-B optimization algorithm, 
which takes into account second-order derivatives and selects the appropriate step size for changes in the control landscape. 
This approach helps the algorithm to converge in fewer iterations. We also utilized this algorithm in a way that allowed us to use a modulated control 
Hamiltonian, which we refer to as modulated GRAPE.

For further information, you can read our research paper.(paper link)

We set the maximum control amplitude which is the magnetic field in our case in the range of mT and we ran the algorithm
for 1, 10, 100 microseconds. For each evolution time, the minimum time resolution was maintained at 10 picoseconds by
allocating the appropriate number of time slots and to be able to run the algorithm for such a small time resolution, we executed it on the 
HPC -High Performance Computing- cluster. As a result, we could have the perfect gate fidelity for 10 microseconds evolution time for the XX gate. 
We then tried to find a more accurate solution for the XX gate and we decided to design these gate as a nuclear gate so that we can use the higher control 
amplitude and we can have a lower evolution time for our XX gate. We set the evolution time to 4 microseconds for the XX gate. 
We set upper and lower bounds to 40 mT for the control field and we used 400000 time slots which correspond to every 10 picoseconds for a pulse. We run
the algorithm for every frequency in the system and the algorithm was able to find XX parity check gates with a single pulse with an accuracy of 0.9999.

For citation:

@misc{üstün2023singlestep,
      title={Single-Step Parity Check Gate Set for Quantum Error Correction}, 
      author={Gözde Üstün and Andrea Morello and Simon Devitt},
      year={2023},
      eprint={2306.08849},
      archivePrefix={arXiv},
      primaryClass={quant-ph}
}

and

@software{Ustun_Implementing-XX-Single-Step-Parity-Check-Gate-on-Silicon-Spin-System_2023,
author = {Üstün, Gözde and Morello, Andrea and Devitt, Simon and Pitchford,  Alexander},
month = apr,
title = {{Implementing-XX-Single-Step-Parity-Check-Gate-on-Silicon-Spin-System}},
url = {https://github.com/apassenger/Implementing-XX-Single-Step-Parity-Check-Gate-on-Silicon-Spin-System},
version = {1.0.0},
year = {2023}

Please cite the sotware and the paper, citation will not take anything from your work!
}
