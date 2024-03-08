Implementing XX gate with non-modulated control Hamiltonian

The output files include the results of the code, as well as the simulated output fidelities 
for the implemented desired unitary (XX gate) by the single pulse, which is found by the algorithm

You can explore the code by changing
- number of time slots: 
This is a crucial parameter of the algorithm; we need to choose fine enough time slots. 
A very high or very low number of time slots will not allow the algorithm to work.
- gate evolution time (try 1, 10, 100 mucroseconds.)

- the pulse type, as defined in the code:

    A ZERO pulse will give a smoother pulse.
    A SINE type pulse is better for one- and two-qubit gates if we have a good prediction for the possible pulse time.
    An RND type pulse is good for more sophisticated pulse and exploring the space widely. 
    
However, we need to carefully define the experimental parameters, such as the minimum and maximum amplitude, because the RND pulse can be really wild.

Regarding experimental parameters, these include offset, minimum and maximum control amplitude, 
frequencies of the system, driving off and on resonance, phase option, and so on. 
All of these things can change and yield different results.
