
import sys
import numpy as np
import matplotlib.pyplot as plt
import datetime
import dill
from qutip import Qobj, identity, sigmax, sigmaz, sigmay, tensor,qeye
from qutip.qip import hadamard_transform
import qutip.logging_utils as logging
logger = logging.get_logger()
#QuTiP control modules
import qutip.control.pulseoptim as cpo
from qutip import *

example_name = 'CNOT_modulated_100us'
log_level = logging.INFO


# # THE DANGER ZONE - NOTHING CAN CHANGE HERE 



# Spin Hamiltonian
nSpins = 3
[px, py, pz]= jmat(1/2)
p0 = qeye(2)
B0=1.33#tesla
B1=1e-4
gamma_E=-27.97e9*2*np.pi #Hz/tesla
gamma_n=17.23e6*2*np.pi  #GHz/tesla
A_1 = 95e6*2*np.pi
A_2 = 9e6*2*np.pi
H_c = [-gamma_E*tensor(px, p0, p0)- gamma_n*tensor(p0, px, p0)- gamma_n*tensor(p0, p0, px),
       -gamma_E*tensor(py, p0, p0)- gamma_n*tensor(p0, py, p0)- gamma_n*tensor(p0, p0, py)] #control Hamiltonain
print(H_c[0])
print(H_c[1])
#assert 1==2
n_ctrls = len(H_c) #Number of ctrls
Sz = tensor(pz,p0,p0)
Sy = tensor(py,p0,p0)
Sx = tensor(px,p0,p0)
S=[Sx,Sy,Sz]
I_z1 = tensor(p0,pz,p0)
I_z2 = tensor(p0,p0,pz)
I_x1 = tensor(p0,px,p0)
I_x2 = tensor(p0,p0,px)
I_y1 = tensor(p0,py,p0)
I_y2 = tensor(p0,p0,py)

#DRIFT HAMILTONIAN
H_d = -gamma_E*B0*Sz-gamma_n*B0*(I_z1+I_z2)+A_1*(Sx*I_x1+Sy*I_y1+Sz*I_z1) + A_2*(Sx*I_x2+Sy*I_y2+Sz*I_z2)


U_0 = identity(2**nSpins) #Initial unitary: U0

# Target unitary that we want to implement: The important part here is that: we need to specify dims. Although my target unitary 8*8 matrix, because I have 3 spin system, the dims is 2:2:2
#which is 3 qubit system. If you do not specify the dims correctly, you will get error in this part. Also the target unitary has to be QObj type

U_targ=tensor(cnot(),identity(2)) #Target unitary

### Time evolution and resolution parameters################
# Time allowed for the evolution
evo_time = 1e-4 #1 microsecond 
# Number of time slots
n_ts = 100000 # this will make the time resolution in the picoseconds range. In this case: 1e-6/100000=10 picoseconds
              # The minimum number_of_time_slot for this hamiltonain is 10000 and then we can use different number of time slots depending on the gate evolution time. 
              # Less than 10000 number of time slot will not give us the correct solution for this Hamiltonian/example
              # For 100 us, we can increase the n_ts up to 1000000
times = np.linspace(0, evo_time, n_ts) 
#initail_pulse=B1*np.array([ np.sin(gamma_n*B0*times*np.pi/4.0+phi*np.pi/2.0)  for phi in range(len(H_c))]) #this is my initial pulse guess

################ Here, we define the frequencies of the system ##############################################
w=B0*gamma_E  
E0 = H_d.eigenenergies()
print(len(E0))
w2 = E0[1]-E0[0]
w3 = E0[2]-E0[0]
w4 = E0[3]-E0[0]
w5 = E0[4]-E0[0]
w6 = E0[5]-E0[0]
w7 = E0[6]-E0[0]
w8 = E0[7]-E0[0]
###
w10 = E0[1]-E0[0]
w21 = E0[2]-E0[1]
w32 = E0[3]-E0[2]
w43 = E0[4]-E0[3]
###
w41 = E0[4]-E0[1]
w31 = E0[3]-E0[1]
w42 = E0[4]-E0[2]
print(w,w2,w3,w4,w5,w6,w7,w8,w10,w21,w32,w43,w21,w32,w42)
#assert(1==2)
#H_c_coef1 =[]    
#H_c_coef2 =[]      
####################################### Here, we defined the modulated control Hamiltonian###########################

def modulate(t):
    return np.cos(w*t)

def modulate2(t):
    return np.sin(w*t)
    
cntrls =[] #The important point is that cntrls is the list of list. We modulate the control hamiltonian for each time slots and then add them into the list which is called as cntrls.
           #The GRAPE receives control hamiltonian as input. This input type should be a list.
for k in range(n_ts):
    cntrls.append([modulate(times[k])*H_c[0]])   # since we choose the radial frequency in the x direction, we modulate the control 
                                                 # hamiltonian with cosine function. Here, the control Hamiltonian is the first element of H_c 
                                                 # which is the non-modulated control hamiltonian list.  
                                             
initail_pulse=B1*np.array([ np.sin(gamma_E*B0*times+phi*np.pi/2.0)  for phi in range(len(H_c))])   #-phi*np.pi/2.0 # In general, there is no single pulse for impleneting the desired uniatry. 
                                                                                                  # There are several pulses can implement the same desired unitary
                                                                                                  # The decison about which pulse can be used depend on the insturments we have in the lab.
                                                                                                  # For the GRAPE algorithm, the resultant pulse is highly depend on our initial guess
                                                                                                  # This why, we can either give a totally random initial guess or we can carefully think about our
                                                                                                  # initial guess. In this example, we know the solution for implementing
                                                                                                  # CNOT gate on our system and it is a sin wave. 
                                                                                                  # So, we give our initial pulse as a sine wave with a proper coefficient. Because our initial
                                                                                                  # guess is not a random pulse but a good prediction for the algorithm, we then use SINE wave
                                                                                                  # as input for pulse_type. 

#The below, we are just printing our parameters for checking    
#print(len(initail_pulse))

U_targ


U_0


H_d


# # OUT OF DANGER ZONE

# End of parameters for our physical system
# We are now defining the GRAPE parameters****************************************************************

# ****************************************************************

# ***** Define the termination conditions *****
# Fidelity error target
fid_err_targ = 1e-15
# Maximum iterations for the optisation algorithm
max_iter = 5000
# Maximum (elapsed) time allowed in seconds
max_wall_time = 1500000000
# Minimum gradient (sum of gradients squared)
# as this tends to 0 -> local minima has been found
min_grad = 1e-20


# Initial pulse type
# pulse type alternatives: RND|ZERO|LIN|SINE|SQUARE|SAW|TRIANGLE|
# Here, if we are implementing a relatively simple gate such as one qubit gate: X,Y,Z gate or 
# if we are implementing the gate that we know(predict) the pulse shape, we can give SINE wave as a parameter for pulse shape and then we can give a very good guess as initial pulse.
#Again, this also depends on how we generally implement the gates in the lab and what kind of pulse we generally use in the lab.
# if we want to explore more different possible pulses, then we can try RND as pulse type
# if we want smoother pulse, we should give ZERO pulse type

p_type = 'SINE'
# *************************************************************
# File extension for output files

f_ext = "{}_n_ts{}_ptype{}.txt".format(example_name, n_ts, p_type)

# Run the optimisation
print("\n***********************************")
print("Starting pulse optimisation")
result = cpo.optimize_pulse_unitary(H_d, cntrls, U_0, 
                                    U_targ, n_ts, evo_time,
                fid_err_targ=fid_err_targ, min_grad=min_grad,
                max_iter=max_iter, max_wall_time=max_wall_time,
                                    alg='GRAPE',
                amp_ubound=B1*20,amp_lbound=-B1*20, #This is the parameter that depends on our limitations in the lab.
                                                    #If we are designing an electron gate, because of the gyromagnetic ratio is very high for electron, the B1 field can be less than 10 mT
                                                    #Actually, it should not be more than 2-5 mT
                                                    #However, if we are designing the nuclei gates, because the gyromagnetic ratio of nucleus is small, we can go up to a few tens of mili Tesla   
                #If we want to fully explore the space of pulses, we can close this line. IF we are sure for our prediction we can open the line. Our solution will be similar to the our initial guess
                init_pulse_params=initail_pulse,                   
                fid_params={'phase_option':'PSU'},#pulse_scaling=0.001,
                out_file_ext=f_ext, init_pulse_type=p_type,
                log_level=log_level, gen_stats=True)

print("\n***********************************")
print("Optimising complete. Stats follow:")
result.stats.report()
print("\nFinal evolution\n{}\n".format(result.evo_full_final))

print("********* Summary *****************")
print("Initial fidelity error {}".format(result.initial_fid_err))
print("Final fidelity error {}".format(result.fid_err))
print("Final gradient normal {}".format(result.grad_norm_final))
print("Terminated due to {}".format(result.termination_reason))
print("Number of iterations {}".format(result.num_iter))
#print("wall time: ", result.wall_time
print("Completed in {} HH:MM:SS.US".        format(datetime.timedelta(seconds=result.wall_time)))
print("***********************************")
# In[ ]:

dill.dump(result,open('/scratch/tk29/gu2385/grape/result_git.dill', 'wb'))#load wb->rb



