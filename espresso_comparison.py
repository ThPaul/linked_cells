import espressomd
import numpy as np
lj_eps = 1.0
lj_sig = 1.0
lj_cut = 2.5
n_part = 1000000
np.random.seed(seed=42)
system = espressomd.System(box_l=[50,50,50])
system.time_step = 0.01
system.cell_system.skin=0
#system.periodicity = [False]*3

partpos=[]
for i in range(n_part):
    partpos.append(np.random.random(3)*system.box_l)
system.part.add(pos=partpos,type=[0]*n_part)
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")
np.savetxt("partPos",partpos,delimiter=",")
system.integrator.run(0)
print("LJ-parameters:")
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())
allpart=[]
for p in system.part:
    allpart.append(p)
allpart.sort(key=lambda x: x.pos[0]*x.pos[0]+x.pos[1]*x.pos[1]+x.pos[2]*x.pos[2])
forces=[]
for p in allpart:
    #print(p.f)
    forces.append(p.f)
    #print(p.pos)
np.savetxt("forces_esp",forces,delimiter=",")

