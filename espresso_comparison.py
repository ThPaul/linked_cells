from espressomd import electrostatics
import espressomd
import numpy as np
lj_eps = 1.0
lj_sig = 1.0
lj_cut = 2.5
n_part = 1000
np.random.seed(seed=10)
system = espressomd.System(box_l=[50,50,50])
system.time_step = 0.01
system.cell_system.skin=0
#system.periodicity = [False]*3
coulombC=1

partpos=[]
charge=[]
for i in range(n_part):
    partpos.append(np.random.random(3)*system.box_l)
    charge.append(2*(i%2)-1)
system.part.add(pos=partpos,type=[0]*n_part,q=charge)
#system.non_bonded_inter[0, 0].lennard_jones.set_params(
#   epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")
#p3m = electrostatics.P3M(prefactor=coulombC,r_cut=lj_cut,accuracy=1e-2)
p3m = electrostatics.P3M(prefactor=coulombC,r_cut=0,mesh=[30,30,30],alpha=0.666,accuracy=1,tune=0,cao=6)
system.actors.add(p3m)
np.savetxt("partPos",partpos,delimiter=",")
np.savetxt("charges",charge,delimiter=",")
system.integrator.run(0)
print("LJ-parameters:")
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

print("PM3")
print(p3m.get_params())
allpart=[]
for p in system.part:
    allpart.append(p)
allpart.sort(key=lambda x: x.pos[0]*x.pos[0]+x.pos[1]*x.pos[1]+x.pos[2]*x.pos[2])
forces=[]
charges=[]
possorted=[]
for p in allpart:
    #print(p.f)
    forces.append(p.f)
    possorted.append(p.pos)

    #print(p.pos)
np.savetxt("forces_esp",forces,delimiter=",")
np.savetxt("possorted",possorted,delimiter=",")

