from espressomd import electrostatics
import espressomd
import numpy as np
lj_eps = 1.0
lj_sig = 1.0
lj_cut = 2.5
n_part = 10000
np.random.seed(seed=100)
system = espressomd.System(box_l=[100,100,100])
system.time_step = 0.01
system.cell_system.skin=0
#system.periodicity = [False]*3
coulombC=1

partpos=[]
charge=[]
i=0
while i<n_part:
    par=np.random.random(3)*system.box_l
    md=2*lj_eps;

    for b in partpos:
        dist=par-b
        md=min(md,dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2])
    if md>lj_eps :
        print(i)
        partpos.append(par)
        charge.append(2*(i%2)-1)
        i=i+1
    else:
        print("no")
system.part.add(pos=partpos,type=[0]*n_part,q=charge)
#system.non_bonded_inter[0, 0].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")
p3m = electrostatics.P3M(prefactor=coulombC,accuracy=0.001)
#short range only
#p3m = electrostatics.P3M(prefactor=coulombC,r_cut=2*lj_cut,accuracy=1,alpha=0.666,tune=0,cao=1,mesh=[2,2,2])
#long range only
#p3m = electrostatics.P3M(prefactor=coulombC,r_cut=0,mesh=[20,20,20],alpha=0.26840984934898815,accuracy=1,tune=0,cao=6)
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

