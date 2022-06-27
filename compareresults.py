import numpy as np

forces_esp=np.loadtxt("forces_esp",delimiter=",")
forces_lc=np.loadtxt("forces_lc",delimiter=",")
for i in range (0, len(forces_esp)-1):
    if((np.allclose(forces_esp[i],forces_lc[i]))==False):
        print(i)
        print(forces_esp[i])
        print(forces_lc[i])
print(np.allclose(forces_esp,forces_lc))

