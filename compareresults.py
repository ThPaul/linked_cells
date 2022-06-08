import numpy as np

forces_esp=np.loadtxt("forces_esp",delimiter=",")

forces_lc=np.loadtxt("forces_lc",delimiter=",")
print(np.allclose(forces_esp,forces_lc))

