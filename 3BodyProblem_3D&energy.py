from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

data4 = pd.read_csv("3Corpi_RK4_2.dat", delimiter="\t")
fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot3D(data4["x1"], data4["y1"], data4["z1"], label=r'$m_{1}$')
ax.plot3D(data4["x2"], data4["y2"], data4["z2"], label=r'$m_{2}$')
ax.plot3D(data4["x3"], data4["y3"], data4["z3"], label=r'$m_{3}$')

plt.legend()

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.savefig("3Body_3D.png")

def energy(data, m):
    t = data["t"]
    k = []
    e = []
    
    for i in range(len(t)):
        r12 = np.sqrt((data["x1"][i] - data["x2"][i]) ** 2 +
                      (data["y1"][i] - data["y2"][i]) ** 2 +
                      (data["z1"][i] - data["z2"][i]) ** 2)
        r23 = np.sqrt((data["x3"][i] - data["x2"][i]) ** 2 +
                      (data["y3"][i] - data["y2"][i]) ** 2 +
                      (data["z3"][i] - data["z2"][i]) ** 2)
        r13 = np.sqrt((data["x1"][i] - data["x3"][i]) ** 2 +
                      (data["y1"][i] - data["y3"][i]) ** 2 +
                      (data["z1"][i] - data["z3"][i]) ** 2)

        # Kinetic Energy
        k1 = (1./2.) * m[0] * (data["vx1"][i]**2 + data["vy1"][i]**2 + data["vz1"][i]**2)
        k2 = (1./2.) * m[1] * (data["vx2"][i]**2 + data["vy2"][i]**2 + data["vz2"][i]**2)
        k3 = (1./2.) * m[2] * (data["vx3"][i]**2 + data["vy3"][i]**2 + data["vz3"][i]**2)
        kinetic_energy = k1 + k2 + k3
        k.append(kinetic_energy)
        
        # Potential Energy 
        potential_energy = - (G * m[0] * m[1] / r12) - (G * m[0] * m[2] / r13) - (G * m[1] * m[2] / r23)
        
        # Total Energy
        e.append(kinetic_energy + potential_energy)
        
    return e

m = [1.6, 0.4, 0.4]
G = 1

plt.figure()
plt.plot(data4["t"], energy(data4, m))

plt.grid(True)
plt.xlabel("t")
plt.ylabel("E")
plt.ylim(-5, 5)
plt.savefig("3Body_energy.png")