from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

data_E = pd.read_csv("data_Attrattore_E_0.010000.dat", delimiter="\t")
data_RK2 = pd.read_csv("data_Attrattore_RK2_0.010000.dat", delimiter="\t")
data_RK4 = pd.read_csv("data_Attrattore_RK4_0.010000.dat", delimiter="\t")
fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot3D(data_E["x_n"], data_E["y_n"], data_E["z_n"], label='E')
ax.plot3D(data_RK2["x_n"], data_RK2["y_n"], data_RK2["z_n"], label='RK2')
# ax.plot3D(data_RK4["x_n"], data_RK4["y_n"], data_RK4["z_n"], label='RK4')

plt.legend()

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.savefig("LorenzAttractor_3D.png")