import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

radJupiter = 778e9

db = pd.read_csv("Test_Data.csv",dtype = float).to_numpy(np.double)

arrTotalE = []
arrTotalGPE = []
arrTotalKE = []

numParticles = int((int(len(db[0])) - int(3)) / int(9))
print(f"numParticles = {numParticles}")
xPositions = np.zeros((numParticles,len(db))) # matrix where the i'th row is the i'th particle and the j'th column is it's position at the j'th timestep
yPositions = np.zeros((numParticles,len(db))) #
zPositions = np.zeros((numParticles,len(db))) #

for n,i in enumerate(db):
    arrTotalE.append(i[-1])
    arrTotalGPE.append(i[-2])
    arrTotalKE.append(i[-3])
    k = 0
    for j in range(numParticles):
        if j == 0:
            k = 0
        else :
            k = 1
        xPositions[j][n] = i[(j*9) + 1]
        yPositions[j][n] = i[(j*9) + 2]
        zPositions[j][n] = i[(j*9) + 3]


print(f"First energy input = {arrTotalE[0]}")
arrTotalE = [i for i in map(lambda x: 100*(arrTotalE[0] - x)/arrTotalE[0],arrTotalE)] # This is dumb but I had already gone through the trouble of using map.
energy_difference = arrTotalE[-1]
print(f"Energy Difference = {energy_difference}%")
plt.plot(arrTotalE)
plt.xlabel("iteration(s)")
plt.ylabel("Energy difference (%)")
plt.savefig("TotalE.pdf",bbox_inches = "tight")
plt.show()


fig = plt.figure(1)
ax = plt.axes(projection = '3d')
for n in range(len(xPositions)):

    ax.plot3D(xPositions[n],yPositions[n],zPositions[n])
    # ax.plot3D(xPositions[1],yPositions[1],zPositions[1])
    # ax.plot3D(xPositions[2],yPositions[2],zPositions[2])
a = 10
ax.set_xlim(-a*radJupiter,a*radJupiter)
ax.set_ylim(-a*radJupiter,a*radJupiter)
ax.set_zlim(-a*radJupiter,a*radJupiter)
plt.show()

print(xPositions[1])