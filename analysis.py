import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
db = pd.read_csv("output.csv").to_numpy(np.double)

arrTotalE = []
arrTotalGPE = []
arrTotalKE = []


for i in db:
    arrTotalE.append(i[-1])
    arrTotalGPE.append(i[-2])
    arrTotalKE.append(i[-3])
print(f"First energy input = {arrTotalE[0]}")
arrTotalE = [i for i in map(lambda x: 100*(arrTotalE[0] - x)/arrTotalE[0],arrTotalE)] # This is dumb but I had already gone through the trouble of using map.
energy_difference = arrTotalE[-1]
print(f"Energy Difference = {energy_difference}%")
plt.plot(arrTotalE)
plt.xlabel("iteration(s)")
plt.ylabel("Energy difference (%)")
plt.savefig("TotalE.pdf",bbox_inches = "tight")
plt.show()
