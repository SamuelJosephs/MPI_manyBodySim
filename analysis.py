import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
db = pd.read_csv("output.csv").to_numpy(np.double)
arr = []

for i in db:
    arr.append(i[-1])
energy_difference = ((arr[-1] - arr[0])/arr[0]) * 100
print(f"Energy Difference = {energy_difference}%")
plt.plot(arr)
plt.show()
