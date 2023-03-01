import pandas as pd
import numpy as np
import csv
# from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
db = pd.read_csv("Test_Data.csv", dtype = float).to_numpy()

with open("Test_Data.csv",'r') as f:
    reader = csv.reader(f)
    fig = plt.figure()
    ax = fig.add_subplot(111,projection = "3d")
    xcoords = []
    ycoords = []
    zcoords = []
    for i, line in enumerate(reader):
        if float(line[0]) != 0:
            print(i,line)
            xcoords.append(float(line[1]))
            ycoords.append(float(line[2]))
            zcoords.append(float(line[3]))
    ax.scatter(xcoords,ycoords,zcoords)
    plt.show()

            
        
            
