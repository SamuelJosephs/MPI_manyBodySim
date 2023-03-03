import csv as csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
file = open("Test2_data.csv",'r')
reader = csv.reader(file)

xs = []
ys = []
zs = []
ds = []

xcoords = []
ycoords = []
zcoords = []
densities = []


for line in reader:
    line.pop()
    
    for i in range(0,len(line),3):
        x = float(line[i])
        y = float(line[i+1])
        z = float(line[i+2])
        xcoords.append(x)
        ycoords.append(y)
        zcoords.append(z)
    xs.append(xcoords)
    ys.append(ycoords)
    zs.append(zcoords)
    xcoords = []
    ycoords = []
    zcoords = []


fig = plt.figure()
scat = fig.add_subplot(projection = "3d")
scatter = scat.scatter(xs[0],ys[0],zs[0])
def updatePlot(i):
    xtemp = xs[i]
    ytemp = ys[i]
    ztemp = zs[i]
    print(f"Frawing frame {i}")
    scatter._offsets3d = (xtemp,ytemp,ztemp)
    

ani = FuncAnimation(fig,updatePlot,repeat=True,frames=[i for i in range(0,len(xs))])
writer = animation.FFMpegWriter(fps = 10)
ani.save("animation.mp4",writer=writer)
plt.show()
# plt.show()



        



