import csv as csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
file = open("Test2_data.csv",'r')
reader = csv.reader(file)
densityFile = open("Test2_dataDensity.csv",'r')
densityReader = csv.reader(densityFile)
xs = []
ys = []
zs = []
ds = []

xcoords = []
ycoords = []
zcoords = []
densities = []
px = []
py = []
pz = []
pxmag = 0
pymag = 0
pzmag = 0

totalP = []
for line in reader:
    
    for i in range(0,len(line) - 3,6):
        pxTemp = float(line[i+3])
        pyTemp = float(line[i+4])
        pzTemp = float(line[i+5])
        pxTemp = pxTemp**2
        pyTemp = pyTemp**2
        pzTemp = pzTemp**2
        totalP.append(pxTemp + pyTemp + pzTemp)
        x = float(line[i])
        y = float(line[i+1])
        z = float(line[i+2])
        xcoords.append(x)
        ycoords.append(y)
        zcoords.append(z)
        
    px.append(float(line[len(line)-3]))
    py.append(float(line[len(line)-2]))
    pz.append(float(line[len(line)-1]))
    pxmag += float(line[len(line)-3])
    print(f"px,py,pz = {px[-1]}, {py[-1]}, {pz[-1]}")
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
    print(f"Drawing frame {i}")
    scatter._offsets3d = (xtemp,ytemp,ztemp)
    # scat.view_init(elev=10.0,azim=45.0)
    

ani = FuncAnimation(fig,updatePlot,repeat=True,frames=[i for i in range(0,len(xs))])
writer = animation.FFMpegWriter(fps = 10)
ani.save("animation.mp4",writer=writer)
pxDifference = (float(px[len(px)-1]) - float(px[0]))
pyDifference = (float(py[len(px)-1]) - float(py[0]))
pzDifference = (float(pz[len(px)-1]) - float(pz[0]))
print(f"px%diff = {pxDifference}, py%diff = {pyDifference}, pz%diff = {pzDifference}")
plt.show()
fig, ax = plt.subplots()
ax.plot(px,label = "$p_x$")
ax.plot(py, label = "$p_y$")
ax.plot(pz, label = "$p_z$")
# ax.plot(totalP,label = "$P^4$",color = "purple")

# plt.gca().set_ylim((min(px),max(pz)))
# plt.gca().set_xlim((0,len(px)))
yticks = np.linspace(float(min(px + py + pz)),float(max(px + py + pz)),10)

# plt.gca().set_yticks(yticks)
ax.legend()
ax.set_xlabel("Iteration(s)")
ax.set_ylabel("Total Momentum (Kg ms$^{-1})$")
# ax2 = ax.twinx()
# ax2.plot(totalP,label = "$P^4$",color = "purple")
# ax2.set_ylabel("Total $P^2$")
# ax2.legend()
plt.savefig("Momentum.pdf",bbox_inches = "tight")
plt.show()
print(px[0])

# Work out density profile
densities = []
density = []
xcoords = []
ycoords = []
zcoords = []
xs = []
ys = []
zs = []
count = 0
densitiesList = []
total_densities = []
for line in densityReader:
    line.pop() 
    count += 1
    for i in range(0,len(line),4):
        d = float(line[i])
        # x = float(line[i+1])
        # y = float(line[i+2])
        # z = float(line[i+3])
        # xcoords.append(x)
        # ycoords.append(y)
        # zcoords.append(z)
        density.append(d)
    print(f"count = {count}")
    densities.append(np.mean(density))
    total_densities.append(np.sum(np.abs(densities)))

    xcoords = []
    ycoords = []
    zcoords = []
    density = []

plt.plot(densities)
plt.xlabel("Iteration(s)")
plt.ylabel("Mean Density $(Kgm^{-3})$")
plt.savefig("density.pdf",bbox_inches = "tight")
plt.show()
plt.plot(total_densities)
plt.xlabel("Iteration(s)")
plt.ylabel("Total Density $(Kgm^{-3})$")
plt.savefig("total_density.pdf",bbox_inches = "tight")
plt.show()

print(f"Gradient of total densities: {(total_densities[-1] - total_densities[0])/len(total_densities)}")
print(f"Gradient of meandensities: {(densities[-1] - densities[0])/len(densities)}")
print(f"Mradient of meandensities: {(densities[-1] - densities[0])/len(densities)}")


# plt.show()



        



