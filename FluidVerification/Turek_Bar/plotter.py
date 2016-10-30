import os
import matplotlib.pyplot as plt
import numpy as np
count = 1
get = ["Re", "DOF", "T", "dt", "solver", "theta_scheme", "Discretization", "time", "mesh",]

while os.path.exists("./experiments/cfd3/"+str(count)):
    count += 1

for i in range(1,count):
    info = []
    info.append("Case %d" % i)
    counter = 0
    with open("./experiments/cfd3/"+str(i)+"/report.txt", 'r') as f:
        data = f.readlines()

        for line in data:
            words = line.split()
            for i in words:
                if any(i in s for s in get):
                    #print i, words[-1]
                    info.append(i); info.append(words[-1])

            #print words
    print '[%s]' % ', '.join(map(str, info))
    #print count

selected = [4]
def plot(selected):
    count = 1
    for i in selected:
        time = np.loadtxt("./experiments/cfd3/"+str(i)+"/time.txt", delimiter=',')
        Lift = np.loadtxt("./experiments/cfd3/"+str(i)+"/Lift.txt", delimiter=',')
        Drag = np.loadtxt("./experiments/cfd3/"+str(i)+"/Drag.txt", delimiter=',')
        plt.figure(count)
        count+= 1
        #plt.title("LIFT CFD3 \n Re = %.1f, dofs = %d, cells = %d \n T = %g, dt = %g"
        #% (Re, U_dof, mesh_cells, T, dt) )
        plt.xlabel("Time Seconds")
        plt.ylabel("Lift force Newton")
        plt.plot(time, Lift)
        plt.axis([4, 12, -500, 500])
        #plt.legend(loc=4)
        #plt.savefig("./experiments/cfd3/"+str(count)+"/lift.png")

        #plt.title("DRAG CFD3 \n Re = %.1f, dofs = %d, cells = %d \n T = %g, dt = %g"
        #% (Re, U_dof, mesh_cells, T, dt) )
        plt.figure(count)
        plt.xlabel("Time Seconds")
        plt.ylabel("Drag force Newton")
        plt.plot(time, Drag)
        plt.axis([11, 12, 425, 450])
        #plt.legend(loc=4)
        #plt.savefig("./experiments/cfd3/"+str(count)+"/lift.png")
        count += 1
plot(selected)
plt.show()
