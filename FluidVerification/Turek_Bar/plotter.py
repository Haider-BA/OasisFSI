import os, sys
import matplotlib.pyplot as plt
import numpy as np

def display_runs(count):
    data_list = []
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
        data_list.append(info)
    return data_list
        #print count
def plot(sel, data_list):
    count = 0
    for i in range(len(sel) / 2 ):
        print "Getting stuff from", sel[count]
        print "Count is", count
        print "GETTING VALUES AT PLACE %d in data_list" %sel[count]
        time = np.loadtxt("./experiments/cfd3/"+str(sel[count])+"/time.txt", delimiter=',')
        Lift = np.loadtxt("./experiments/cfd3/"+str(sel[count])+"/Lift.txt", delimiter=',')
        Drag = np.loadtxt("./experiments/cfd3/"+str(sel[count])+"/Drag.txt", delimiter=',')

        plt.figure(i)
        plt.subplot(221)
        plt.title("LIFT CFD3, case = %d, discr = %s\n Re=%s, Dof=%d, dt=%.3f, theta=%.1f" \
        % (sel[count],  data_list[sel[count]-1][-1], data_list[sel[count]-1][2], int(data_list[sel[count]-1][6]), float(data_list[sel[count]-1][10]), float(data_list[sel[count]-1][14]) ))
        plt.xlabel("Time Seconds")
        plt.ylabel("Lift force Newton")
        plt.plot(time, Lift)
        plt.axis([4, 12, -540, 540])

        plt.subplot(223)
        plt.title("Drag CFD3, case = %d"  % (sel[count]))
        plt.xlabel("Time Seconds")
        plt.ylabel("Drag force Newton")
        plt.plot(time, Drag)
        plt.axis([11, 12, 425, 445])

        count +=1

        print "Getting stuff from", sel[count]
        print "Count is", count
        print "GETTING VALUES AT PLACE %d in data_list" %sel[count]

        time = np.loadtxt("./experiments/cfd3/"+str(sel[count])+"/time.txt", delimiter=',')
        Lift = np.loadtxt("./experiments/cfd3/"+str(sel[count])+"/Lift.txt", delimiter=',')
        Drag = np.loadtxt("./experiments/cfd3/"+str(sel[count])+"/Drag.txt", delimiter=',')


        plt.subplot(222)
        plt.title("LIFT CFD3, case = %d, discr = %s\n Re=%s, Dof=%d, dt=%.3f, theta=%.1f" \
        % (sel[count],  data_list[sel[count]-1][-1], data_list[sel[count]-1][2], int(data_list[sel[count]-1][6]), float(data_list[sel[count]-1][10]), float(data_list[sel[count]-1][14]) ))
        plt.xlabel("Time Seconds")
        plt.ylabel("Lift force Newton")
        plt.plot(time, Lift)
        plt.axis([4, 12, -540, 540])
        #plt.legend(loc=4)
        #plt.savefig("./experiments/cfd3/"+str(count)+"/lift.png")

        #plt.title("DRAG CFD3 \n Re = %.1f, dofs = %d, cells = %d \n T = %g, dt = %g"
        #% (Re, U_dof, mesh_cells, T, dt) )
        plt.subplot(224)
        plt.title("Drag CFD3, case = %d  " % sel[count])
        plt.xlabel("Time Seconds")
        plt.ylabel("Drag force Newton\nTEst\nthis")
        plt.plot(time, Drag)
        plt.axis([11, 12, 425, 445])
        #plt.legend(loc=4)
        #plt.savefig("./experiments/cfd3/"+str(count)+"/lift.png"

        count +=1


    if len(sel) % 2 != 0:
        time = np.loadtxt("./experiments/cfd3/"+str(sel[-1])+"/time.txt", delimiter=',')
        Lift = np.loadtxt("./experiments/cfd3/"+str(sel[-1])+"/Lift.txt", delimiter=',')
        Drag = np.loadtxt("./experiments/cfd3/"+str(sel[-1])+"/Drag.txt", delimiter=',')
        plt.figure(3)
        count = -1
        plt.subplot(221)
        plt.title("LIFT CFD3, case = %d, discr = %s\n Re=%s, Dof=%d, dt=%.3f, theta=%.1f" \
        % (sel[count],  data_list[sel[count]-1][-1], data_list[sel[count]-1][2], int(data_list[sel[count]-1][6]), float(data_list[sel[count]-1][10]), float(data_list[sel[count]-1][14]) ))
        plt.xlabel("Time Seconds")
        plt.ylabel("Lift force Newton")
        plt.plot(time, Lift)
        plt.axis([4, 12, -540, 540])

        plt.subplot(222)
        plt.title("Drag CFD3, case = %d  " % sel[0])
        plt.xlabel("Time Seconds")
        plt.ylabel("Drag force Newton")
        plt.plot(time, Drag)
        plt.axis([11, 12, 425, 445])


count = 1
get = ["Re", "DOF", "T", "dt", "solver", "theta_scheme", "Discretization", "time", "mesh",]

while os.path.exists("./experiments/cfd3/"+str(count)):
    count += 1

data_list = display_runs(count)

var = raw_input("Enter one or two cases, separated by space: ")
print "Plotting case", var
selected = []

print "LSIT VAR", var
for i in range(len(var)/2 + 1):
    selected.append(int(var[2*i]))
print selected
plot(selected, data_list)
plt.show()
