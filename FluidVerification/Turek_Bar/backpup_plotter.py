import os, sys
import matplotlib.pyplot as plt
import numpy as np

import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description="########################################"
"\nProgram to help data visualization\n######################################## \n \n"
"-  The program prints all available numerical experiment with paramters \n"
"   The user is asked to select experiments of interest with space between choices")
args = parser.parse_args()

def display_runs(count, filter):
    data_list = []
    data_dict = {}
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

                        info.append(i); info.append(words[-1])
                        data_dict[i] = words[-1]

        if len(filt) != 0:
            check = 0
            for i in filt:
                if(str(filt[i]) == str(data_dict[i])):
                    check += 1
            if check == len(filt):
                print '[%s]' % ', '.join(map(str, info))
                data_list.append(data_dict)
            #break
        elif len(filt) == 0:
            print '[%s]' % ', '.join(map(str, info))
            data_list.append(data_dict)
    return data_list
            #print count


def plot(sel, data_list):
    count = 0
    """
    if len(sel) % 2 == 0:
        for i in range(len(sel) / 2 ):
            print "Getting stuff from", sel[count]

            time = np.loadtxt("./experiments/cfd3/"+str(sel[count])+"/time.txt", delimiter=',')
            Lift = np.loadtxt("./experiments/cfd3/"+str(sel[count])+"/Lift.txt", delimiter=',')
            Drag = np.loadtxt("./experiments/cfd3/"+str(sel[count])+"/Drag.txt", delimiter=',')

            plt.figure(i)
            plt.subplot(221)
            plt.title("LIFT CFD3, case = %d, discr = %s\n Re=%s, Dof=%d, dt=%.3f, theta=%.1f" \
            % (sel[count],  data_list[sel[count]]["Discretization"], data_list[sel[count]]["Re"],\
            int(data_list[sel[count]]["DOF"]), float(data_list[sel[count]]["dt"]), float(data_list[sel[count]]["theta_scheme"]) ))
            #% (sel[count],  data_list[sel[count]-1][-1], data_list[sel[count]-1][2],\
            #int(data_list[sel[count]-1][6]), float(data_list[sel[count]-1][10]), float(data_list[sel[count]-1][14]) ))
            plt.xlabel("Time Seconds")
            plt.ylabel("Lift force Newton")
            plt.plot(time, Lift)
            if float(data_list[sel[count]]["dt"]) == 0.01 and float(data_list[sel[count]]["theta_scheme"]) == 1.0:
                print
            else:
                plt.axis([4, 12, -540, 540])

            plt.subplot(223)
            plt.title("Drag CFD3, case = %d"  % (sel[count]))
            plt.xlabel("Time Seconds")
            plt.ylabel("Drag force Newton")
            plt.plot(time, Drag)
            if float(data_list[sel[count]-1]["dt"]) == 0.01 and float(data_list[sel[count]-1]["theta_scheme"]) == 1.0:
                print
            else:
                plt.axis([11, 12, 425, 445])

            count +=1

            print "Getting stuff from", sel[count]

            time = np.loadtxt("./experiments/cfd3/"+str(sel[count])+"/time.txt", delimiter=',')
            Lift = np.loadtxt("./experiments/cfd3/"+str(sel[count])+"/Lift.txt", delimiter=',')
            Drag = np.loadtxt("./experiments/cfd3/"+str(sel[count])+"/Drag.txt", delimiter=',')


            plt.subplot(222)
            plt.title("LIFT CFD3, case = %d, discr = %s\n Re=%s, Dof=%d, dt=%.3f, theta=%.1f" \
            % (sel[count],  data_list[sel[count]-1]["Discretization"], data_list[sel[count]-1]["Re"],\
            int(data_list[sel[count]-1]["DOF"]), float(data_list[sel[count]-1]["dt"]), float(data_list[sel[count]-1]["theta_scheme"]) ))
            plt.xlabel("Time Seconds")
            plt.ylabel("Lift force Newton")
            plt.plot(time, Lift)
            if float(data_list[sel[count]-1]["dt"]) == 0.01 and float(data_list[sel[count]-1]["theta_scheme"]) == 1.0:
                print
            else:
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
            if float(data_list[sel[count]-1]["dt"]) == 0.01 and float(data_list[sel[count]-1]["theta_scheme"]) == 1.0:
                print
            else:
                plt.axis([11, 12, 425, 445])

            #plt.legend(loc=4)
            #plt.savefig("./experiments/cfd3/"+str(count)+"/lift.png"

            count +=1

    """
    if len(sel) % 2 != 0:
        time = np.loadtxt("./experiments/cfd3/"+str(sel[-1])+"/time.txt", delimiter=',')
        Lift = np.loadtxt("./experiments/cfd3/"+str(sel[-1])+"/Lift.txt", delimiter=',')
        Drag = np.loadtxt("./experiments/cfd3/"+str(sel[-1])+"/Drag.txt", delimiter=',')
        plt.figure(3)
        count = -1
        plt.subplot(221)
        plt.title("LIFT CFD3, case = %d, discr = %s\n Re=%s, Dof=%d, dt=%.3f, theta=%.1f" \
        % (sel[count],  data_list[sel[count]-1]["Discretization"], data_list[sel[count]-1]["Re"],\
        int(data_list[sel[count]-1]["DOF"]), float(data_list[sel[count]-1]["dt"]), float(data_list[sel[count]-1]["theta_scheme"]) ))
        plt.xlabel("Time Seconds")
        plt.ylabel("Lift force Newton")
        plt.plot(time, Lift)
        if float(data_list[sel[count]-1]["dt"]) == 0.01 and float(data_list[sel[count]-1]["theta_scheme"]) == 1.0:
            print
        else:
            plt.axis([4, 12, -540, 540])


        plt.subplot(222)
        plt.title("Drag CFD3, case = %d  " % sel[0])
        plt.xlabel("Time Seconds")
        plt.ylabel("Drag force Newton")
        plt.plot(time, Drag)
        if float(data_list[sel[count]-1]["dt"]) == 0.01 and float(data_list[sel[count]-1]["theta_scheme"]) == 1.0:
            print
        else:
            plt.axis([11, 12, 425, 445])


count = 1
get = ["Re", "DOF", "T", "dt", "solver", "theta_scheme", "Discretization", "time", "mesh"]
print "################################################"
print "Getting data from all runs, please wait ...."
print "################################################\n"

while os.path.exists("./experiments/cfd3/"+str(count)):
    count += 1

s = raw_input("Welcome to experiment plotter. \n"
"To see all experiments type: all\n"
"To filter experiments type: filter\n"
"Eks: from an experiment dataset [T, 14, dt, 0.01, theta_scheme, 1.0]\n"
"Filter can be used by: T 14 dt 0.01\n"
"[Re, DOF, T, dt, solver, theta_scheme, Discretization  ]\n "
"-------:")
if s == "filter":
    filt = {}
    s = raw_input("Apply filters: ")
    selected = map(str, s.split())

    for i in range(len(selected) /2):
        filt[str(selected[2*i])] = selected[2*i + 1]
    data_list = display_runs(count, filt)
    s = raw_input("Enter one or several cases, separated by space: ")
    selected = map(int, s.split())
    plot(selected, data_list)

if s == "all":
    filt = {}
    data_list = display_runs(count, filt)
    s = raw_input("Enter one or several cases, separated by space: ")
    selected = map(int, s.split())
    plot(selected, data_list)


plt.show()
