import sys, os
import matplotlib.pyplot as plt
import numpy as np
from argpar import parse



def postpro(Lift, Drag, time, Re, m, U_dof, run_time, mesh_cells, case):
    args = parse()
    T = args.T; dt = args.dt
    v_deg = args.v_deg; p_deg = args.p_deg
    solver = args.solver
    theta = args.theta
    fig = False
    count = 1

    path = "./experiments/cfd"+case+"/"+str(count)
    while os.path.exists(path):
        count+= 1
        path = "./experiments/cfd"+case+"/"+str(count)

    os.makedirs(path)

    print "Creating report file" + path + "/report.txt"
    name = path+"/report.txt"  # Name of text file coerced with +.txt
    f = open(name, 'w')
    f.write("""cfd%(case)s Turek parameters\n
    \nRe = %(Re)g \nmesh = %(m)s\nDOF = %(U_dof)d\nT = %(T)g\ndt = %(dt)g\nv_deg = %(v_deg)g\np_deg = %(p_deg)g\n"""
    """solver = %(solver)s\ntheta_scheme = %(theta).1f\n""" % vars())
    f.write("""Runtime = %f \n\n""" % run_time)

    f.write("Steady Forces:\nLift Force = %g \n"
            "Drag Force = %g \n" % (Lift[-1], Drag[-1]))

    f.close()

    np.savetxt(path+"/Lift.txt", Lift, delimiter=',')
    np.savetxt(path+"/Drag.txt", Drag, delimiter=',')
    np.savetxt(path+"/time.txt", time, delimiter=',')

    plt.figure(1)
    plt.title("LIFT CFD%s \n Re = %.1f, dofs = %d, cells = %d \n T = %g, dt = %g"
    % (case, Re, U_dof, mesh_cells, T, dt) )
    plt.xlabel("Time Seconds")
    plt.ylabel("Lift force Newton")
    plt.plot(time, Lift, label='dt  %g' % dt)
    plt.legend(loc=4)
    plt.savefig("./experiments/cfd1/"+str(count)+"/lift.png")

    plt.figure(2)
    plt.title("DRAG CFD%s\n Re = %.1f, dofs = %d, cells = %d \n T = %g, dt = %g"
    % (case, Re, U_dof, mesh_cells, T, dt) )
    plt.xlabel("Time Seconds")
    plt.ylabel("Drag force Newton")
    plt.plot(time, Drag, label='dt  %g' % dt)
    plt.legend(loc=4)
    plt.savefig("./experiments/cfd1/"+str(count)+"/drag.png")
    #plt.show()
