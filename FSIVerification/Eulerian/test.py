from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description="Implementation of Turek test case CFD1\n"
"For details: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.550.1689&rep=rep1&type=pdf",\
 formatter_class=RawTextHelpFormatter, \
  epilog="############################################################################\n"
  "Example --> python fsi1.py\n"
  "############################################################################")
group = parser.add_argument_group('Parameters')
group.add_argument("-p_deg",  type=int, help="Set degree of pressure                     --> Default=1", default=1)
group.add_argument("-v_deg",  type=int, help="Set degree of velocity                     --> Default=2", default=2)
group.add_argument("-d_deg",  type=int, help="Set degree of velocity                     --> Default=1", default=1)
group.add_argument("-T",  type=float, help="Set end time                                 --> Default=0.1", default=0.1)
group.add_argument("-dt",  type=float, help="Set timestep                                --> Default=0.001", default=0.001)
group.add_argument("-theta",  type=float, help="Explicit, Implicit, Cranc-Nic (0, 1, 0.5)  --> Default=1", default=1)
group.add_argument("-discr",  help="Write out or keep tensor in variational form --> Default=keep", default="keep")
group.add_argument("-r", "--refiner", action="count", help="Mesh-refiner using built-in FEniCS method refine(Mesh)")
group2 = parser.add_argument_group('Solvers')
group2.add_argument("-solver", help="Newton   -- Fenics built-in module (DEFAULT SOLVER) \n"
"Newton2  -- Manuell implementation\n"
"Piccard  -- Manuell implementation\n", default="Newton")

args = parser.parse_args()

v_deg = args.v_deg
p_deg = args.p_deg
d_deg = args.d_deg
solver = args.solver
theta = args.theta
discr = args.discr
T = args.T
dt = args.dt
fig = False

parameters['allow_extrapolation']=True

#mesh = Mesh("von_karman_street_FSI_fluid.xml")
mesh = Mesh("fluid_new.xml")
m = "fluid_new.xml"
#plot(mesh,interactive=True)
if args.refiner != None:
    for i in range(args.refiner):
        mesh = refine(mesh)

for coord in mesh.coordinates():
    if coord[0]==0.6 and (0.199<=coord[1]<=0.2001): # to get the point [0.2,0.6] end of bar
        print coord
        break

V1 = VectorFunctionSpace(mesh, "CG", v_deg) # Velocity
Q  = FunctionSpace(mesh, "CG", p_deg)       # Fluid Pressure

VVQ = MixedFunctionSpace([V1,Q])

#Dofs and cells
U_dof = V1.dim()
mesh_cells = mesh.num_cells()


#NOS = AutoSubDomain(lambda x: "on_boundary" and( near(x[1],0) or near(x[1], 0.41)))
Inlet = AutoSubDomain(lambda x: "on_boundary" and near(x[0],0))
Outlet = AutoSubDomain(lambda x: "on_boundary" and (near(x[0],2.5)))
Wall =  AutoSubDomain(lambda x: "on_boundary" and (near(x[1], 0.41) or near(x[1], 0)))
#Bar = AutoSubDomain(lambda x: "on_boundary" and (near(x[1], 0.21)) or near(x[1], 0.19) or near(x[0], 0.6 ) )
Bar = AutoSubDomain(lambda x: near(x[1], 0.21) or near(x[1], 0.19) or near(x[0], 0.6 ) )
Circle =  AutoSubDomain(lambda x: "on_boundary" and (( (x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.2)*(x[1] - 0.2)  < 0.0505*0.0505 )  ))
Barwall =  AutoSubDomain(lambda x: "on_boundary" and (( (x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.2)*(x[1] - 0.2)  < 0.0505*0.0505 )  and x[1]>=0.19 and x[1]<=0.21 and x[0]>0.2 ))

Allboundaries = DomainBoundary()

boundaries = FacetFunction("size_t",mesh)
boundaries.set_all(0)
Wall.mark(boundaries, 2)
Inlet.mark(boundaries, 3)
Outlet.mark(boundaries, 4)
Bar.mark(boundaries, 5)
Circle.mark(boundaries, 6)
Barwall.mark(boundaries, 7)

# All boundaries execpt flag, for the surface integrals
Neumann = FacetFunction("size_t",mesh, 0)
DomainBoundary().mark(Neumann, 1)
Bar.mark(Neumann, 0)
Barwall.mark(Neumann, 0)
neu = File("neu.pvd")
neu << Neumann

# Flag boundary, for balance of momentum on interface
interface = FacetFunction("size_t", mesh)
interface.set_all(0)
Bar.mark(interface, 1)

inter = File("inter.pvd")
inter << interface

# Full inner geometry(flag + circle) for lift/drag integral
geometry = FacetFunction("size_t", mesh, 0)
Bar.mark(geometry, 1)
Circle.mark(geometry, 1)
Barwall.mark(geometry, 0)

geom = File("geom.pvd")
geom << geometry

# Area functions, to set pressure = 0 in structure
Bar_area = AutoSubDomain(lambda x: (0.19 <= x[1] <= 0.21) and 0.24<= x[0] <= 0.6) # only the "flag" or "bar"
boundary_parts = FacetFunction("size_t", mesh, 0)
Bar_area.mark(boundary_parts, 1)
bararea = File("bararea.pvd")
bararea << boundary_parts

#Are functions, for variational form
#domains = CellFunction("size_t", mesh)
domains = MeshFunction("size_t", mesh, mesh.topology().dim())
domains.set_all(1)
Bar_area.mark(domains, 2) #Overwrites structure domain

dx = Measure("dx",subdomain_data = domains)
ds = Measure("ds", subdomain_data = Neumann)
#ds = Measure("ds", subdomain_data = boundaries)
dS = Measure("dS", subdomain_data = interface) # For interface (interior)
n = FacetNormal(mesh)

# Parameters
Um = 0.2
H = 0.41
L = 2.5
D = 0.1
k = Constant(dt)
t = 0

#Fluid properties
rho_f   = Constant(1000.0)
mu_f    = Constant(1.0)
nu_f = mu_f/rho_f

inlet = Expression(("1.5*Um*x[1]*(H - x[1]) / pow((H/2.0), 2) * (1 - cos(t*pi/2))/2"\
                        ,"0"), t = 0.0, Um = Um, H = H)

#Structure properties
#FSI 1
rho_s = 1E3
mu_s = 0.5E6
nu_s = 0.4
E_1 = 1.4E6
lamda = nu_s*2*mu_s/(1-2*nu_s)
g = Constant((0,-2*rho_s))


# velocity conditions
u_inlet   = DirichletBC(VVQ.sub(0), inlet, boundaries, 3)
u_wall    = DirichletBC(VVQ.sub(0), ((0, 0)), boundaries, 2)
u_geo = DirichletBC(VVQ.sub(0), ((0, 0)), geometry, 1)
u_barwall = DirichletBC(VVQ.sub(0), ((0, 0)), boundaries, 7)

# Pressure Conditions
p_out     = DirichletBC(VVQ.sub(1), 0, boundaries, 4)
p_barwall = DirichletBC(VVQ.sub(1), 0, boundaries, 7)

# Assemble boundary conditions
bcs = [u_inlet, u_wall, u_barwall, u_geo,\
       p_out, p_barwall]

# Functions
psi,eta = TestFunctions(VVQ)

udp = Function(VVQ)
u, p  = split(udp)

udp0 = Function(VVQ)
u0, p0  = split(udp0)


Re = Um*D/nu_f
print "SOLVING FOR Re = %f" % Re #0.1 Cylinder diameter

def integrateFluidStress(p, u, geo):
    ds_g = Measure("ds", subdomain_data = geo) # surface of geometry
    #test3 = File("geo3.pvd")
    #test3 << geometry

    eps   = 0.5*(grad(u) + grad(u).T)
    sig   = -p*Identity(2) + 2.0*mu_f*eps

    traction  = dot(sig, -n)

    forceX  = traction[0]*ds_g(1)
    forceY  = traction[1]*ds_g(1)
    fX      = assemble(forceX)
    fY      = assemble(forceY)

    return fX, fY

def sigma_f(p, u):
  return - p*Identity(2) + mu_f*(grad(u) + grad(u).T)

def eps(v):
    return 0.5*(grad(v).T + grad(v))

# Fluid variational form
Fluid_momentum = (rho_f/k)*inner(u - u0, psi)*dx(1) \
    + rho_f*(theta*inner(dot(grad(u), u), psi) + (1 - theta)*inner(dot(grad(u0), u0), psi) )*dx(1) \
    + inner(theta*sigma_f(p, u) + (1 - theta)*sigma_f(p0, u0), eps(psi))*dx(1) \
    - inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds \
    - inner(sigma_f(p('+'), u('+'))*n('+') ,psi('+'))*dS(1)
    #- inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds(2) \
    #- inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds(3) \
    #- inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds(4) \
    #- inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds(6) \
    #- inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*dS(1)

Fluid_continuity = eta*div(u)*dx(1)


F = Fluid_momentum + Fluid_continuity \

time = []; dis_x = []; dis_y = []
Lift = []; Drag = []

#info(prm,True)  #get full info on the parameters
#list_linear_solver_methods()

vel = File("velocity/vel.pvd")
tic()
d_up = TrialFunction(VVQ)
J = derivative(F, udp, d_up)
udp_res = Function(VVQ)

#Solver parameters
atol, rtol = 1e-7, 1e-7             # abs/rel tolerances
lmbda = 1.0                         # relaxation parameter
residual   = 1                      # residual (To initiate)
rel_res    = residual               # relative residual
max_it    = 100                     # max iterations
Iter = 0                            # Iteration counter

while t <= T:
    print "Time %f" % t
    time.append(t)

    if t < 2:
        inlet.t = t;
    if t >= 2:
        inlet.t = 2;

    while rel_res > rtol and residual > atol and Iter < max_it:
        A = assemble(J, keep_diagonal=True)
        b = assemble(-F)
        A.ident_zeros()

        [bc.apply(A, b, udp.vector()) for bc in bcs]

        solve(A, udp_res.vector(), b, "mumps")

        udp.vector().axpy(1., udp_res.vector())
        [bc.apply(udp.vector()) for bc in bcs]
        rel_res = norm(udp_res, 'l2')
        residual = b.norm('l2')

        if MPI.rank(mpi_comm_world()) == 0:
            print "Newton iteration %d: r (atol) = %.3e (tol = %.3e), r (rel) = %.3e (tol = %.3e) " \
        % (Iter, residual, atol, rel_res, rtol)
        Iter += 1

    u_, p_  = udp.split(True)

    vel << u_

    u0, p0  = udp0.split(True)

    #d_disp.vector()[:] = d_.vector()[:] - d0.vector()[:]
    #ALE.move(mesh, d_disp)
    #mesh.bounding_box_tree().build(mesh)

    drag, lift =integrateFluidStress(p_, u_, geometry)
    Drag.append(drag)
    Lift.append(lift)
    if MPI.rank(mpi_comm_world()) == 0:
        print "Time: ",t ," drag: ",drag, "lift: ",lift


    udp0.assign(udp)

    #Reset counters
    Iter      = 0
    residual   = 1
    rel_res    = residual
    t += dt

run_time = toc()

if MPI.rank(mpi_comm_world()) == 0:
    count = 1
    while os.path.exists("./experiments/fsi1/"+str(count)):
        count+= 1

    os.makedirs("./experiments/fsi1/"+str(count))

    print("Creating report file ./experiments/fsi1/"+str(count)+"/report.txt")
    name = "./experiments/fsi1/"+str(count)+"/report.txt"  # Name of text file coerced with +.txt
    f = open(name, 'w')
    f.write("""FSI1 Turek parameters\n"""
            """Re = %(Re)g \nmesh = %(m)s\nDOF = %(U_dof)d\nT = %(T)g\ndt = %(dt)g\nv_deg = %(v_deg)g\nd_deg%(d_deg)g\np_deg = %(p_deg)g\n"""
            """solver = %(solver)s\ntheta_scheme = %(theta).1f\nDiscretization = %(discr)s\n""" % vars())
    f.write("""Runtime = %f \n\n""" % run_time)

    f.write("Steady Forces:\nLift Force = %g\n"
            "Drag Force = %g\n\n" % (Lift[-1], Drag[-1]))
    f.close()

    np.savetxt("./experiments/fsi1/"+str(count)+"/Lift.txt", Lift, delimiter=',')
    np.savetxt("./experiments/fsi1/"+str(count)+"/Drag.txt", Drag, delimiter=',')
    np.savetxt("./experiments/fsi1/"+str(count)+"/time.txt", time, delimiter=',')

    plt.figure(1)
    plt.title("LIFT \n Re = %.1f, dofs = %d, cells = %d" % (Re, U_dof, mesh_cells))
    plt.xlabel("Time Seconds")
    plt.ylabel("Lift force Newton")
    plt.plot(time, Lift, label='dt  %g' % dt)
    plt.legend(loc=4)
    plt.savefig("./experiments/fsi1/"+str(count)+"/lift.png")

    plt.figure(2)
    plt.title("DRAG \n Re = %.1f, dofs = %d, cells = %d" % (Re, U_dof, mesh_cells))
    plt.xlabel("Time Seconds")
    plt.ylabel("Drag force Newton")
    plt.plot(time, Drag, label='dt  %g' % dt)
    plt.legend(loc=4)
    plt.savefig("./experiments/fsi1/"+str(count)+"/drag.png")
    #plt.show()

    #vel_file << u
    print "DOFS", VVQ.dim()
    print "Cells", mesh.num_vertices()
    print "Discretization theta = %g" % theta
    print "Lift %g" % Lift[-1]
    print "Drag %g" % Drag[-1]
