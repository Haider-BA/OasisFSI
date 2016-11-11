from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import os
from argpar import parse
from postproses import postpro
from solvers import Newton_manual
parameters['allow_extrapolation']=True

args = parse()
v_deg = args.v_deg
p_deg = args.p_deg
d_deg = args.d_deg
theta = args.theta
discr = args.discr
T = args.T
dt = args.dt
fig = False

def Venant_Kirchhof(d):
    I = Identity(2)
    F = I - grad(d)
    J = det(F)
    E = 0.5*((inv(F.T)*inv(F))-I)
    return inv(F)*(2.*mu_s*E + lamda*tr(E)*I)*inv(F.T)

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

mesh = Mesh("fluid_new.xml")
m = "fluid_new.xml"

# For refining of mesh
if args.refiner != None:
    for i in range(args.refiner):
        mesh = refine(mesh)

for coord in mesh.coordinates():
    if coord[0]==0.6 and (0.199<=coord[1]<=0.2001): # to get the point [0.2,0.6] end of bar
        print coord
        break

# VectorFunctionSpaces
V1 = VectorFunctionSpace(mesh, "CG", v_deg) # Velocity
V2 = VectorFunctionSpace(mesh, "CG", d_deg) # Structure deformation
Q  = FunctionSpace(mesh, "CG", p_deg)       # Fluid Pressure
VVQ = MixedFunctionSpace([V1,V2,Q])

#Dofs and cells
U_dof = VVQ.dim()
mesh_cells = mesh.num_cells()

#Surfaces
Inlet = AutoSubDomain(lambda x: "on_boundary" and near(x[0],0))
Outlet = AutoSubDomain(lambda x: "on_boundary" and (near(x[0],2.5)))
Wall =  AutoSubDomain(lambda x: "on_boundary" and (near(x[1], 0.41) or near(x[1], 0)))
Bar = AutoSubDomain(lambda x: near(x[1], 0.21) or near(x[1], 0.19) or near(x[0], 0.6 ) )
Circle =  AutoSubDomain(lambda x: "on_boundary" and (( (x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.2)*(x[1] - 0.2)  < 0.0505*0.0505 )  ))
Barwall =  AutoSubDomain(lambda x: "on_boundary" and (( (x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.2)*(x[1] - 0.2)  < 0.0505*0.0505 )  and x[1]>=0.19 and x[1]<=0.21 and x[0]>0.2 ))

#Areas
Bar_area = AutoSubDomain(lambda x: (0.19 <= x[1] <= 0.21) and 0.24<= x[0] <= 0.6) # only the "flag" or "bar"

# FacetFunction for surfaces
boundaries = FacetFunction("size_t",mesh)
boundaries.set_all(0)
Wall.mark(boundaries, 2)
Inlet.mark(boundaries, 3)
Outlet.mark(boundaries, 4)
Bar.mark(boundaries, 5)
Circle.mark(boundaries, 6)
Barwall.mark(boundaries, 7)

# Flag boundary, for balance of momentum on interface
#interface = FacetFunction("size_t", mesh)
#interface.set_all(0)
#Bar.mark(interface, 1)

# Full inner geometry(flag + circle) for lift/drag integral
geometry = FacetFunction("size_t", mesh, 0)
Bar.mark(geometry, 1)
Circle.mark(geometry, 1)
Barwall.mark(geometry, 0)

#Are functions, for variational form
#domains = CellFunction("size_t", mesh)
domains = MeshFunction("size_t", mesh, mesh.topology().dim())
domains.set_all(1)
Bar_area.mark(domains, 2) #Overwrites structure domain

dx = Measure("dx",subdomain_data = domains)
ds = Measure("ds", subdomain_data = boundaries)
#ds1 = Measure("ds", subdomain_data = interface) # For interface (interior)
n = FacetNormal(mesh)

# Parameters
Um = 0.2
H = 0.41
L = 2.5
D = 0.1
k = Constant(dt)
t = 0

# List for storing parameters
time = []; dis_x = []; dis_y = []
Lift = []; Drag = []

#Fluid properties
rho_f   = Constant(1000.0)
mu_f    = Constant(1.0)
nu_f = mu_f/rho_f

inlet = Expression(("1.5*Um*x[1]*(H - x[1]) / pow((H/2.0), 2) * (1 - cos(t*pi/2))/2"\
                        ,"0"), t = 0.0, Um = Um, H = H)

#Structure properties
rho_s = 1E3
mu_s = 0.5E6
nu_s = 0.4
E_1 = 1.4E6
lamda = nu_s*2*mu_s/(1-2*nu_s)
g = Constant((0,-2*rho_s))


# velocity conditions
u_inlet   = DirichletBC(VVQ.sub(0), inlet,    boundaries, 3)
u_wall    = DirichletBC(VVQ.sub(0), ((0, 0)), boundaries, 2)
u_circ    = DirichletBC(VVQ.sub(0), ((0, 0)), boundaries, 6) #No slip on geometry in fluid
u_barwall = DirichletBC(VVQ.sub(0), ((0, 0)), boundaries, 7)

# Deformation conditions
d_inlet   = DirichletBC(VVQ.sub(1), ((0, 0)), boundaries, 3)
d_wall    = DirichletBC(VVQ.sub(1), ((0, 0)), boundaries, 2)
d_circ    = DirichletBC(VVQ.sub(1), ((0, 0)), boundaries, 6) #No slip on geometry in fluid
d_barwall = DirichletBC(VVQ.sub(1), ((0, 0)), boundaries, 7)
d_outlet  = DirichletBC(VVQ.sub(1), ((0, 0)), boundaries, 4)

# Pressure Conditions
p_out     = DirichletBC(VVQ.sub(2), 0, boundaries, 4)
p_barwall = DirichletBC(VVQ.sub(2), 0, boundaries, 7)

# Assemble boundary conditions
bcs = [u_inlet, u_wall, u_circ, u_barwall, \
       d_inlet, d_wall, d_circ, d_barwall, d_outlet,\
       p_out, p_barwall]

# Functions
psi, gamma, eta = TestFunctions(VVQ)

udp = Function(VVQ)
u, d, p  = split(udp)

udp0 = Function(VVQ)
u0, d0, p0  = split(udp0)

d_disp = Function(V2)

# Fluid variational form
Fluid_momentum = (rho_f/k)*inner(u - u0, psi)*dx(1) \
                + rho_f*(theta*inner(dot(grad(u), u), psi) + (1 - theta)*inner(dot(grad(u0), u0), psi) )*dx(1) \
                + inner(theta*sigma_f(p, u) + (1 - theta)*sigma_f(p0, u0), eps(psi))*dx(1) \
                - inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds(2) \
                - inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds(3) \
                - inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds(4) \
                - inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds(5) \
                - inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds(6) \

Fluid_continuity = eta*div(u)*dx(1)

# Structure Variational form
I = Identity(2)
F_ = I - grad(d)
J_ = det(F_)
F_1 = I - grad(d0)
J_1 = det(F_1)

Solid_momentum = ( J_*rho_s/k*inner(u - u0, psi) \
                + rho_s*( J_*theta*inner(dot(grad(u), u), psi) + J_1*(1 - theta)*inner(dot(grad(u0), u0), psi) ) \
                + inner(J_*theta*Venant_Kirchhof(d) + (1 - theta)*J_1*Venant_Kirchhof(d0) , grad(psi))  \
                - (theta*J_*inner(g, psi) + (1-theta)*J_1*inner(g, psi) ) ) * dx(2) \
                - inner(J_*theta*Venant_Kirchhof(d)*n + (1 - theta)*J_1*Venant_Kirchhof(d0)*n , psi)*ds(7)

Solid_deformation = dot(d - d0 + k*(theta*dot(grad(d), u) + (1-theta)*dot(grad(d0), u0) ) \
                    - k*(theta*u + (1 -theta)*u0 ), gamma)  * dx(2)

# Mesh velocity function in fluid domain
d_smooth = inner(grad(d),grad(gamma))*dx(1) - inner(grad(d)*n, gamma)*ds(5)

#Conservation of dynamics (CHECK SIGNS!!)
#dynamic = - inner(Venant_Kirchhof(d('-'))*n('-'), psi('-'))*dS(1) - inner(sigma_f(p('+'), u('+'))*n('+') ,psi('+'))*dS(5)
#dynamic = - inner(Venant_Kirchhof(d)*n, gamma)*ds1(1) - inner(sigma_f(p, u)*n ,gamma)*ds1(1)

F = Fluid_momentum + Fluid_continuity \
  + Solid_momentum + Solid_deformation \
  + d_smooth #+ dynamic

#vel_file = File("./velocity/velocity.pvd")
#vel_file << u0



#Reset counters
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

#[bc.apply(udp0.vector()) for bc in bcs]
Re = Um*D/nu_f
print "SOLVING FOR Re = %f" % Re #0.1 Cylinder diameter
tic()
while t <= T:
    time.append(t)

    if t < 2:
        #inlet.t = t;
        inlet.t = 2;
    if t >= 2:
        inlet.t = 2;

    Newton_manual(F, udp, bcs, J, atol, rtol, max_it, lmbda\
                     , udp_res)

    u_, d_, p_  = udp.split(True)
    #vel << u_

    u0, d0, p0  = udp0.split(True)

    d_disp.vector()[:] = d_.vector()[:] - d0.vector()[:]
    ALE.move(mesh, d_disp)
    mesh.bounding_box_tree().build(mesh)

    drag, lift =integrateFluidStress(p_, u_, geometry)
    Drag.append(drag)
    Lift.append(lift)
    if MPI.rank(mpi_comm_world()) == 0:
        print "Time: ",t ," drag: ",drag, "lift: ",lift

    udp0.assign(udp)

    #Reset
    residual   = 1
    rel_res    = residual
    Iter = 0

    #vel_file << u

    dis_x.append(d_(coord)[0])
    dis_y.append(d_(coord)[1])

    t += dt

run_time = toc()

if MPI.rank(mpi_comm_world()) == 0:
    postpro(Lift, Drag, dis_x, dis_y, time, Re, m, U_dof, run_time, mesh_cells, case = 1)
