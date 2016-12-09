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

#Bar_area_inner_facet = AutoSubDomain(lambda x: (0.19 < x[1] < 0.21) and 0.24< x[0] < 0.6) # only the "flag" or "bar"
#test = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
#test.set_all(0)
#Bar_area_inner_facet.mark(test, 1)
#test_file = File("test.pvd")
#test_file << test

# FacetFunction for surfaces
boundaries = FacetFunction("size_t",mesh)
boundaries.set_all(0)
Wall.mark(boundaries, 2)
Inlet.mark(boundaries, 3)
Outlet.mark(boundaries, 4)
Bar.mark(boundaries, 5)
Circle.mark(boundaries, 6)
Barwall.mark(boundaries, 7)
#bou = File("bound.pvd")
#bou << boundaries

# Full inner geometry(flag + circle) for lift/drag integral
geometry = FacetFunction("size_t", mesh, 0)
Bar.mark(geometry, 1)
Circle.mark(geometry, 1)
Barwall.mark(geometry, 0)

#Are functions, for variational form
domains = CellFunction("size_t", mesh)
#domains = MeshFunction("size_t", mesh, mesh.topology().dim())
domains.set_all(1)
Bar_area.mark(domains, 2) #Overwrites to structure domain

dx = Measure("dx", subdomain_data = domains)
ds = Measure("ds", subdomain_data = boundaries)
dS = Measure("ds", subdomain_data = boundaries) # For interface (interior)
n = FacetNormal(mesh)

# List for storing parameters
time = []; dis_x = []; dis_y = []
Lift = []; Drag = []

#Fluid properties
rho_f   = Constant(1000.0)
mu_f    = Constant(1.0)
nu_f = mu_f/rho_f

# Parameters
Um = 0.2
H = 0.41
L = 2.5
D = 0.1
k = Constant(dt)
t = dt

#Structure properties
rho_s = 1.0E3
mu_s = 0.5E6
nu_s = 0.4
E_1 = 1.4E6
lamda = nu_s*2*mu_s/(1 - 2*nu_s)

# velocity conditions
inlet = Expression(("1.5*Um*x[1]*(H - x[1]) / pow((H/2.0), 2) * (1 - cos(t*pi/2))/2"\
                        ,"0"), t = 0.0, Um = Um, H = H)

u_inlet   = DirichletBC(VVQ.sub(0), inlet,    boundaries, 3)
u_wall    = DirichletBC(VVQ.sub(0), ((0, 0)), boundaries, 2)
u_circ    = DirichletBC(VVQ.sub(0), ((0, 0)), boundaries, 6)
u_barwall = DirichletBC(VVQ.sub(0), ((0, 0)), boundaries, 7)

# Deformation conditions
d_barwall = DirichletBC(VVQ.sub(1), ((0, 0)), boundaries, 7)

# Pressure Conditions
p_out     = DirichletBC(VVQ.sub(2), 0, boundaries, 4)
#p_inner   = DirichletBC(VVQ.sub(2), 0, test, 1)

# Assemble boundary conditions
bcs = [u_inlet, u_wall, u_circ, u_barwall, \
       d_barwall, \
       p_out]

# Functions
psi, gamma, eta = TestFunctions(VVQ)

udp = Function(VVQ)
u, d, p  = split(udp)

udp0 = Function(VVQ)
u0, d0, p0  = split(udp0)

# Fluid variational form
Fluid_momentum = (rho_f/k)*inner(u - u0, psi)*dx(1) \
                + rho_f*(theta*inner(dot(grad(u), u), psi) + (1 - theta)*inner(dot(grad(u0), u0), psi) )*dx(1) \
                + inner(theta*sigma_f(p, u)   + (1 - theta)*sigma_f(p0, u0), eps(psi))*dx(1) \
                - inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds(2) \
                - inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds(3) \
                - inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds(4) \
                - inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds(6)


Fluid_continuity = eta*div(u)*dx(1)

# Structure Variational form
I = Identity(2)
F_ = I - grad(d)
J_ = det(F_)
F_1 = I - grad(d0)
J_1 = det(F_1)

Solid_momentum = ( J_*rho_s/k*inner(u - u0, psi) \
                + rho_s*( J_*theta*inner(dot(grad(u), u), psi) + J_1*(1 - theta)*inner(dot(grad(u0), u0), psi) ) \
                + inner(J_*theta*Venant_Kirchhof(d) + (1 - theta)*J_1*Venant_Kirchhof(d0) , grad(psi))) * dx(2) \

Solid_deformation = inner(d - d0 + k*(theta*dot(grad(d), u) + (1-theta)*dot(grad(d0), u0) ) \
                    - k*(theta*u + (1 -theta)*u0 ), gamma)  * dx(2)


F = Fluid_momentum + Fluid_continuity \
  + Solid_momentum + Solid_deformation

#Reset counters
d_up = TrialFunction(VVQ)
J = derivative(F, udp, d_up)
udp_res = Function(VVQ)

W = VectorFunctionSpace(mesh, 'CG', 1)
w = TrialFunction(W)
chi = TestFunction(W)
w_ = Function(W)
w0 = Function(W)

# mesh deformation
w_inlet   = DirichletBC(W, ((0, 0)), boundaries, 3)
w_wall    = DirichletBC(W, ((0, 0)), boundaries, 2)
w_out     = DirichletBC(W, ((0, 0)), boundaries, 4)
w_circ    = DirichletBC(W, ((0, 0)), boundaries, 6)
w_barwall = DirichletBC(W, ((0, 0)), boundaries, 7)

bcs_w = [w_inlet, w_wall, w_out, w_circ, w_barwall]

F_W = inner(w - (d - d0), chi)*dx(2) + 1./k*inner(w - w0, chi)*dx(1)+ inner(grad(w), grad(chi))*dx(1) #+ 1./k*inner(w - w0, chi)*dx(1)

#Solver parameters
atol, rtol = 1e-10, 1e-10             # abs/rel tolerances
lmbda = 1.0                         # relaxation parameter
residual   = 1                      # residual (To initiate)
rel_res    = residual               # relative residual
max_it    = 15                      # max iterations
Iter = 0                            # Iteration counter

vel_file = File("./velocity/vel.pvd")
def_file = File("./deformation/def.pvd")

#[bc.apply(udp0.vector()) for bc in bcs]
u0, d0, p0  = udp0.split(True)
#u0.rename("u", "velocity")
#d0.rename("d", "deformation")
#vel_file << u0
#def_file << d0

Re = Um*D/nu_f
print "SOLVING FOR Re = %f" % Re #0.1 Cylinder diameter
tic()
while t <= T:
    time.append(t)

    if t < 2:
        inlet.t = t;
    if t >= 2:
        inlet.t = 2;

    Newton_manual(F, udp, bcs, J, atol, rtol, max_it, lmbda\
                     , udp_res)

    u, d, p  = udp.split(True)
    #u.rename("u", "velocity")
    #vel_file << u
    #d.rename("d", "deformation")
    #def_file << d


    u0, d0, p0  = udp0.split(True)

    # Mesh deformation function in fluid domain

    solve(lhs(F_W) == rhs(F_W), w_, bcs_w)
    w0.assign(w_)
    ALE.move(mesh, w_)
    mesh.bounding_box_tree().build(mesh)

    drag, lift =integrateFluidStress(p, u, geometry)
    Drag.append(drag)
    Lift.append(lift)
    if MPI.rank(mpi_comm_world()) == 0:
        print "Time: ",t ," drag: ",drag, "lift: ",lift

    udp0.assign(udp)

    dis_x.append(d(coord)[0])
    dis_y.append(d(coord)[1])

    t += dt

run_time = toc()

if MPI.rank(mpi_comm_world()) == 0:
    postpro(Lift, Drag, dis_x, dis_y, time, Re, m, U_dof, run_time, mesh_cells, case = 1)
