from dolfin import *
from solvers import Newtonsolver, Newton_manual
from postproses import postpro
from argpar import parse

args = parse()
case = args.case
T = args.T; dt = args.dt
v_deg = args.v_deg; p_deg = args.p_deg
solver = args.solver
theta = args.theta
fig = False

Um = 0
print "HERE"
print case
if case == "1":
    Um = 0.2
if case == "2":
    Um = 1.0
if case == "3":
    Um = 2.0

t = 0.0
k = Constant(dt)
H = 0.41
D = 0.1
R = D/2.
nu = 0.001
rho = 1000.
mu = rho*nu

# Refines mesh
m = "Mesh/turek1.xml"
mesh = Mesh(m)
Drag = []; Lift = []; time = []
if args.refiner != None:
    for i in range(args.refiner):
        mesh = refine(mesh)

# Fluid stress tensor
def sigma_f(p, u):
    return - p*Identity(2) + mu*(grad(u) + grad(u).T)

def eps(v):
    return 0.5*(grad(v).T + grad(v))

# Fluid stress forces
def integrateFluidStress(p, u):
  eps   = 0.5*(grad(u) + grad(u).T)
  sig   = -p*Identity(2) + 2.0*mu*eps

  traction  = dot(sig, -n)

  forceX = traction[0]*ds(1)
  forceY = traction[1]*ds(1)
  fX = assemble(forceX)
  fY = assemble(forceY)

  return fX, fY

V = VectorFunctionSpace(mesh, "CG", v_deg) # Fluid velocity
Q  = FunctionSpace(mesh, "CG", p_deg)      # Fluid Pressure

mesh_cells = mesh.num_cells()
VQ = MixedFunctionSpace([V, Q])
U_dof = VQ.dim()

## Facetfunctions ########
Inlet  = AutoSubDomain(lambda x: "on_boundary" and near(x[0], 0))
Outlet = AutoSubDomain(lambda x: "on_boundary" and near(x[0], 2.5))
Walls  = AutoSubDomain(lambda x: "on_boundary" and near(x[1], 0) or near(x[1], 0.41))

boundaries = FacetFunction("size_t",mesh)
boundaries.set_all(0)
DomainBoundary().mark(boundaries, 1)
Inlet.mark(boundaries, 2)
Outlet.mark(boundaries, 3)
Walls.mark(boundaries, 4)
ds = Measure("ds", subdomain_data = boundaries)
n = FacetNormal(mesh)

## Boundary conditions
inlet = Expression(("1.5*Um*x[1]*(H - x[1]) / pow((H/2.0), 2) * (1 - cos(t*pi/2))/2"\
                        ,"0"), t = 0.0, Um = Um, H = H)

u_inlet = DirichletBC(VQ.sub(0), inlet, boundaries, 2)
nos_geo = DirichletBC(VQ.sub(0), ((0, 0)), boundaries, 1)
nos_wall = DirichletBC(VQ.sub(0), ((0, 0)), boundaries, 4)
p_out = DirichletBC(VQ.sub(1), 0, boundaries, 3)

bcs = [u_inlet, nos_geo, nos_wall, p_out]

# Functions
phi, eta = TestFunctions(VQ)
up = Function(VQ)
u, p = split(up)

up0 = Function(VQ)
u0, p0 = split(up0)

d_up = TrialFunction(VQ) #For manual Newton solver
up_res = Function(VQ)

# Variational form
F = rho/k*inner(u - u0, phi) *dx \
+ rho*(theta*inner(dot(grad(u), u), phi) + (1 - theta)*inner(dot(grad(u0), u0), phi) )*dx \
+ inner(theta*sigma_f(p, u) + (1 - theta)*sigma_f(p0, u0), eps(phi))*dx \
- inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, phi)*ds \
- eta*div(u)*dx

J = derivative(F, up, d_up)

#Solver parameters
atol, rtol = 1e-7, 1e-7             # abs/rel tolerances
lmbda = 1.0                         # relaxation parameter
residual   = 1                      # residual (To initiate)
rel_res    = residual               # relative residual
max_it    = 100                     # max iterations
Iter = 0                            # Iteration counter

Re = Um*D*rho/mu
if MPI.rank(mpi_comm_world()) == 0:
    print "Starting Newton iterations \nComputing for t = %g" % ( dt)
    print "SOLVING FOR Re = %f" % Re #0.1 Cylinder diameter
    print "DOF = %f,  cells = %f" % (U_dof, mesh_cells)

vel_file = File("velocity/vel1y.pvd")
tic()

#Main time iterator
while t <= T:
    time.append(t)

    if t < 2:
        inlet.t = t;
    if t >= 2:
        inlet.t = 2;

    if solver == "Newton":
        Newtonsolver(F, up, bcs, J, atol, rtol, max_it, lmbda)

    if solver == "Newton2":
        Newton_manual(F, up, bcs, J, atol, rtol, max_it, lmbda, up_res)

    u_, p_ = up.split(True)
    up0.assign(up)

    #vel_file << u_

    drag, lift =integrateFluidStress(p_, u_)
    Drag.append(drag)
    Lift.append(lift)

    if MPI.rank(mpi_comm_world()) == 0:
      print "Time: ",t ," drag: ",drag, "lift: ",lift
    t += dt

run_time = toc()

if MPI.rank(mpi_comm_world()) == 0:
    postpro(Lift, Drag, time, Re, m, U_dof, run_time, mesh_cells, case)
