from dolfin import *
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
group.add_argument("-d_deg",  type=int, help="Set degree of velocity                     --> Default=2", default=1)
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
fig = False

parameters['allow_extrapolation']=True

#mesh = Mesh("von_karman_street_FSI_fluid.xml")
mesh = Mesh("fluid_new.xml")
#plot(mesh,interactive=True)
if args.refiner != None:
    for i in range(args.refiner):
        mesh = refine(mesh)

for coord in mesh.coordinates():
    if coord[0]==0.6 and (0.199<=coord[1]<=0.2001): # to get the point [0.2,0.6] end of bar
        print coord
        break

V1 = VectorFunctionSpace(mesh, "CG", v_deg) # Velocity
V2 = VectorFunctionSpace(mesh, "CG", d_deg) # Structure deformation
Q  = FunctionSpace(mesh, "CG", p_deg)       # Fluid Pressure

VVQ = MixedFunctionSpace([V1,V2,Q])

# BOUNDARIES

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

# All boundaries execpt flag, for
Neumann = FacetFunction("size_t",mesh, 0)
DomainBoundary().mark(Neumann, 1)
Barwall.mark(Neumann, 0)

# Flag boundary, for balance of momentum on interface
interface = FacetFunction("size_t", mesh)
interface.set_all(0)
Bar.mark(interface, 5)

# Full geometry for force integral
geometry = FacetFunction("size_t",mesh, 0)
Bar.mark(geometry, 1)
Circle.mark(geometry, 1)
Barwall.mark(geometry, 0)

ds = Measure("ds", subdomain_data = Neumann)
dS = Measure("dS", subdomain_data = interface) # For interface (interior)
n = FacetNormal(mesh)

# Parameters
Um = 0.2
H = 0.41
L = 2.5
D = 0.1

#Fluid properties
rho_f   = Constant(1000.0)
mu_f    = Constant(1.0)
nu_f = mu_f/rho_f

#Structure properties
#FSI 1
rho_s = 1E3
mu_s = 0.5E6
nu_s = 0.4
E_1 = 1.4E6
lamda = nu_s*2*mu_s/(1-2*nu_s)
g = Constant((0,-2*rho_s))

# AREAS
Bar_area = AutoSubDomain(lambda x: (0.19 <= x[1] <= 0.21) and 0.24<= x[0] <= 0.6) # only the "flag" or "bar"
boundary_parts = FacetFunction("size_t", mesh, 0)
Bar_area.mark(boundary_parts, 1)

domains = CellFunction("size_t",mesh)
domains.set_all(1)
Bar_area.mark(domains,2) #Overwrites structure domain
dx = Measure("dx",subdomain_data=domains)

# Boundary conditions
inlet = Expression(("1.5*Um*x[1]*(H - x[1]) / pow((H/2.0), 2)" ,"0"), Um = Um, H = H)

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
p_struc   = DirichletBC(VVQ.sub(2), 0, boundary_parts, 1)

# Assemble boundary conditions
bcs = [u_inlet, u_wall, u_circ, u_barwall, \
       d_inlet, d_wall, d_circ, d_barwall, d_outlet,\
       p_out, p_barwall, p_struc]



# TEST TRIAL FUNCTIONS
psi, gamma, eta = TestFunctions(VVQ)

udp = Function(VVQ)
u, d, p  = split(udp)

udp0 = Function(VVQ)
u0, d0, p0  = split(udp0)

d_disp = Function(V2)

Re = Um*D/nu_f
print "SOLVING FOR Re = %f" % Re #0.1 Cylinder diameter

def Venant_Kirchhof(d):
    I = Identity(2)
    F = I - grad(d)
    J = det(F)
    E = 0.5*((inv(F.T)*inv(F))-I)
    return inv(F)*(2.*mu_s*E + lamda*tr(E)*I)*inv(F.T)

def integrateFluidStress(p, u, geo, n):
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
Fluid_momentum = rho_f*(theta*inner(dot(grad(u), u), psi) + (1 - theta)*inner(dot(grad(u0), u0), psi) )*dx(1) \
    + inner(theta*sigma_f(p, u) + (1 - theta)*sigma_f(p0, u0), eps(psi))*dx \
    - inner(theta*sigma_f(p, u)*n + (1 - theta)*sigma_f(p0, u0)*n, psi)*ds(1) \

Fluid_continuity = eta*div(u)*dx(1)

# Structure Variational form

I = Identity(2)
F_ = I - grad(d)
J_ = det(F_)
F_1 = I - grad(d0)
J_1 = det(F_1)



Solid_momentum = ( inner(J_*theta*Venant_Kirchhof(d) + (1 - theta)*J_1*Venant_Kirchhof(d0) , grad(psi))  \
                - (theta*J_*inner(g, psi) + (1-theta)*J_1*inner(g, psi) ) ) * dx(2)

Solid_deformation = inner(theta*u + (1 -theta)*u0, gamma)  * dx(2)

# Mesh velocity function in fluid domain
d_smooth = inner(grad(d),grad(gamma))*dx(1) #- inner(grad(u("-"))*n("-"),phi("-"))*dS(5)

#Conservation of dynamics (CHECK SIGNS!!)
#dynamic = - inner(Venant_Kirchhof(d('+'))*n('+'), psi('+'))*dS(5) - inner(sigma_f(p('-'), u('-'))*n('-') ,psi('-'))*dS(5)
dynamic = - inner(Venant_Kirchhof(d('-'))*n('-'), psi('-'))*dS(5) - inner(sigma_f(p('+'), u('+'))*n('+') ,psi('+'))*dS(5)
#dynamic = - inner(Venant_Kirchhof(d)*n, psi)*dS(5) - inner(sigma_f(p, u)*n ,psi)*dS(5)

F = Fluid_momentum + Fluid_continuity \
  + Solid_momentum + Solid_deformation \
  + d_smooth + dynamic

t = 0

time = []; dis_x = []; dis_y = []
#vel_file = File("./velocity/velocity.pvd")
#vel_file << u0

J = derivative(F, udp)

problem = NonlinearVariationalProblem(F, udp, bcs, J)
sol  = NonlinearVariationalSolver(problem)

prm = sol.parameters
#info(prm,True)  #get full info on the parameters
prm['nonlinear_solver'] = 'newton'
prm['newton_solver']['absolute_tolerance'] = 1E-10
prm['newton_solver']['relative_tolerance'] = 1E-10
prm['newton_solver']['maximum_iterations'] = 10
prm['newton_solver']['relaxation_parameter'] = 1.0
prm['newton_solver']['linear_solver'] = 'mumps'


sol.solve()

u, d, p  = udp.split(True)

u0, d0, p0  = udp0.split(True)

d_disp.vector()[:] = d.vector()[:] - d0.vector()[:]

ALE.move(mesh, d_disp)
mesh.bounding_box_tree().build(mesh)
udp0.assign(udp)

drag, lift = integrateFluidStress(p, u, geometry, n)

#vel_file << u
print "DOFS", VVQ.dim()
print "Cells", mesh.num_vertices()
print "Discretization theta = %g" % theta
print "Displacement x", d(coord)[0]
print "Displacement y", d(coord)[1]
print "Lift %g" % lift
print "Drag %g" % drag
