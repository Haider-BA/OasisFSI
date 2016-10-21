from dolfin import *
#mesh = Mesh("von_karman_street_FSI_fluid.xml")
mesh = Mesh("fluid_new.xml")
#plot(mesh,interactive=True)

V1 = VectorFunctionSpace(mesh, "CG", 2) # Velocity
V2 = VectorFunctionSpace(mesh, "CG", 1) # Mesh movement
V3 = VectorFunctionSpace(mesh, "CG", 1) # Structure deformation
Q  = FunctionSpace(mesh, "CG", 1)       # Fluid Pressure

VVQ = MixedFunctionSpace([V1,V2,V3,Q])

# BOUNDARIES

#NOS = AutoSubDomain(lambda x: "on_boundary" and( near(x[1],0) or near(x[1], 0.41)))
Inlet = AutoSubDomain(lambda x: "on_boundary" and near(x[0],0))
Outlet = AutoSubDomain(lambda x: "on_boundary" and (near(x[0],2.5)))
Wall =  AutoSubDomain(lambda x: "on_boundary" and (near(x[1], 0.41) or near(x[1], 0)))
Bar = AutoSubDomain(lambda x: "on_boundary" and (near(x[1], 0.21)) or near(x[1], 0.19) or near(x[0], 0.6 ) )
Circle =  AutoSubDomain(lambda x: "on_boundary" and (( (x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.2)*(x[1] - 0.2)  < 0.0505*0.0505 )  ))
Barwall =  AutoSubDomain(lambda x: "on_boundary" and (( (x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.2)*(x[1] - 0.2)  < 0.0505*0.0505 )  and x[1]>=0.19 and x[1]<=0.21 and x[0]>0.2 ))

Allboundaries = DomainBoundary()

boundaries = FacetFunction("size_t",mesh)
boundaries.set_all(0)
#Allboundaries.mark(boundaries, 1)
Wall.mark(boundaries, 2)
Inlet.mark(boundaries, 3)
Outlet.mark(boundaries, 4)
Bar.mark(boundaries, 5)
Circle.mark(boundaries, 6)
Barwall.mark(boundaries, 7)
#plot(boundaries,interactive=True)


ds = Measure("ds", subdomain_data = boundaries)
dS = Measure("dS", subdomain_data = boundaries)
n = FacetNormal(mesh)

#BOUNDARY CONDITIONS

################################################


Um = 0.2
H = 0.41
L = 2.5
D = 0.1

# t = 2.0 implies steady flow
inlet = Expression(("1.5*Um*x[1]*(H - x[1]) / pow((H/2.0), 2) * (1 - cos(t*pi/2))/2"\
,"0"), t = 2.0, Um = Um, H = H)

#velocity conditions
u_inlet  = DirichletBC(VVQ.sub(0), inlet, boundaries, 3)
u_wall   = DirichletBC(VVQ.sub(0), ((0, 0)), boundaries, 2)
u_circ   = DirichletBC(VVQ.sub(0), ((0, 0)), boundaries, 6) #No slip on geometry in fluid
u_barwall = DirichletBC(VVQ.sub(0), ((0, 0)), boundaries, 7)

#Mesh velocity conditions
w_wall    = DirichletBC(VVQ.sub(1), ((0, 0)), boundaries, 2)
w_inlet   = DirichletBC(VVQ.sub(1), ((0, 0)), boundaries, 3)
w_outlet  = DirichletBC(VVQ.sub(1), ((0, 0)), boundaries, 4)
w_cirlce  = DirichletBC(VVQ.sub(1), ((0, 0)), boundaries, 6)
w_barwall = DirichletBC(VVQ.sub(1), ((0, 0)), boundaries, 7)

#Deformation conditions
d_barwall = DirichletBC(VVQ.sub(2), ((0, 0)), boundaries, 7)

#Pressure Conditions
p_out = DirichletBC(VVQ.sub(3), 0, boundaries, 4)

#Assemble boundary conditions
bcs = [u_inlet, u_wall, u_circ, u_barwall, \
       w_inlet, w_wall, w_cirlce, w_outlet, w_barwall, \
       d_barwall]
# AREAS

Bar_area = AutoSubDomain(lambda x: (0.19 <= x[1] <= 0.21) and 0.24<= x[0] <= 0.6) # only the "flag" or "bar"

domains = CellFunction("size_t",mesh)
domains.set_all(1)
Bar_area.mark(domains,2) #Overwrites structure domain
dx = Measure("dx",subdomain_data=domains)
plot(domains,interactive = True)


# TEST TRIAL FUNCTIONS
phi, psi, gamma, eta = TestFunctions(VVQ)
#u, w, p = TrialFunctions(VVQ)
uwdp = Function(VVQ)
u, w, d, p  = split(uwdp)
uwdp0 = Function(VVQ)
u0, w0, d0, p0  = split(uwdp0)
#u1 = Function(V1) Piccard

#EkPa = '62500
#E = Constant(float(EkPa))

#Fluid properties
rho_f   = Constant(1000.0)
mu_f    = Constant(1.0)
nu_f = mu_f/rho_f

#Structure properties

rho_s = 1.0E3
mu_s = 0.5E6
nu_s = 0.4
E_1 = 1.4E6
lamda = nu_s*2*mu_s/(1-2*nu_s)
g = Constant((0,-2*rho_s))
dt = 0.01
k = Constant(dt)

Re = Um*D/nu_f
print "SOLVING FOR Re = %f" % Re #0.1 Cylinder diameter



def Venant_Kirchhof(d):
    I = Identity(2)
    F = I - grad(d)
    J = det(F)
    E = 0.5*((inv(F.T)*inv(F))-I)
    return inv(F)*(2.*mu_s*E + lamda*tr(E)*I)*inv(F.T)

def sigma_f(p, u):
  return - p*Identity(2) + mu_f*(grad(u) + grad(u).T)


# Fluid variational form
Fluid_momentum = (rho_f/k)*inner(u - u0, phi)*dx + rho_f*inner(grad(u)*u, phi)*dx+ \
    inner(sigma_f(p, u), grad(phi))*dx

Fluid_continuity = eta*div(u)*dx

#############################
I = Identity(2)
F_ = I - grad(d)
J_ = det(F_)
F_1 = I - grad(d0)
J_1 = det(F_1)
theta = Constant(0.593)

# Structure Variational form
Solid_momentum = ( J_*rho_s/k*inner(w - w0, psi) \
    + rho_s*( J_*theta*inner(dot(grad(w), w), psi) + J_1*(1 - theta)*inner(dot(grad(w0), w0), psi) ) \
    + inner(J_*theta*Venant_Kirchhof(d) + (1 - theta)*J_1*Venant_Kirchhof(d0) , grad(psi))  \
    - (theta*J_*inner(g, psi) + (1-theta)*J_1*inner(g, psi) ) ) * dx(2)

Solid_deformation = dot(d - d0 + k*(theta*dot(grad(d), w) + (1-theta)*dot(grad(d0), w0) ) \
    - k*(theta*w + (1 -theta)*w0 ), gamma)  * dx(2)

# Mesh velocity function in fluid domain
d_smooth = inner(grad(u),grad(phi))*dx(1) - inner(grad(u("-"))*n("-"),phi("-"))*dS(5)

F = Fluid_momentum + Fluid_continuity \
  + Solid_momentum + Solid_deformation \
  + d_smooth

#u_file = File("mvelocity/velocity.pvd")
#u_file << u0
t = 0
T = 1
time = []; Lift = []; Drag = []
vel_file = File("velocity/velocity.pvd")
while t <= T:
    time.append(t)

    if t < 2:
        inlet.t = t;
    if t >= 2:
        inlet = inlet_steady;

    J = derivative(F, uwdp)

    problem = NonlinearVariationalProblem(F, uwdp, bcs, J)
    solver  = NonlinearVariationalSolver(problem)

    prm = solver.parameters
    prm['newton_solver']['absolute_tolerance'] = 1E-10
    prm['newton_solver']['relative_tolerance'] = 1E-10
    prm['newton_solver']['maximum_iterations'] = 5
    prm['newton_solver']['relaxation_parameter'] = 1.0


    solver.solve()

    u_, w_, d_, p_ = uwdp.split(True)
    u0, w0, d0, p0 = uwdp0.split(True)

    d_disp.vector()[:] = d_.vector()[:] - d0.vector()[:]
    ALE.move(mesh, d_disp)
    mesh.bounding_box_tree().build(mesh)
    uwdp0.assign(uwdp)

    u0, w0, d0, p0 = uwdp0.split(True)
    dis_x.append(d0(coord)[0])
    dis_y.append(d0(coord)[1])


    drag, lift =integrateFluidStress(p_, u_)
    if MPI.rank(mpi_comm_world()) == 0:
      print "Time: ",t ," drag: ",drag, "lift: ",lift
    Drag.append(drag)
    Lift.append(lift)

    t += dt



plt.figure(1)
plt.title("LIFT \n Re = %.1f, dofs = %d, cells = %d" % (Re, U_dof, mesh_cells))
plt.xlabel("Time Seconds")
plt.ylabel("Lift force Newton")
plt.plot(time, Lift, label='dt  %g' % dt)
plt.legend(loc=4)
plt.savefig("lift.png")

plt.figure(2)
plt.title("DRAG \n Re = %.1f, dofs = %d, cells = %d" % (Re, U_dof, mesh_cells))
plt.xlabel("Time Seconds")
plt.ylabel("Drag force Newton")
plt.plot(time, Drag, label='dt  %g' % dt)
plt.legend(loc=4)
plt.savefig("drag.png")
#plt.show()


#drag, lift = integrateFluidStress(p_, u_)
#print('Drag = %f, Lift = %f' \
	#%  (drag, lift))
