from dolfin import *
parameters['allow_extrapolation']=True

#mesh = Mesh("von_karman_street_FSI_fluid.xml")
mesh = Mesh("fluid_new.xml")
#plot(mesh,interactive=True)

for coord in mesh.coordinates():
    if coord[0]==0.6 and (0.199<=coord[1]<=0.2001): # to get the point [0.2,0.6] end of bar
        print coord
        break

V1 = VectorFunctionSpace(mesh, "CG", 2) # Velocity
V2 = VectorFunctionSpace(mesh, "CG", 1) # Structure deformation
Q  = FunctionSpace(mesh, "CG", 1)       # Fluid Pressure

VVQ = MixedFunctionSpace([V1,V2,Q])

# BOUNDARIES

#NOS = AutoSubDomain(lambda x: "on_boundary" and( near(x[1],0) or near(x[1], 0.41)))
Inlet = AutoSubDomain(lambda x: "on_boundary" and near(x[0],0))
Outlet = AutoSubDomain(lambda x: "on_boundary" and (near(x[0],2.5)))
Wall =  AutoSubDomain(lambda x: "on_boundary" and (near(x[1], 0.41) or near(x[1], 0)))
Bar = AutoSubDomain(lambda x: "on_boundary" and (near(x[1], 0.21)) or near(x[1], 0.19) or near(x[0], 0.6 ) )
Circle =  AutoSubDomain(lambda x: "on_boundary" and (( (x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.2)*(x[1] - 0.2)  < 0.0505*0.0505 )  ))
Barwall =  AutoSubDomain(lambda x: "on_boundary" and (( (x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.2)*(x[1] - 0.2)  < 0.0505*0.0505 )  and x[1]>=0.19 and x[1]<=0.21 and x[0]>0.2 ))

Allboundaries = DomainBoundary()

boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
Wall.mark(boundaries, 2)
Inlet.mark(boundaries, 3)
Outlet.mark(boundaries, 4)
Bar.mark(boundaries, 5)
Circle.mark(boundaries, 6)
Barwall.mark(boundaries, 7)
#test = File("Facet.pvd")
#test << boundaries
#plot(boundaries,interactive=True)

#BOUNDARY CONDITIONS

################################################


Um = 1.0
H = 0.41
L = 2.5
D = 0.1

# t = 2.0 implies steady flow
inlet = Expression(("1.5*Um*x[1]*(H - x[1]) / pow((H/2.0), 2) * (1 - cos(t*pi/2))/2"\
,"0"), t = 0.0, Um = Um, H = H)

#velocity conditions
#u_inlet  = DirichletBC(VVQ.sub(0), inlet, boundaries, 3)
#u_wall   = DirichletBC(VVQ.sub(0), ((0, 0)), boundaries, 2)
#u_circ   = DirichletBC(VVQ.sub(0), ((0, 0)), boundaries, 6) #No slip on geometry in fluid
#u_barwall = DirichletBC(VVQ.sub(0), ((0, 0)), boundaries, 7)


ds = Measure("ds", subdomain_data = boundaries)
dS = Measure("dS", subdomain_data = boundaries)
n = FacetNormal(mesh)

Bar_area = AutoSubDomain(lambda x: (0.19 <= x[1] <= 0.21) and 0.24<= x[0] <= 0.6) # only the "flag" or "bar"

#domains = CellFunction("size_t",mesh)
#domains.set_all(1)
#Bar_area.mark(domains,2) #Overwrites structure domain
#dx = Measure("dx",subdomain_data=domains)
#plot(domains,interactive = True)

#test = File("cellfunc.pvd")
#test << domains

u_inlet  = DirichletBC(V1, ((1, 0)), boundaries, 3)
u_wall   = DirichletBC(V1, ((0, 0)), boundaries, 2)
u_circ   = DirichletBC(V1, ((0, 0)), boundaries, 6) #No slip on geometry in fluid
u_bar = DirichletBC(V1, ((0, 0)), boundaries, 5)
u_barwall = DirichletBC(V1, ((0, 0)), boundaries, 7)

#u_inlet  = DirichletBC(Q, 1, boundaries, 3)
#u_wall   = DirichletBC(Q, 0, boundaries, 2)
#u_circ   = DirichletBC(Q, 0, boundaries, 6) #No slip on geometry in fluid
#u_bar = DirichletBC(Q, 0, boundaries, 5)

bcs = [u_inlet, u_wall, u_circ, u_bar, u_barwall]

u = TrialFunction(V1)
v = TestFunction(V1)

u_ = Function(V1)

a = inner(grad(u), grad(v))*dx
L = inner(Constant((0,0)), v) *dx

solve(a == L, u_, bcs)

test = File("veri.pvd")
test << u_
