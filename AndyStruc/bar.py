from dolfin import *


Nx = 50; Ny = 10
mesh = RectangleMesh(Point(0, 0), Point(2.5, 0.41), 10, 10, "right")
print "Plotting a RectangleMesh"

#FunctionSpaces and Functions
V = VectorFunctionSpace(mesh, "CG", 1)
W = V*V

psi, eta = TestFunctions(W)

up = Function(W)
v, u = split(up)

up1 = Function(W)
v1, u1 = split(up1)
#v1 = Function(V)
#u1 = Function(V)




#Parameters
dt = 0.01; t = 0
k = Constant(dt)
T = 1

rho_s = 1.0E3; nu_s = 0.4; mu_s = 0.5E6
E_fac = 1.4E6
lamb_s = nu_s*E_fac / ((1. + mu_s)*(1. - 2.*mu_s) )

#First Piola Stress Tensor
def P(u):
    F = Identity(2) + grad(u) #Deformation gradient
    E = 0.5*(F.T*F - Identity(2)) #Lagrangian Green strain
    J = det(F) #Determinant of Deformation gradient
    #Cauchy Stress tensor
    Sigma = 1./J * F * (lamb_s*tr(E)*Identity(2) + 2*mu_s*E) * F.T
    # P = J*Sigma*F^(-T) Gives
    P_tensor = F * (lamb_s*tr(E)*Identity(2) + 2*mu_s*E)
    return P_tensor

B = Constant((0, -2)) #2 is set in Turek Paper!!!!!!

lhs = ( rho_s/k*inner((v - v1),psi) + inner(P(0.5*(u + u1)), grad(psi)) + \
      inner((u - u1),eta)- k*inner((v + v1)/2., eta) )*dx
rhs = rho_s*inner(B, psi)*dx
F = lhs - rhs

#Boundary Conditions
#nos_V = DirichletBC(W.sub(0), ((10, 0)), "x[0] < DOLFIN_EPS")
nos_U = DirichletBC(W.sub(1), ((0, 0)), "x[0] < DOLFIN_EPS")
bcs = [nos_U]

print "Starting Newton iterations \nComputing for t = %g" % ( dt)
vel_file = File("velocity/vel.pvd")
def_file = File("deformation/def.pvd")

while t <= T:
    print "iteration for time %.2f" % t
    #time.append(t)

    J = derivative(F, up)

    problem = NonlinearVariationalProblem(F, up, bcs, J)
    solver  = NonlinearVariationalSolver(problem)

    prm = solver.parameters
    prm['newton_solver']['absolute_tolerance'] = 1E-6
    prm['newton_solver']['relative_tolerance'] = 1E-6
    prm['newton_solver']['maximum_iterations'] = 20
    prm['newton_solver']['relaxation_parameter'] = 1.0

    v_, u_ = up.split(True)

    up1.assign(up)
    #v1.assign(v_)
    #u1.assign(u_)

    def_file << u_
    vel_file << v_
    t += dt
