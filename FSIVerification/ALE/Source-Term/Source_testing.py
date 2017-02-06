from ufl import *
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
set_log_active(False)

def Run(N, dt):
    mesh = UnitSquareMesh(N,N)
    x = SpatialCoordinate(mesh)
    W = VectorFunctionSpace(mesh, "CG", 2)

    u = TrialFunction(W)
    u0 = Function(W)
    u_sol = Function(W)
    v = TestFunction(W)

    u_e = Expression(("x[0]*x[0]*t","x[1]*x[1]*t"), t = 0)
    #u_e = Expression(("sin(pow(x[0],2))*t","cos(pow(x[1],2))*t"),t=0)
    #u_e = Expression(("sin(pow(x[0],2))*t","cos(pow(x[1],2))*t"),t=0)

    #u0 = interpolate(u_e,W)

    #f = Expression(("-x[0]*x[0]*sin(t) - 2*cos(t)","0"), t = 0)
    k = Constant(dt)
    t = variable(Constant(dt))
    u_x = x[0]*x[0]*t
    u_y = x[1]*x[1]*t
    #u_x = sin(x[0]**2)*t
    #u_y = cos(x[1]**2)*t
    u_vec = as_vector([u_x, u_y])

    #Sourceterm for heatequation
    F = diff(u_vec, t) - div(grad(u_vec))
    bcs = DirichletBC(W, u_e, "on_boundary")
    #bc1 = DirichletBC(W, u_e, "x[0] < DOLFIN_EPS")
    #bc2 = DirichletBC(W, u_e, "x[0] > 1 - DOLFIN_EPS")
    #bcs = [bc1, bc2]

    a = 1./k*inner(u, v)*dx + inner(grad(u), grad(v))*dx
    L = 1./k*inner(u0, v)*dx + inner(F, v)*dx
    #L = 1./k*inner(u0, v)*dx + inner(F, v)*dx

    T = 0.1
    t_step = dt
    error=[]
    while t_step <= T:
        #print "Solving for time %f" % t_step
        u_e.t = t_step
        #f.t = t_step
        t = variable(Constant(t_step))

        solve(a == L, u_sol, bcs)
        u0.assign(u_sol)
        #print "ERRORNORM", errornorm(u_e, u_sol, norm_type="l2", degree_rise = 3)
        error.append(errornorm(u_e, u_sol, norm_type="l2", degree_rise = 3))

        t_step += dt

    L2 = np.mean(error)
    E.append(L2)
    h.append(mesh.hmin())

#N = [8]
#dt = [0.001, 0.0005, 0.00025, 0.000125, 0.000125/2]
dt = [0.0001]
N = [2**i for i in range(1, 5)]
E = []; h = []
for n in N:
    print "Solving for %d" % n
    for t in dt:
        Run(n, t)

print "ERROR"
for i in E:
    print i

print "Convergence rate"
for i in range(len(E) - 1):
    r = np.log(E[i+1]/E[i] ) / np.log(h[i+1]/h[i])
    #r = np.log(E[i+1]/E[i] ) / np.log(dt[i+1]/dt[i])
    print r


"""
mesh_points=mesh.coordinates()
x0=mesh_points[:,0]
y = u_sol.vector().array()[:]
exact = u_e.vector().array()[:]
print len(x0), len(y), len(exact)
plt.figure(1)
plt.plot(x0, y, label="Numerical")
plt.plot(x0, exact, label="Exact")
plt.legend()
plt.show()
"""
