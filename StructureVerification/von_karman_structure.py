from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import sys
set_log_active(False)
implementation = sys.argv[1]
print implementation

mesh = Mesh("von_karman_street_FSI_structure.xml")

for coord in mesh.coordinates():
    if coord[0]==0.6 and (0.199<=coord[1]<=0.2001): # to get the point [0.2,0.6] end of bar
        print coord
        break

V = VectorFunctionSpace(mesh, "CG", 1)
if implementation == "2":
    VV=V*V
print "Dofs: ",V.dim(), "Cells:", mesh.num_cells()

BarLeftSide =  AutoSubDomain(lambda x: "on_boundary" and (( (x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.2)*(x[1] - 0.2)  < 0.0505*0.0505 )  and x[1]>=0.19 and x[1]<=0.21 and x[0]>0.2 ))

boundaries = FacetFunction("size_t",mesh)
boundaries.set_all(0)
BarLeftSide.mark(boundaries,1)
#plot(boundaries,interactive=True)

#PARAMETERS:
rho_s = 1.0E3
mu_s = 0.5E6
nu_s = 0.4
E_1 = 1.4E6
lamda = nu_s*2*mu_s/(1-2*nu_s)
g = Constant((0,-2*rho_s))
dt = 0.02
k = Constant(dt)

#Hookes
def sigma_structure(d):
    return 2*mu_s*sym(grad(d)) + lamda*tr(sym(grad(d)))*Identity(2)

#Second piola stress
def s_s_n_l(d):
    I = Identity(2)
    F = I + grad(d)
    E = 0.5*((F.T*F)-I)
    #J = det(F)
    return lamda*tr(E)*I + 2*mu_s*E




#Variational form with single space, just solving displacement
if implementation =="1":
    bc1 = DirichletBC(V, ((0,0)),boundaries, 1)
    bcs = [bc1]
    psi = TestFunction(V)
    d = Function(V)
    d0 = Function(V)
    d1 = Function(V)
    G =rho_s*((1./k**2)*inner(d - 2*d0 + d1,psi))*dx \
    + inner(s_s_n_l(0.5*(d+d1)),grad(psi))*dx - inner(g,psi)*dx

#Variational form with double spaces
elif implementation =="2":
    bc1 = DirichletBC(VV.sub(0), ((0,0)),boundaries, 1)
    bc2 = DirichletBC(VV.sub(1), ((0,0)),boundaries, 1)
    bcs = [bc1,bc2]
    psi,phi = TestFunctions(VV)
    wd = Function(VV)
    w,d = split(wd)
    w0d0 = Function(VV)
    w0,d0 = split(w0d0)

    G =rho_s*((1./k)*inner(w-w0,psi))*dx + rho_s*inner(dot(grad(0.5*(w+w0)),0.5*(w+w0)),psi)*dx + inner(s_s_n_l(0.5*(d+d0)),grad(psi))*dx \
    - inner(g,psi)*dx - dot(d-d0,phi)*dx + k*dot(0.5*(w+w0),phi)*dx
elif implementation == "3":
    bc1 = DirichletBC(V, ((0,0)),boundaries, 1)
    bcs = [bc1]
    psi = TestFunction(V)
    w = Function(V)
    w0 = Function(V)
    d0 = Function(V)
    d = d0 + w*k
    G =rho_s*((1./k)*inner(w-w0,psi))*dx + rho_s*inner(dot(grad(0.5*(w+w0)),0.5*(w+w0)),psi)*dx \
    + inner(s_s_n_l(0.5*(d+d0)),grad(psi))*dx - inner(g,psi)*dx




dis_x = []
dis_y = []

T = 10.0
t = 0
time = np.linspace(0,T,(T/dt)+1)
print "time len",len(time)

from time import sleep
while t < T:
    print "Time: ",t
    if implementation == "1":
        solve(G == 0, d, bcs, solver_parameters={"newton_solver": \
        {"relative_tolerance": 1E-9,"absolute_tolerance":1E-9,"maximum_iterations":100,"relaxation_parameter":1.0}})
        d1.assign(d0)
        d0.assign(d)
        #plot(d,mode="displacement")
        dis_x.append(d(coord)[0])
        dis_y.append(d(coord)[1])
    elif implementation == "2":
        solve(G == 0, wd, bcs, solver_parameters={"newton_solver": \
        {"relative_tolerance": 1E-9,"absolute_tolerance":1E-9,"maximum_iterations":100,"relaxation_parameter":1.0}})
        w0d0.assign(wd)
        w,d = wd.split()
        #plot(d,mode="displacement")
        dis_x.append(d(coord)[0])
        dis_y.append(d(coord)[1])
    elif implementation == "3":
        solve(G == 0, w, bcs, solver_parameters={"newton_solver": \
        {"relative_tolerance": 1E-9,"absolute_tolerance":1E-9,"maximum_iterations":100,"relaxation_parameter":1.0}})
        w0.assign(w)
        d0.assign(d)
        #plot(d0,mode="displacement")
        dis_x.append(d0(coord)[0])
        dis_y.append(d0(coord)[1])
    t += dt

print len(dis_x), len(time)
if implementation == "1":
    title = plt.title("Single space solving d")
elif implementation == "2":
    title = plt.title("Double space")
elif implementation == "3":
    title = plt.title("Single space solving w")

plt.plot(time,dis_x,);title; plt.ylabel("Displacement x");plt.xlabel("Time");plt.grid();
plt.show()
plt.plot(time,dis_y);title;plt.ylabel("Displacement y");plt.xlabel("Time");plt.grid();
plt.show()
