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
if implementation == "2" or "3":
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
dt = float(sys.argv[2]) #FOR FASTER TESTNING!! REMOVE AFTER TESTS!!!!
#dt = 0.02
k = Constant(dt)

#Hookes
def sigma_structure(d):
    return 2*mu_s*sym(grad(d)) + lamda*tr(sym(grad(d)))*Identity(2)

#Cauchy Stress Tensor sigma
def Cauchy(d):
    I = Identity(2)
    F = I + grad(d)
    J = det(F)
    E = 0.5*((F.T*F)-I)
    return 1./J*F*(lamda*tr(E)*I + 2*mu_s*E )*F.T

def Venant_Kirchhof(d):
    I = Identity(2)
    F = I - grad(d)
    J = det(F)
    E = 0.5*((inv(F.T)*inv(F))-I)
    return J*inv(F)*(2.*mu_s*E + lamda*tr(E)*I)*inv(F.T)

#First Piola stress tensor P = J*Sigma*F^(-T)
def Piola1(d):
    I = Identity(2)
    F = I + grad(d)
    J = det(F)
    E = 0.5*((F.T*F)-I)
    return F*(lamda*tr(E)*I + 2*mu_s*E )



#Second piola stress S = J*F^(-1)*sigma*F^(-T)
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

#Full Eulerian description  Full Implicit
elif implementation == "3":
    bc1 = DirichletBC(VV.sub(0), ((0,0)),boundaries, 1)
    bc2 = DirichletBC(VV.sub(1), ((0,0)),boundaries, 1)
    bcs = [bc1,bc2]
    psi, phi = TestFunctions(VV)
    wd = Function(VV)
    w, d = split(wd)
    w0d0 = Function(VV)
    w0, d0 = split(w0d0)
    d_disp = Function(V)


    #Fully implicit, gives significant damping to the system
    I = Identity(2)
    F_ = I - grad(d) #d here must be the same as in variational formula
    J_ = det(F_)

    G = ( J_*rho_s/k*inner(w - w0, psi) + J_*rho_s*inner(dot(grad(w), w), psi) \
    + inner(Venant_Kirchhof(d), grad(psi)) - J_*inner(g, psi) ) * dx \
    + inner((d - d0) + k*dot(grad(d), w) - k*w, phi) * dx

#Full Eulerian description Full Cran-Nic
elif implementation == "4":
        bc1 = DirichletBC(VV.sub(0), ((0,0)),boundaries, 1)
        bc2 = DirichletBC(VV.sub(1), ((0,0)),boundaries, 1)
        bcs = [bc1,bc2]
        psi, phi = TestFunctions(VV)
        wd = Function(VV)
        w, d = split(wd)
        w0d0 = Function(VV)
        w0, d0 = split(w0d0)
        d_disp = Function(V)

        #Full CrancNic Diverge after 2 period
        I = Identity(2)
        F_ = I - grad(0.5*(d+d0)) #d here must be the same as in variational formula
        J_ = det(F_)

        ### TEST CODE
        G = ( J_*rho_s/k*inner(w - w0, psi) + J_*rho_s*inner(dot(grad(0.5*(w+w0)), 0.5*(w + w0)), psi) \
        + inner(0.5*(Venant_Kirchhof(d) + Venant_Kirchhof(d0)), grad(psi)) - J_*inner(g, psi) ) * dx \
        + (dot(d - d0, phi) + k*dot(dot(grad(d), w), phi ) \
        - k*dot(0.5*(w+w0), phi) ) *dx
        #+ inner(d - d0 + k*dot(grad(0.5*(d+d0)), 0.5*(w+w0)) - k*w) , phi) * dx

        ### Middled d inside sigma arguent, q(0.5*(d+d0))  t > 2
        #G = ( J_*rho_s/k*inner(w - w0, psi) + J_*rho_s*inner(dot(grad(0.5*(w+w0)), 0.5*(w + w0)), psi) \
        #+ inner(Venant_Kirchhof(0.5*(d+d0)), grad(psi)) - J_*inner(g, psi) ) * dx \
        #+ (dot(d - d0, phi) + k*dot(dot(grad(0.5*(d+d0)), 0.5*(w+w0)), phi ) \
        #- k*dot(w, phi) ) *dx
        #+ inner(d - d0 + k*dot(grad(0.5*(d+d0)), 0.5*(w+w0)) - k*w) , phi) * dx

        ### Middled d by 0.5*(q(d) + q(d0) ) t > 2
        #G = ( J_*rho_s/k*inner(w - w0, psi) + J_*rho_s*inner(dot(grad(0.5*(w+w0)), 0.5*(w + w0)), psi) \
        #+ inner(0.5*(Venant_Kirchhof(d) + Venant_Kirchhof(d0)), grad(psi)) - J_*inner(g, psi) ) * dx \
        #+ (dot(d - d0, phi) + k*dot(dot(grad(0.5*(d+d0)), 0.5*(w+w0)), phi ) \
        #- k*dot(w, phi) ) *dx
        #+ inner(d - d0 + k*dot(grad(0.5*(d+d0)), 0.5*(w+w0)) - k*w) , phi) * dx

        # FULL CRAC DIVERGES,
        #G = ( J_*rho_s/k*inner(w - w0, psi) + J_*rho_s*inner(dot(grad(0.5*(w+w0)), 0.5*(w + w0)), psi) \
        #+ inner(0.5*(Venant_Kirchhof(d) + Venant_Kirchhof(d0)), grad(psi)) - J_*inner(g, psi) ) * dx \
        #+ inner(d - d0 + k*dot(grad(0.5*(d+d0)), 0.5*(w+w0)) - k*0.5*(w+w0) , phi) * dx


#Full Eulerian description Forward Euler
elif implementation == "5":
        bc1 = DirichletBC(VV.sub(0), ((0,0)),boundaries, 1)
        bc2 = DirichletBC(VV.sub(1), ((0,0)),boundaries, 1)
        bcs = [bc1,bc2]
        psi, phi = TestFunctions(VV)
        wd = Function(VV)
        w, d = split(wd)
        w0d0 = Function(VV)
        w0, d0 = split(w0d0)
        d_disp = Function(V)

        #Full CrancNic Diverge after 2 period
        I = Identity(2)
        F_ = I - grad(d0) #d here must be the same as in variational formula
        J_ = det(F_)


        G = ( J_*rho_s/k*inner(w - w0, psi) + J_*rho_s*inner(dot(grad(w0), w0), psi) \
        + inner(Venant_Kirchhof(d0), grad(psi)) - J_*inner(g, psi) ) * dx \
        + inner(d - d0 + k*dot(grad(d0), w0) - k*w0 , phi) * dx

elif implementation == "6":
    bc1 = DirichletBC(V, ((0,0)),boundaries, 1)
    bcs = [bc1]
    psi = TestFunction(V)
    d = Function(V)
    d0 = Function(V)
    d1 = Function(V)
    d_disp = Function(V)
    """
    G = rho_s*((1./k**2)*inner(d - 2*d0 + d1,psi))*dx \
    + rho_s*inner(dot(grad(0.5*(d + d1) ), (d - d1)/(2.*k)), psi)*dx \
    + inner(Venant_Kirchhof(0.5*(d + d1)), grad(psi))*dx - inner(g,psi)*dx
    """
    I = Identity(2)
    F_ = I - grad(0.5*(d + d1)) #d here must be the same as in variational formula
    J_ = det(F_)

    G = J_*rho_s*((1./k**2)*inner(d - 2*d0 + d1, psi))*dx \
    + J_*rho_s*inner(dot(grad(0.5*(d + d1) ), (d - d1)/(2.*k)) ,psi)*dx \
    + inner(Venant_Kirchhof(0.5*(d + d1) ),grad(psi))*dx - J_*inner(g,psi)*dx


dis_x = []
dis_y = []
time = []
T = 5;
t = 0

print "time len",len(time)
###### REMOVE BEFORE MERGE
dis_1 = File("./displace_1/dis.pvd")
dis_2 = File("./displace_2/dis.pvd")
dis_3 = File("./displace_3/dis.pvd")
dis_4 = File("./displace_4/dis.pvd")
dis_5 = File("./displace_5/dis.pvd")
dis_6 = File("./displace_6/dis.pvd")

from time import sleep
count = 0
while t < T:
    print "Time: ",t
    time.append(t)
    if implementation == "1":
        solve(G == 0, d, bcs, solver_parameters={"newton_solver": \
        {"relative_tolerance": 1E-9,"absolute_tolerance":1E-9,"maximum_iterations":100,"relaxation_parameter":1.0}})
        d1.assign(d0)
        d0.assign(d)
        dis_1 << d
        #plot(d,mode="displacement")
        dis_x.append(d(coord)[0])
        dis_y.append(d(coord)[1])

    elif implementation == "2":
        solve(G == 0, wd, bcs, solver_parameters={"newton_solver": \
        {"relative_tolerance": 1E-9,"absolute_tolerance":1E-9,"maximum_iterations":100,"relaxation_parameter":1.0}})
        w0d0.assign(wd)
        w,d = wd.split()
        dis_2 << d
        #plot(d,mode="displacement")
        dis_x.append(d(coord)[0])
        dis_y.append(d(coord)[1])

    elif implementation == "3":
        solve(G == 0, wd, bcs, solver_parameters={"newton_solver": \
        {"relative_tolerance": 1E-9,"absolute_tolerance":1E-9,"maximum_iterations":100,"relaxation_parameter":1.0}})
        w0d0.assign(wd)
        w0,d0 = w0d0.split(True)
        dis_x.append(d0(coord)[0])
        dis_y.append(d0(coord)[1])
        d_disp.vector()[:] = w0.vector()[:]*float(k)
        dis_3 << d0
        ALE.move(mesh, d_disp)
        mesh.bounding_box_tree().build(mesh)


    elif implementation == "4":
        solve(G == 0, wd, bcs, solver_parameters={"newton_solver": \
        {"relative_tolerance": 1E-9,"absolute_tolerance":1E-9,"maximum_iterations":100,"relaxation_parameter":1.0}})
        w0d0.assign(wd)
        w0,d0 = w0d0.split(True)
        dis_x.append(d0(coord)[0])
        dis_y.append(d0(coord)[1])
        d_disp.vector()[:] = w0.vector()[:]*float(k)
        dis_4 << d0
        ALE.move(mesh, d_disp)
        mesh.bounding_box_tree().build(mesh)

    elif implementation == "5":
        solve(G == 0, wd, bcs, solver_parameters={"newton_solver": \
        {"relative_tolerance": 1E-9,"absolute_tolerance":1E-9,"maximum_iterations":100,"relaxation_parameter":1.0}})
        w0d0.assign(wd)
        w0,d0 = w0d0.split(True)
        dis_x.append(d0(coord)[0])
        dis_y.append(d0(coord)[1])
        d_disp.vector()[:] = w0.vector()[:]*float(k)
        dis_5 << d0
        ALE.move(mesh, d_disp)
        mesh.bounding_box_tree().build(mesh)

    elif implementation == "6":
        solve(G == 0, d, bcs, solver_parameters={"newton_solver": \
        {"relative_tolerance": 1E-9,"absolute_tolerance":1E-9,"maximum_iterations":100,"relaxation_parameter":1.0}})
        d1.assign(d0)
        d0.assign(d)
        dis_6 << d
        d_disp.vector()[:] = d0.vector()[:] - d1.vector()[:]
        ALE.move(mesh, d_disp)
        mesh.bounding_box_tree().build(mesh)
        dis_x.append(d(coord)[0])
        dis_y.append(d(coord)[1])


    count += 1
    t += dt


print len(dis_x), len(time)
if implementation == "1":
    title = plt.title("Single space solving d")
elif implementation == "2":
    title = plt.title("Double space")
elif implementation == "3":
    title = plt.title("Eulerian")
elif implementation == "4":
    title = plt.title("Eulerian cranc-adam/bash combined Oasis style")


plt.figure(1)
plt.plot(time,dis_x,);title; plt.ylabel("Displacement x");plt.xlabel("Time");plt.grid();
plt.savefig("Run_x.jpg")
#plt.show()
plt.figure(2)
plt.plot(time,dis_y);title;plt.ylabel("Displacement y");plt.xlabel("Time");plt.grid();
plt.savefig("Run_y.jpg")
#plt.show()
