from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import sys
set_log_active(False)
"""if len(sys.argv) != 2:
    print "Usage: %s [implementation (1/2/3)]" % sys.argv[0]
    sys.exit(0)"""
implementation = sys.argv[1]
print implementation

mesh = Mesh("von_karman_street_FSI_structure.xml")

for coord in mesh.coordinates():
    if coord[0]==0.6 and (0.199<=coord[1]<=0.2001): # to get the point [0.2,0.6] end of bar
        print coord
        break

V = VectorFunctionSpace(mesh, "CG", 1)
if implementation == "1":
    VV=V*V
print "Dofs: ",V.dim(), "Cells:", mesh.num_cells()

BarLeftSide =  AutoSubDomain(lambda x: "on_boundary" and (( (x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.2)*(x[1] - 0.2)  < 0.0505*0.0505 )  and x[1]>=0.19 and x[1]<=0.21 and x[0]>0.2 ))

boundaries = FacetFunction("size_t",mesh)
boundaries.set_all(0)
BarLeftSide.mark(boundaries,1)
#plot(boundaries,interactive=True)

#PARAMETERS:
#CSM 1
rho_s = 1.0E3
mu_s = 0.5E6
nu_s = 0.4
E_1 = 1.4E6
lamda = nu_s*2*mu_s/(1-2*nu_s)

#CSM 1
rho_s = 1.0E3
mu_s = 2.0E6
nu_s = 0.4
E_1 = 1.4E6
lamda = nu_s*2*mu_s/(1-2*nu_s)


g = Constant((0,-2*rho_s))




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
    return inv(F)*(2.*mu_s*E + lamda*tr(E)*I)*inv(F.T)
for i in [0.593]:
    print "TESTING TOF THETA = %.2f" % i
    #Full Eulerian formulation with theta-rule scheme
    if implementation == "1":
        bc1 = DirichletBC(VV.sub(0), ((0,0)),boundaries, 1)
        bc2 = DirichletBC(VV.sub(1), ((0,0)),boundaries, 1)
        bcs = [bc1,bc2]
        psi, phi = TestFunctions(VV)
        wd = Function(VV)
        w, d = split(wd)
        w0d0 = Function(VV)
        w0, d0 = split(w0d0)
        d_disp = Function(V)

        I = Identity(2)
        F_ = I - grad(d) #d here must be the same as in variational formula
        J_ = det(F_)
        F_1 = I - grad(d0)
        J_1 = det(F_1)
        #theta = Constant(0.593)
        theta = Constant(i)


        G = (rho_s*( J_*theta*inner(dot(grad(w), w), psi) + J_1*(1 - theta)*inner(dot(grad(w0), w0), psi) ) \
        + inner(J_*theta*Venant_Kirchhof(d) + (1 - theta)*J_1*Venant_Kirchhof(d0) , grad(psi))  \
        - (theta*J_*inner(g, psi) + (1-theta)*J_1*inner(g, psi) ) ) * dx \
        + (dot((theta*dot(grad(d), w) + (1-theta)*dot(grad(d0), w0) ) \
        - (theta*w + (1 -theta)*w0 ), phi) ) * dx

    dis_x = []
    dis_y = []

    from time import sleep

    d_csm1 = File("csm1.pvd")

    if implementation == "1":
        solve(G == 0, wd, bcs, solver_parameters={"newton_solver": \
        {"relative_tolerance": 1E-9,"absolute_tolerance":1E-9,"maximum_iterations":100,"relaxation_parameter":1.0}})
        w,d = wd.split(True)
        w0,d0 = w0d0.split(True)
        d_csm1 << d0

        dis_x.append(d0(coord)[0])
        dis_y.append(d0(coord)[1])

        d_disp.vector()[:] = d.vector()[:] - d0.vector()[:]
        ALE.move(mesh, d_disp)
        mesh.bounding_box_tree().build(mesh)
        w0d0.assign(wd)
        d_csm1 << d0

        w0,d0 = w0d0.split(True)
        dis_x.append(d0(coord)[0])
        dis_y.append(d0(coord)[1])

    print dis_x, dis_y

    plt.figure(1)
    plt.title("Eulerian Mixed, schewed Crank-Nic")
    plt.plot([0, 1],dis_x,); plt.ylabel("Displacement x");plt.xlabel("Time");plt.grid();
    #plt.savefig("run_x.jpg")
    #plt.savefig("./csm1/XDIR_theta%2.f.jpg" % i)
    #plt.show()
    plt.figure(2)
    plt.title("Eulerian Mixed, schewed Crank-Nic")
    plt.plot([0, 1],dis_y);plt.ylabel("Displacement y");plt.xlabel("Time");plt.grid();
    #plt.savefig("./csm1/YDIR_theta%2.f.jpg" % i)
    #plt.savefig("run_y.jpg")
    plt.show()
