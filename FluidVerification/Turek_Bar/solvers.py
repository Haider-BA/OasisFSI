from fenics import NonlinearVariationalProblem, NonlinearVariationalSolver, assemble, solve, \
norm, MPI, mpi_comm_world

def Newtonsolver(F, up, bcs, J, atol, rtol, max_it, lmbda):

    problem = NonlinearVariationalProblem(F, up, bcs, J)
    sol  = NonlinearVariationalSolver(problem)

    prm = sol.parameters
    #info(prm,True)  #get full info on the parameters
    prm['nonlinear_solver'] = 'newton'
    prm['newton_solver']['absolute_tolerance'] = atol
    prm['newton_solver']['relative_tolerance'] = rtol
    prm['newton_solver']['maximum_iterations'] = max_it
    prm['newton_solver']['relaxation_parameter'] = lmbda
    prm['newton_solver']['linear_solver'] = 'mumps'

    sol.solve()

def Newton_manual(F, up, bcs, J, atol, rtol, max_it, lmbda\
                 , up_res):
    #Reset counters
    Iter      = 0
    residual   = 1
    rel_res    = residual
    while rel_res > rtol and residual > atol and Iter < max_it:

        A = assemble(J)
        b = assemble(-F)

        [bc.apply(A, b, up.vector()) for bc in bcs]

        solve(A, up_res.vector(), b)

        up.vector().axpy(1., up_res.vector())
        [bc.apply(up.vector()) for bc in bcs]
        rel_res = norm(up_res, 'l2')
        residual = b.norm('l2')

        if MPI.rank(mpi_comm_world()) == 0:
            print "Newton iteration %d: r (atol) = %.3e (tol = %.3e), r (rel) = %.3e (tol = %.3e) " \
        % (Iter, residual, atol, rel_res, rtol)
        Iter += 1

    return up
