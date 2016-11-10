def sigma_f(p, u):
    return - p*Identity(2) + mu*(grad(u) + grad(u).T)

def eps(v):
    return 0.5*(grad(v).T + grad(v))

#MY WAY
def integrateFluidStress(p, u):
  eps   = 0.5*(grad(u) + grad(u).T)
  sig   = -p*Identity(2) + 2.0*mu*eps

  traction  = dot(sig, -n)

  forceX = traction[0]*ds(1)
  forceY = traction[1]*ds(1)
  fX = assemble(forceX)
  fY = assemble(forceY)

  return fX, fY
