import matplotlib.pyplot as plt
import numpy as np 
from fenics import *

m = 10
p = 1
kap = 5.

fig1 = plt.figure(1)

for kk in range(0, 3): 

    mesh = IntervalMesh(m, 0, 1)
    xm = mesh.coordinates()
    ym = np.zeros((m+1), "float") 
    
    V = FunctionSpace(mesh, "CG", p)
    n = V.dim()-1
    
    u = TrialFunction(V)
    v = TestFunction(V)
    
    def boundary(x, on_boundary):
        return on_boundary
    bc = DirichletBC(V, Constant("0."), boundary)
    
    f = Expression("1", degree=p+2)
    kk = Expression("x[0] < 0.5 ? 1 : kap", kap = kap, degree=p+2)
    a = kk*dot(grad(u), grad(v))*dx
    L = f*v*dx

    u = Function(V)
    solve(a == L, u, bc)
    
    N = 500
    xx = np.linspace(0., 1., N) 
    yy = np.linspace(0., 1., N)  
    
    for i in range(0, N): 
        yy[i]  = u(Point(xx[i]))
    s = "$p = $" + str(p)
    plt.plot(xx, yy, label = s) 
    p = p+1

ss = "$m = $" + str(m) + "$, \ \\kappa = $" + str(kap)
plt.title(ss)
plt.scatter(xm, ym)  
plt.xlabel('$x$') 
plt.legend(loc=0)
plt.grid(True)     

plt.show()
