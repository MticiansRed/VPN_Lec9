import streamlit as st
 
menu = st.sidebar.radio('***',
    (
    "Тестовая задача", 
    "Ключевые фрагменты кода (FEniCS)", 
    "Влияние разрывного коэффициента", 
    "Влияние шага сетки",
    "Влияние порядка полинома",
    )
)

if menu == "Тестовая задача":
    r"""
##### Тестовая задача

Уравнение

$\begin{aligned} - 
 \frac{d}{dx} \left (k(x) \frac{du}{dx} \right ) = 1,
  \quad 0 < x < 1
\end{aligned}$

с переменным коэффициентом

$\begin{aligned}
 k(x) = \left \{ \begin{array}{cc}
  1 ,  &  0 \leq x < 0.5 \\
  \kappa ,  &  0.5 \leq x \leq 1 \\
\end{array}
\right .
\end{aligned}$

Граничные условия 

$\begin{aligned}
  u(0) = 0,
  \quad u(1) = 0
\end{aligned}$

    """     
if menu == "Ключевые фрагменты кода (FEniCS)":
    r"""
##### Ключевые фрагменты кода (FEniCS)

**Сетка**
    """    
    code = """  
mesh = IntervalMesh(m, 0, 1)
    """ 
    st.code(code, language="python")  
    r"""    
**Конечно-элементное пространство**
    """    
    code = """  
V = FunctionSpace(mesh, "CG", p)
u = TrialFunction(V)
v = TestFunction(V)
    """ 
    st.code(code, language="python")    
    r"""        
**Граничные условия Дирихле**
    """    
    code = """  
def boundary(x, on_boundary):
    return on_boundary
bc = DirichletBC(V, Constant("0."), boundary)
    """ 
    st.code(code, language="python")   
    r"""    
**Вариационная формулировка задачи**
    """    
    code = """  
f = Expression("1", degree=p+2)
kk = Expression("x[0] < 0.5 ? 1 : kap", kap = kap, degree=p+2)
a = kk*dot(grad(u),grad(v))*dx
L = f*v*dx
    """ 
    st.code(code, language="python")    
    r"""        
**Решение задачи**
    """    
    code = """  
w = Function(V)
solve(a == L, w, bc)    
    """ 
    st.code(code, language="python")       
         
                   
if menu == "Влияние разрывного коэффициента":
    r"""
##### Влияние разрывного коэффициента

Коэффициент $\kappa$ увеличивается в 10 и 100 раз

    """

    import matplotlib.pyplot as plt
    import numpy as np 
    from fenics import *
    
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Число ячеек $m$")        
    m = c2.slider("", 5, 50, 10)
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Порядок полинома $p$")        
    p = c2.slider("", 1, 3, 1)  
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Коэффициент $\kappa$")        
    kap = c2.slider("", 0.01, 10., 1., 0.1)     
      
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
    
        w = Function(V)
        solve(a == L, w, bc)
        
        N = 500
        xx = np.linspace(0., 1., N) 
        yy = np.linspace(0., 1., N)  
        
        for i in range(0, N): 
            yy[i]  = w(Point(xx[i]))
        s = "$\\kappa = $" + str(kap)
        plt.plot(xx, yy, label = s) 
        kap = kap*10
    
    ss = "$m = $" + str(m) + "$, \ p = $" + str(p)
    plt.title(ss)
    plt.scatter(xm, ym)  
    plt.xlabel('$x$') 
    plt.legend(loc=0)
    plt.grid(True)      
    
    c1, c2, = st.columns([3,1]) 
    c1.pyplot(fig1)  
 
    
if menu == "Влияние шага сетки":
    r"""
##### Влияние шага сетки

Число ячеек $m$ увеличивается в 2 и 4 раза

    """

    import matplotlib.pyplot as plt
    import numpy as np 
    from fenics import *
    
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Число ячеек $m$")        
    m = c2.slider("", 5, 50, 10)
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Порядок полинома $p$")        
    p = c2.slider("", 1, 3, 1)  
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Коэффициент $\kappa$")        
    kap = c2.slider("", 0.01, 10., 5., 0.1)     
      
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
    
        w = Function(V)
        solve(a == L, w, bc)
        
        N = 500
        xx = np.linspace(0., 1., N) 
        yy = np.linspace(0., 1., N)  
        
        for i in range(0, N): 
            yy[i]  = w(Point(xx[i]))
        s = "$h = $" + str(1./m)
        plt.plot(xx, yy, label = s) 
        m = m*2
    
    ss = "$p = $" + str(p) + "$, \ \\kappa = $" + str(kap)
    plt.title(ss)
    plt.xlabel('$x$') 
    plt.legend(loc=0)
    plt.grid(True)      
    
    c1, c2, = st.columns([3,1]) 
    c1.pyplot(fig1)  
    
if menu == "Влияние порядка полинома":
    r"""
##### Влияние порядка полинома

Порядок лагранжевого полинома $p$ увеличивается в 1 и на 2

    """

    import matplotlib.pyplot as plt
    import numpy as np 
    from fenics import *
    
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Число ячеек $m$")        
    m = c2.slider("", 5, 50, 10)
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Порядок полинома $p$")        
    p = c2.slider("", 1, 3, 1)  
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Коэффициент $\kappa$")        
    kap = c2.slider("", 0.01, 10., 5., 0.1)     
      
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
    
        w = Function(V)
        solve(a == L, w, bc)
        
        N = 500
        xx = np.linspace(0., 1., N) 
        yy = np.linspace(0., 1., N)  
        
        for i in range(0, N): 
            yy[i]  = w(Point(xx[i]))
        s = "$p = $" + str(p)
        plt.plot(xx, yy, label = s) 
        p = p+1
    
    ss = "$m = $" + str(m) + "$, \ \\kappa = $" + str(kap)
    plt.title(ss)
    plt.scatter(xm, ym)  
    plt.xlabel('$x$') 
    plt.legend(loc=0)
    plt.grid(True)      
    
    c1, c2, = st.columns([3,1]) 
    c1.pyplot(fig1)  