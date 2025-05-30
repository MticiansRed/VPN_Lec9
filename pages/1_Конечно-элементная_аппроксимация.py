import streamlit as st
 
menu = st.sidebar.radio('***',
    (
    "Задачи Дирихле", 
    "Аппроксимация расчетной области", 
    "Лагранжевые конечные элементы",
    "Пробные функции",
    "Приближенное решение задачи Дирихле",
    )
)
  
if menu == "Задачи Дирихле":
    r"""
##### Задача Дирихле

Обыкновенное дифференциальное уравнение второго порядка

$\begin{aligned} - 
 \frac{d}{dx} \left (k(x) \frac{du}{dx} \right ) + q(x) u = f(x),
  \quad a < x < b
\end{aligned}$

с переменными коэффициентами

$\begin{aligned}
  k(x) \geq \kappa > 0,
  \quad q(x) \geq 0
\end{aligned}$

Граничные условия Дирихле

$\begin{aligned}
  u(a) = \mu_1,
  \quad u(b) = \mu_2
\end{aligned}$

    """    
    
if menu == "Аппроксимация расчетной области":
    r"""
##### Аппроксимация расчетной области

Расчетная область 

$\begin{aligned}
 \Omega = [a,b]
\end{aligned}$ 

Сетка  

$\begin{aligned}
   \overline{\omega} = \{x~ |~ x = \bar{x}_\alpha , \ \alpha  = 0,1,\ldots,m, \ a = \bar{x}_0  < \bar{x}_1 < \cdots < \bar{x}_m = b \}
\end{aligned}$

* $\omega$ - множество внутренних узлов ($\bar{x}_\alpha, \ \alpha  = 1,2,\ldots,m-1$) 
* $\partial\omega$ - множество граничных узлов ($\bar{x}_0, \bar{x}_m$)

Ячейки

$\begin{aligned}
 \Omega = \bigcup_{\alpha=1}^{m} \Omega_\alpha  
\end{aligned}$
 
$\begin{aligned}
 \Omega_\alpha = \{ x \ | \ \bar{x}_{\alpha -1} \leq x \leq \bar{x}_\alpha\}, \ \alpha = 1, 2, \ldots , m 
\end{aligned}$

    """  
    
if menu == "Лагранжевые конечные элементы":
    r"""
##### Лагранжевые конечные элементы

Область

$\begin{aligned}
 \Omega = \bigcup_{\alpha=1}^{m} \Omega_\alpha  
\end{aligned}$

Ячейки
 
$\begin{aligned}
 \Omega_\alpha = \{ x \ | \ \bar{x}_{\alpha -1} \leq x \leq \bar{x}_\alpha\}, \ \alpha = 1, 2, \ldots , m 
\end{aligned}$

Полином степени $p$ на ячейке

$y(x) = a_0 + a_1 x + \cdots + a_p x^p$

Узлы аппроксимации

* два узла сетки
* $p-1$ внутренних узла (равномерное разбиение ячейки)

    """   
    
if menu == "Пробные функции":
    r"""
##### Пробные функции

Узлы аппроксимации

$\begin{aligned}
 x_i,  \quad i = 0,1, \ldots, n 
\end{aligned}$

Пробные функции

$\begin{aligned}
 \varphi_i(x),
 \quad i = 0,1, \ldots, n 
\end{aligned}$

Основное свойство

$\begin{aligned}
 \varphi_i(x_j) = \left \{ \begin{array}{cc}
  1 ,  &  i = j \\
  0 ,  &  i \neq j \\
\end{array}
\right .
\end{aligned}$

**Визуализация**
    """
        
    import matplotlib.pyplot as plt
    import numpy as np 
    from fenics import *

    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Число ячеек")        
    m = c2.slider("m", 2, 20, 5)
    c1.write("$~~$")
    c1.write("Порядок полинома")        
    p = c2.slider("p", 1, 5, 2)     
    
    a = 0.
    b = 1.
    
    mesh = IntervalMesh(m, a, b)
    xm = mesh.coordinates()
    ym = np.zeros((m+1), "float") 
    
    V = FunctionSpace(mesh, "CG", p)
    n = V.dim()-1

    c1.write("Пробная функция")        
    k = c2.slider("k", 0, n, 0)   
    u = Function(V)
    u.vector()[k] = 1 
    
    xn = np.linspace(a, b, n+1) 
    yn = np.zeros((n+1), "float")   
    
    N = 500
    xx = np.linspace(a, b, N) 
    yy = np.linspace(a, b, N)  
    
    for i in range(0, N): 
        yy[i]  = u(Point(xx[i]))
              
    fig1 = plt.figure(1)
    ss = "$m = $" + str(m) + "$, \ p = $" + str(p)
    plt.title(ss)
    
    plt.scatter(xn, yn)  
    plt.scatter(xm, ym)  
    plt.plot(xx, yy)  
        
    plt.xlabel('$x$') 
    plt.grid(True) 
    
    c1, c2, = st.columns([3,1]) 
    c1.pyplot(fig1)   
    
if menu == "Приближенное решение задачи Дирихле":
    r"""
##### Приближенное решение задачи Дирихле

Краевая задача

$\begin{aligned} - 
 \frac{d}{dx} \left (k(x) \frac{du}{dx} \right ) + q(x) u = f(x),
  \quad a < x < b
\end{aligned}$

$\begin{aligned}
  u(a) = \mu_1,
  \quad u(b) = \mu_2
\end{aligned}$

Приближенное решение

$\begin{aligned}
 u_h(x) \approx u(x),
  \quad a < x < b
\end{aligned}$

Конечно-элементная аппроксимация

$\begin{aligned}
   u_h(x) = \sum_{k=0}^{n} c_k \varphi_k(x) 
\end{aligned}$

    """