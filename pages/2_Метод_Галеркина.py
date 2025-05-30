import streamlit as st
 
menu = st.sidebar.radio('***',
    (
    "Задачи приближения", 
    "Краевая задача", 
    "Вариационная формулировка задачи",
    "1D краевая задача",
    "Конечно-элементное решение",
    "Линейные лагранжевые элементы",    
    )
)
  
if menu == "Задачи приближения":
    r"""
##### Задача приближения

Для заданной $f(x)$

$\begin{aligned} 
 u_h(x) \approx f(x),
  \quad a < x < b
\end{aligned}$

Аппроксимирующая функция

$\begin{aligned}
   u_h(x) = \sum_{k=0}^{n} c_k \varphi_k(x) 
\end{aligned}$

Метод Галеркина

$\begin{aligned}
  (u_h,\varphi_l) = (f,\varphi_l),
  \quad l = 0,1,\ldots,n 
\end{aligned}$

где $(\cdot, \cdot)$ - скалярное произведение в $L_2(a,b)$

$\begin{aligned}
  (u, v) = \int_a^b u(x) v(x) d x
\end{aligned}$

Система уравнений

$\begin{aligned}
  \sum_{k=0}^{n} c_k (\varphi_k,\varphi_l) = (f,\varphi_l),
  \quad l = 0,1,\ldots,n 
\end{aligned}$

    """    
    
if menu == "Краевая задача":
    r"""
##### Краевая задача

1D задача Дирихле

$\begin{aligned} - 
 \frac{d}{dx} \left (k(x) \frac{du}{dx} \right ) + q(x) u = f(x),
  \quad a < x < b
\end{aligned}$

$\begin{aligned}
  u(a) = \mu_1,
  \quad u(b) = \mu_2
\end{aligned}$

Общая запись

$\begin{aligned}
 \mathcal{L} u = f (x),
 \quad  x \in \Omega 
\end{aligned}$

$\begin{aligned}
 \mathcal{B}  u = \psi (x),
 \quad  x \in \partial \Omega 
\end{aligned}$

Конечно-элементная аппроксимация

$\begin{aligned}
   u_h(x) = \sum_{i=1}^{n} c_i \varphi_i(x) 
\end{aligned}$

Проблема

$\begin{aligned}
 \mathcal{L} u_h = \ ?
\end{aligned}$

Пробные функции не имеют нужной гладкости

    """
    
if menu == "Вариационная формулировка задачи":
    r"""
##### Вариационная формулировка задачи

Краевая задача

$\begin{aligned}
 \mathcal{L} u = f (x),
 \quad  x \in \Omega 
\end{aligned}$

$\begin{aligned}
 \mathcal{B}  u = \psi (x),
 \quad  x \in \partial \Omega 
\end{aligned}$

Скалярные произведения в $L_2(\Omega)$ и $L_2(\partial\Omega)$

$\begin{aligned}
  (u, v) = \int_\Omega u(x) v(x) d x
\end{aligned}$

$\begin{aligned}
  (u, v)_{\partial\Omega} = \int_{\partial\Omega} u(x) v(x) d x
\end{aligned}$

Домножаем на достаточно гладкую тестовую функцию 

$\begin{aligned}
 (\mathcal{L} u, v) = (f,v) 
\end{aligned}$

$\begin{aligned}
 (\mathcal{B}  u, v)_{\partial\Omega} = (\psi, v)_{\partial\Omega}
\end{aligned}$

Часть производных в $(\mathcal{L} u, v)$ перебрасываем на тестовую функцию

$\begin{aligned}
 (\mathcal{L} u, v) = a(u,v) + \cdots
\end{aligned}$

Билинейная (относительно $u$ и $v$) форма $a(u,v)$ включает

$\begin{aligned}
 (\cdot u, \cdot v), 
 \quad (\cdot u, \cdot v)_{\partial\Omega}
\end{aligned}$

где $\cdot$ есть некоторая функция, оператор

Линейная форма $l(v)$ не зависит от  $u$ и включает 

$\begin{aligned}
 (\cdot , \cdot v), 
 \quad (\cdot , \cdot v)_{\partial\Omega}
\end{aligned}$

Интегральное тождество (с учетом уравнения и граничных условий)

$\begin{aligned}
 a (u,v) = l(v)
\end{aligned}$

Вариационная задача: найти $u$, такую, что

$\begin{aligned}
 a (u,v) = l(v),
 \quad \forall v
\end{aligned}$

    """  
    
if menu == "1D краевая задача":
    r"""
##### 1D краевая задача

Задача Дирихле

$\begin{aligned} - 
 \frac{d}{dx} \left (k(x) \frac{du}{dx} \right ) + q(x) u = f(x),
  \quad a < x < b
\end{aligned}$

$\begin{aligned}
  u(a) = \mu_1,
  \quad u(b) = \mu_2
\end{aligned}$

Пространства

* $V$ - пространство достаточно гладких функций

* Подпространство $V_D = \{ v \ | \ v \in V, \quad v(a) = \mu_1,
  \quad u(b) = \mu_2 \}$ 
  
* Подпространство $V_0 = \{ v \ | \ v \in V, \quad v(a) = 0,
  \quad u(b) = 0 \}$  
  
Для $u \in V_D$ и $v \in V_0$  
  
$\begin{aligned}
 (\mathcal{L} u, v) &= - \int_{a}^{b} \frac{d}{dx} \left (k(x) \frac{du}{dx} \right ) v(x) d x + \int_{a}^{b}  q(x) u(x) v(x) dx = \\
 & = - k(x) \frac{du}{dx} {\color {red}  v(x)} \Big |_a^b +  \int_{a}^{b} k(x) \frac{du}{dx} \frac{dv}{dx} d x + \int_{a}^{b}  q(x) u(x) v(x) dx
\end{aligned}$

Билинейная форма 

$\begin{aligned}
 a(u,v) = \int_{a}^{b} k(x) \frac{du}{dx} \frac{dv}{dx} d x + \int_{a}^{b}  q(x) u(x) v(x) dx
\end{aligned}$

Линейная форма 

$\begin{aligned}
 l(v) = \int_{a}^{b} f(x) v(x) d x
\end{aligned}$

Вариационная задача: найти  $u \in V_D$ такую, что

$\begin{aligned}
 a(u,v) = l(v) ,
 \quad \forall v \in V_0
\end{aligned}$

    """   
    
if menu == "Конечно-элементное решение":
    r"""
##### Конечно-элементное решение

$V$ - конечно-элементное пространство: для $y(x) \in V$

$\begin{aligned}
   y(x) = \sum_{k=0}^{n} c_k \varphi_k(x) 
\end{aligned}$

Приближенное решение краевой задачи  $u_h(x) \in V_D$

$\begin{aligned}
   u_h(x) = {\color{red} \mu_1 \varphi_0(x)} + \sum_{k=1}^{n-1} c_k \varphi_k(x) + {\color{red} \mu_2 \varphi_m(x)} = 
 {\color{red} \mu(x)}  + y(x) 
\end{aligned}$

Вариационная задача: найти  $y \in V_0$ такую, что

$\begin{aligned}
 a(y,v) = l(v) - a(\mu,v),
 \quad \forall v \in V_0
\end{aligned}$

Метод Галеркина

$\begin{aligned}
 \sum_{k=1}^{n-1} a(\varphi_k, \varphi_l) c_k = l(\varphi_l) - z(\mu, \varphi_l),
 \quad i = 1,2, \ldots , n-1
\end{aligned}$

    """

if menu == "Линейные лагранжевые элементы":
    r"""
##### Линейные лагранжевые элементы

Для $y_i = u_h(x_i)$ при $x_i - x_{i -1} = h $ 

$\begin{aligned} -
 \frac{1}{h} \Big (a_{i+1} \frac{y_{i+1} - y_i}{h}  - a_{i} \frac{y_{i} - y_{i-1}}{h} \Big ) + d_i y_i = \phi_i
\end{aligned}$

$\begin{aligned}
  y_0 = \mu_1,
  \quad y_m = \mu_2 
\end{aligned}$

Коэффициенты

$\begin{aligned}
  a_i = \frac 1h \int\limits_{x_{i-1}}^{x_{i}} k(x) dx -
  \frac 1h \int\limits_{x_{i-1}}^{x_{i}} q(x) (x-x_{i-1})(x_{i}-x)   dx
\end{aligned}$

$\begin{aligned}
  d_i = \frac{1}{h^2} \left (
   \int\limits_{x_{i-1}}^{x_{i}} q(x)(x-x_{i-1}) dx +
   \int\limits_{x_{i}}^{x_{i+1}} q(x)(x_{i+1}-x) dx
   \right )
\end{aligned}$

$\begin{aligned}
  \phi_i = \frac{1}{h^2} \left (
   \int\limits_{x_{i-1}}^{x_{i}} f(x)(x-x_{i-1}) dx +
   \int\limits_{x_{i}}^{x_{i+1}} f(x)(x_{i+1}-x) dx
   \right )
\end{aligned}$

    """
    
