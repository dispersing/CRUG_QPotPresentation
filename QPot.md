QPOT: Quasi-POTential analysis for stochastic differential equations
========================================================
author: Christopher M. Moore
date: 17 August 2016
autosize: true

Cleveland R Users Group

Modeling nature
========================================================
*The physicist Heinz Pagels (1982) once observed that there seem to be two kinds of people in the world.*

1. **There are those who seek
and find deterministic order and meaning, if not purpose,
in every event.** 
2. **And then there are those who believe
events to be influenced, if not dominated, by intrinsically
inscrutable, and meaningless, random chance.**

<!-- *One of the intellectual triumphs of twentieth-century physics was to prove that both views of physical nature are simultaneously true and correct, but on very different spatial and temporal scales.* -->
---S. Hubbell (2001), *The Unified Neutral Theory of Biodiversity and Biogeography* (Princeton, p. 7)

<img src="QPot-figure/unnamed-chunk-1-1.png" title="plot of chunk unnamed-chunk-1" alt="plot of chunk unnamed-chunk-1" width="750px" height="500px" style="display: block; margin: auto;" />



Ecology: a science of dynamics and complexity
========================================================
**The study of the spatial and temporal patterns of the distributions and abundances of organisms, including causes and consequences**

Ecology is hierarchically organized:
  - Individuals
  - Populations
  - Communities
  - Ecosystems

Today's focus, *population ecology*, is strongly quantitative and mathematical

Some R notes, pt. I
========================================================
- Packages I'll be using if you want to follow along

```r
  library(package = "deSolve") # Solver for differential equations
  library(package = "rootSolve") # Finds roots, Jacobian matrices, and estimates steady-state conditions
  library(package = "phaseR") # Tools for graphing 
  library(package = "MASS") # For 2-D density estimator
  library(package = "QPot") # Our package!
  library(package = "viridis") # Great color map
```
- For this presentation, I opened RStudio for the second time, and updated from v. 0.98.447 (current is v. 0.99.902)

### Based primarily on two works
  1. Nolting and Abbott's *Ecology* paper
    - B. C. Nolting and K. C. Abbott. Balls, cups, and quasi-potentials: Quantifying stability in stochastic systems. Ecology, 97(4):850–864, 2016. [Link](http://onlinelibrary.wiley.com/doi/10.1890/15-1047.1/abstract)
  2. Our R package and paper
    - C.M. Moore, C.R. Stieha, B.C. Nolting, M.K. Cameron, and K.C. Abbott. QPot: Quasi-Potential Analysis for Stochastic Differential Equations, 2016. URL [https://www.R-project.org/package=QPot](https://www.R-project.org/package=QPot). R package version 1.2.
    - Moore, C.M., Stieha, C.R., Nolting, B.C., Cameron, M.K. and Abbott, K.C. 2015. QPot: An R Package for Stochastic Differential Equation Quasi-Potential Analysis. [arXiv preprint arXiv:1510.07992](http://arxiv.org/abs/1510.07992) or the most updated version (accapted and final version): [link](http://mutualismecology.com/publications/preprints/Moore_et_al._2016.pdf).

Differential equations: as models
========================================================
$$\frac{\text{change in a variable}}{\text{change in another variable}}$$

if $~x$ is a variable that is dependent on variable $~y$, then we can write

$$\frac{d x}{d y}$$

with time, $~t$, as the independent variable

$$\frac{d x}{d t}$$

Differential equations: some examples, pt. I
========================================================
Growth of $~x$ at a rate, $~\alpha$
$$\frac{dx}{dt} = \alpha x $$

Logistic equation, with growth of $~x$ at a rate $~\alpha$, and negative density-dependent growth (crowding effect) of strength $\gamma$ and interaction between $x \times x$:

$$\frac{dx}{dt} = \alpha x - \gamma x^2$$


Differential equations: some examples, pt. II
========================================================
Competition between species $~N_1$ and $~N_2$

$$
  \frac{dN_1}{dt} = r_1N_1 - \alpha_{11} N_1N_1  - \beta_{21} N_2N_1 \\
  \frac{dN_2}{dt} = r_2N_2 - \alpha_{22} N_2N_2  - \beta_{12} N_1N_2
$$

Disease model of susceptible, $~S$, infected $~I$, and recovered, $~R$

$$
\begin{aligned}
  \frac{dS}{dt} &= \alpha S - \beta SI + \delta R - d_SS \\
  \frac{dI}{dt} &= \beta SI -\gamma I - d_II \\
  \frac{dR}{dt} &= \gamma I -\delta R - d_RR
\end{aligned}
$$


Differential equations in R: deSolve format, pt. I
========================================================
`deSolve` has set the standard, with:

- Function arguments
  - time vector (e.g., `seq(from = 0, to = 10, by = 0.1`)
  - initial conditions as a named atomic vector, `y0`, for each of the variables; e.g., `c(a = 2, b = 3.5)`
  - parameters as a named atomic vector; e.g., `c(alpha = 1, beta = 2)`
- Returns a list of variable values

Differential equations in R: deSolve format, pt. II
========================================================
*The mathematics of uncontrolled growth are frightening.*  **A single cell** *of the bacterium E. coli would, under ideal circumstances,* **divide every twenty minutes**.  *That is not particularly disturbing until you think about it, but the fact is that bacteria multiply geometrically: one becomes two, two become four, four become eight, and so on.  In this way,* **it can be shown that in a single day, one cell of E. coli could produce a super-colony equal in size and weight to the entire planet earth.**---Michael Crichton (1969), *The Andromeda Strain* (Dell, New York, p. 247)

$$\frac{dx}{dt} = \alpha x = (1)x$$

Differential equations in R: deSolve format, pt. III
========================================================
$$\frac{dx}{dt} = \alpha x = (1)x$$

```r
growth <- function(t, y0, parameters){
  with(as.list(c(y0, parameters)), {
    dx = alpha*x
    list(c(dx))
  })
}
```


```r
growth(y0 = c(x = 1), parameters = c(alpha = 1))
```

```
[[1]]
[1] 1
```

Differential equations in R: deSolve::ode
========================================================
$$\frac{dx}{dt} = \alpha x = (1)x$$


```r
start <- c(x = 1)
parameters <- c(alpha = 1)
minutes <- seq(from = 0, to = 3*24, by = 1)

out <- ode(y = start, times = minutes, func = growth, parms = parameters, method = "euler")
out[1:3,]
```

```
     time x
[1,]    0 1
[2,]    1 2
[3,]    2 4
```

Differential equations in R: plot deSolve::ode
========================================================
&nbsp;

$$\frac{dx}{dt} = \alpha x = (1)x$$


```
        time            x 
7.200000e+01 4.722366e+21 
```

```r
plot(out[,1], out[,2], xlab = "Time", ylab = "x", type = "l", lwd = 2)
```

<img src="QPot-figure/unnamed-chunk-6-1.png" title="plot of chunk unnamed-chunk-6" alt="plot of chunk unnamed-chunk-6" width="480px" height="480px" />

***
&nbsp;  
&nbsp;  
&nbsp;  

### Was Michael Crichton correct?

```
[1] "Mass of an E. coli cell = 1e-12 (g)"
```

```
[1] "times 4.7e+21 = 4.7e+09 (g)"
```

```
[1] "The mass of Earth is = 6e+27 (g)"
```

```
[1] "MC was incorrect about mass"
```

```
[1] "Volume of an E. coli cell = 1e-18 (mm^3)"
```

```
[1] "times 4.7e+21 = 4.7e+09 (mm^3)"
```

```
[1] "The volume of Earth is = 1.1e+21 (mm^3)"
```

```
[1] "MC was incorrect about volume"
```

Differential equations in R: 2-D competition
========================================================
Competition between species $~N_1$ and $~N_2$

$$\frac{dN_1}{dt} = r_1N_1 - \alpha_{11} N_1N_1  - \beta_{21} N_2N_1 \\
\frac{dN_2}{dt} = r_2N_2 - \alpha_{22} N_2N_2  - \beta_{12} N_1N_2$$

- $r_i$, the population growth rate of species $i$
- $\alpha_{ii}$, competition between the same species $i$
- $\beta_{ji}$, the competitive effect of species $j$ on species $i$

Differential equations in R: competition in ode
========================================================

```r
comp <- function(t, y0, parameters){
  with(as.list(c(y0, parameters)), {
  dx <- r1*x - a11*x*x - b21*y*x
  dy <- r2*y - a22*y*y - b12*x*y
  list(c(dx, dy))
  })
}
state <- c(x = 1, y = 4)
parameters <- c(r1 = 1, r2 = 1, a11 = 0.1, a22 = 0.1, b12 = 0.05, b21 = 0.05)
time <- seq(from = 0, to = 25, length.out = 100)
out <- ode(y = state, times = time, func = comp, parms = parameters)
```

Differential equations in R: competition over time
========================================================

```r
plot(out[,1], out[,2], type = "l", col = "red", xlab = "Time", ylab = "Density of\nspecies 1 (red) and species 2 (blue)", ylim = c(0, max(out[,c(2,3)], na.rm = T)), lwd = 2)
lines(out[,1], out[,3], col = "blue", lwd = 2)
```

<img src="QPot-figure/unnamed-chunk-9-1.png" title="plot of chunk unnamed-chunk-9" alt="plot of chunk unnamed-chunk-9" width="1000px" height="666px" />

Differential equations in R: competition in phase space
========================================================

```r
plot(out[,2], out[,3], xlab = "Density of species 1", ylab = "Density of species 2", type = "l", lwd = 2, xlim = c(0, 10), ylim = c(0, 10))
points(out[1, 2], out[1, 3], pch = 16)
```

<img src="QPot-figure/unnamed-chunk-10-1.png" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" width="600px" height="600px" />

Differential equations in R: competition trajectories
========================================================

```r
plot(out[,2], out[,3], xlab = "Density of species 1", ylab = "Density of species 2", type = "l", lwd = 3, xlim = c(0, 10), ylim = c(0, 10))
for (i in 1:25){
  state <- c(x = runif(n = 1, min = 0, max = 10), y = runif(n = 1, min = 0, max = 10))
  out <- ode(y = state, times = time, func = comp, parms = parameters)
  lines(out[,2], out[,3], lwd = 2)
  points(out[1, 2], out[1, 3], pch = 16)
}
```

<img src="QPot-figure/unnamed-chunk-11-1.png" title="plot of chunk unnamed-chunk-11" alt="plot of chunk unnamed-chunk-11" width="600px" height="600px" />

Differential equations in R: phaseR::flowField format
========================================================

```r
comp <- function(t, y0, parameters){
  with(as.list(c(y0, parameters)), {
  x <- y0[1]
  y <- y0[2]
  dx <- r1*x - a11*x*x - b21*y*x
  dy <- r2*y - a22*y*y - b12*x*y
  list(c(dx, dy))
  })
}
```

Differential equations in R: phaseR::flowField plot
========================================================

```r
ff <- flowField(deriv = comp, x.lim = c(0, 10), y.lim = c(0, 10), parameters = parameters, points = 20, add = F, col = "black", arrow.type = "equal")
```

<img src="QPot-figure/unnamed-chunk-13-1.png" title="plot of chunk unnamed-chunk-13" alt="plot of chunk unnamed-chunk-13" width="800px" height="800px" />
<img src="QPot-figure/unnamed-chunk-14-1.png" title="plot of chunk unnamed-chunk-14" alt="plot of chunk unnamed-chunk-14" width="800px" height="800px" />

Differential equations in R: phaseR::flowField trajectories
========================================================

```r
ff <- flowField(deriv = comp, x.lim = c(0, 10), y.lim = c(0, 10), parameters = parameters, points = 20, add = F, col = "grey50")
for (i in 1:25){
  xy.start <- c(runif(n = 1, min = 0, max = 10), runif(n = 1, min = 0, max = 10))
  traj <- trajectory(deriv = comp, y0 = xy.start, parameters = parameters, col = "black", t.end = 10, pch = 16, lwd = 2)
}
```

<img src="QPot-figure/unnamed-chunk-15-1.png" title="plot of chunk unnamed-chunk-15" alt="plot of chunk unnamed-chunk-15" width="600px" height="600px" />

Differential equations in R: competitive exclusion
========================================================

```r
parameters.ce <- c(r1 = 1, r2 = 1, a11 = 0.1, a22 = 0.1, b21 = 0.25, b12 = 0.25)
ff <- flowField(deriv = comp, x.lim = c(0, 10), y.lim = c(0, 10), parameters = parameters.ce, points = 20, add = F, col = "grey50")
abline(h = 0, v = 0)
for (i in 1:15){
  xy.start <- c(runif(n = 1, min = 0, max = 10), runif(n = 1, min = 0, max = 10))
  traj <- trajectory(deriv = comp, y0 = xy.start, parameters = parameters.ce, col = "black", t.end = 25, lwd = 2, pch = 16)
}
```

<img src="QPot-figure/CE-1.png" title="plot of chunk CE" alt="plot of chunk CE" width="600px" height="600px" />


Differential equations in R: cyclic behaviors
========================================================

$$\frac{dR}{dt} = bR - cRC \\
\frac{dC}{dt} = dRC - eC$$

```r
ff <- flowField(deriv = prey.predator, x.lim = c(0, 3), y.lim = c(0, 3), parameters = parameters, points = 20, add = F, col = "black", xlab = "R", ylab = "C")
```

<img src="QPot-figure/unnamed-chunk-17-1.png" title="plot of chunk unnamed-chunk-17" alt="plot of chunk unnamed-chunk-17" width="600px" height="600px" />
<img src="QPot-figure/unnamed-chunk-18-1.png" title="plot of chunk unnamed-chunk-18" alt="plot of chunk unnamed-chunk-18" width="600px" height="600px" /><img src="QPot-figure/unnamed-chunk-18-2.png" title="plot of chunk unnamed-chunk-18" alt="plot of chunk unnamed-chunk-18" width="600px" height="600px" />

Differential equations in R: unstable behaviors, pt. I
========================================================

$$\frac{dN_1}{dt} = r_1N_1 - \alpha _{11}N_1^2 + \beta _{21}N_2N_1 \\
\frac{dN_2}{dt} = r_2N_2 - \alpha _{22}N_2^2 + \beta _{12}N_1N_2$$

```r
ff <- flowField(deriv = mut, x.lim = c(0, 2), y.lim = c(0, 2), parameters = parameters, points = 20, add = F, col = "black")
nc <- nullclines(deriv = mut, x.lim = c(0, 2), y.lim = c(0, 2), parameters = parameters, lwd = 2)
```

<img src="QPot-figure/unnamed-chunk-20-1.png" title="plot of chunk unnamed-chunk-20" alt="plot of chunk unnamed-chunk-20" width="600px" height="600px" />
<img src="QPot-figure/unnamed-chunk-21-1.png" title="plot of chunk unnamed-chunk-21" alt="plot of chunk unnamed-chunk-21" width="600px" height="600px" />

Differential equations in R: unstable behaviors, pt. II
========================================================

$$\frac{dN_1}{dt} = r_1N_1 - \alpha _{11}N_1^2 + \beta _{21}N_2N_1 \\
\frac{dN_2}{dt} = r_2N_2 - \alpha _{22}N_2^2 + \beta _{12}N_1N_2$$

```r
ff <- flowField(deriv = mut, x.lim = c(0, 2), y.lim = c(0, 2), parameters = parameters, points = 20, add = F, col = "black")
nc <- nullclines(deriv = mut, x.lim = c(0, 2), y.lim = c(0, 2), parameters = parameters, lwd = 2)
```

<img src="QPot-figure/unnamed-chunk-23-1.png" title="plot of chunk unnamed-chunk-23" alt="plot of chunk unnamed-chunk-23" width="600px" height="600px" />
<img src="QPot-figure/unnamed-chunk-24-1.png" title="plot of chunk unnamed-chunk-24" alt="plot of chunk unnamed-chunk-24" width="600px" height="600px" />

Stochasticity: its role in ecology and stability
========================================================
- Model unexplained variation in a variable
- Ecology, like other complex systems, is riddled with a multiplicity of processes ignored by simple, deterministic models
- Consequently, we can model unexplained variation
- [CRAN Task View: Differential Equations](https://cran.r-project.org/web/views/DifferentialEquations.html)

QPOT: stochastic differential equation and a realization
========================================================
$$dN_1 = \left(r_1N_1 - \alpha_{11} N_1N_1  - \beta_{21} N_2N_1\right)dt + \sigma_1 dW_1\\
dN_2 = \left(r_2N_2 - \alpha_{22} N_2N_2  - \beta_{12} N_1N_2\right)dt + \sigma_2 dW_2$$


```r
parameters <- c(r1 = 1, r2 = 1, a11 = 0.1, a22 = 0.1, b12 = 0.05, b21 = 0.05)
model.deltat <- 0.1
ms <- Model2String(model = comp, parms = parameters, deSolve.form = T, supress.print = T)
ts <- TSTraj(y0 = c(x = 1, y = 1), time = 100, deltat = model.deltat, x.rhs = ms[1], y.rhs = ms[2], sigma = 1, lower.bound = 0)
```

QPOT: stochastic simulations, time series
========================================================


```r
TSPlot(mat = ts, deltat = model.deltat)
```

<img src="QPot-figure/unnamed-chunk-26-1.png" title="plot of chunk unnamed-chunk-26" alt="plot of chunk unnamed-chunk-26" width="1500px" height="750px" />

QPOT: stochastic simulations, realization and density
========================================================

```r
TSPlot(ts, deltat = model.deltat, dim = 2)
TSDensity(mat = ts, deltat = model.deltat, dim = 2, col2d = viridis(100, option = "A"))
```

<img src="QPot-figure/unnamed-chunk-28-1.png" title="plot of chunk unnamed-chunk-28" alt="plot of chunk unnamed-chunk-28" width="800px" height="800px" /><img src="QPot-figure/unnamed-chunk-28-2.png" title="plot of chunk unnamed-chunk-28" alt="plot of chunk unnamed-chunk-28" width="800px" height="800px" />

QPOT: what is the quasi-potential
========================================================
The potential relates to the work needed to move from one point to another, with the areas of least potential at surface minima.

<img src="QPot-figure/unnamed-chunk-29-1.png" title="plot of chunk unnamed-chunk-29" alt="plot of chunk unnamed-chunk-29" width="800px" height="500px" />

**The quasi-potential is a tool that yields information about properties of stochastic systems, such as the expected time to escape a basin of attraction, the expected frequency of transitions between basins, and the stationary probability distribution.**


QPOT: how is the quasi-potential calculated
========================================================
### Ordered Upwind Method
**Calucualtes the minimal cost path on a surface**

Sources for more information

1. Cameron, M.K. (2012). Finding the quasipotential for nongradient SDE’s. *Physica D* 241, 1532-1550.
2. Sethian, J.A. and Vladimirsky, A. (2003). Ordered upwind methods for static Hamilton–Jacobi equations: Theory and algorithms. *SIAM Journal on Numerical Analysis* 41(1), 325-363.
3. Sethian, J.A. and Vladimirsky, A. (2001). Ordered upwind methods for static Hamilton–Jacobi equations. *PNAS* 98(20), 11069-11074.


QPOT: competition example, revisited
========================================================
<img src="QPot-figure/CE-1.png" width="800px" height="800px" />


QPOT: competition example, with stochasticity
========================================================
<img src="QPot-figure/unnamed-chunk-30-1.png" title="plot of chunk unnamed-chunk-30" alt="plot of chunk unnamed-chunk-30" width="500px" height="500px" /><img src="QPot-figure/unnamed-chunk-30-2.png" title="plot of chunk unnamed-chunk-30" alt="plot of chunk unnamed-chunk-30" width="500px" height="500px" /><img src="QPot-figure/unnamed-chunk-30-3.png" title="plot of chunk unnamed-chunk-30" alt="plot of chunk unnamed-chunk-30" width="500px" height="500px" /><img src="QPot-figure/unnamed-chunk-30-4.png" title="plot of chunk unnamed-chunk-30" alt="plot of chunk unnamed-chunk-30" width="500px" height="500px" /><img src="QPot-figure/unnamed-chunk-30-5.png" title="plot of chunk unnamed-chunk-30" alt="plot of chunk unnamed-chunk-30" width="500px" height="500px" /><img src="QPot-figure/unnamed-chunk-30-6.png" title="plot of chunk unnamed-chunk-30" alt="plot of chunk unnamed-chunk-30" width="500px" height="500px" /><img src="QPot-figure/unnamed-chunk-30-7.png" title="plot of chunk unnamed-chunk-30" alt="plot of chunk unnamed-chunk-30" width="500px" height="500px" /><img src="QPot-figure/unnamed-chunk-30-8.png" title="plot of chunk unnamed-chunk-30" alt="plot of chunk unnamed-chunk-30" width="500px" height="500px" />


QPOT: competition example, finding stable equilibrium points
========================================================
## Analytically
 $$
 \begin{aligned}
  \frac{dN_1}{dt} &= r_1N_1 - \alpha_{11} N_1N_1  - \beta_{21} N_2N_1 \\
  0 &= r_1N_1 - \alpha_{11} N_1N_1  - \beta_{21} N_2N_1 \\
  0 &= N_1\left(r_1 - \alpha_{11} N_1  - \beta_{21} N_2\right) \\
  0 &= r_1 - \alpha_{11} N_1  - \beta_{21} N_2 \\
  \alpha_{11} N_1 &= r_1 - \beta_{21} N_2 \\
  N_1^* &= \frac{r_1 - \beta_{21} N_2}{\alpha_{11}}
  \end{aligned}
$$

## Numerically
- deSolve::ode
- rootSolve::steady

QPOT: competition example, local quasi-potential
========================================================

```r
ce <- Model2String(model = comp, parms = parameters.ce, deSolve.form = T)
bnd.x <- c(-5, 15); bnd.y <- c(-5, 15)
local.1 <- QPotential(x.rhs = ce[1], x.start = 0, x.bound = bnd.x, x.num.steps = 2000, y.rhs = ce[2], y.start = 10, y.bound = c(-5, 15), y.num.steps = 2000)
QPContour(surface = local.1, dens = c(500, 500), x.bound = bnd.x, y.bound = bnd.y)
```



<img src="QPot-figure/unnamed-chunk-33-1.png" title="plot of chunk unnamed-chunk-33" alt="plot of chunk unnamed-chunk-33" width="600px" height="600px" /><img src="QPot-figure/unnamed-chunk-33-2.png" title="plot of chunk unnamed-chunk-33" alt="plot of chunk unnamed-chunk-33" width="600px" height="600px" />

QPOT: competition example, global quasi-potentials
========================================================

```r
saddle <- steady(y = c(3,3), func = comp, parms = parameters.ce)
```

```
[1] 2.857143 2.857143
```


```r
global.qp <- QPGlobal(local.surfaces = list(local.1, local.2), unstable.eq.x = c(saddle$y[1]), unstable.eq.y = c(saddle$y[2]), x.bound = bnd.x, y.bound = bnd.y)
QPContour(surface = global.qp, dens = c(1000, 1000), x.bound = bnd.x, y.bound = bnd.y, xlim = c(-1, 11), ylim = c(-1, 11))
```

<img src="QPot-figure/unnamed-chunk-35-1.png" title="plot of chunk unnamed-chunk-35" alt="plot of chunk unnamed-chunk-35" width="800px" height="800px" />

QPOT: direction, gradient, and remainder fields
========================================================


```r
VD <- VecDecomAll(surface = global.qp, x.rhs = ce[1], y.rhs = ce[2], x.bound = bnd.x, y.bound = bnd.y)
VecDecomPlot(x.field = VD[,,1], y.field = VD[,,2], dens = c(20, 20), x.bound = bnd.x, y.bound = bnd.y)
VecDecomPlot(x.field = VD[,,3], y.field = VD[,,4], dens = c(20, 20), x.bound = bnd.x, y.bound = bnd.y)
VecDecomPlot(x.field = VD[,,5], y.field = VD[,,6], dens = c(20, 20), x.bound = bnd.x, y.bound = bnd.y)
```

<img src="QPot-figure/unnamed-chunk-37-1.png" title="plot of chunk unnamed-chunk-37" alt="plot of chunk unnamed-chunk-37" width="650px" height="650px" /><img src="QPot-figure/unnamed-chunk-37-2.png" title="plot of chunk unnamed-chunk-37" alt="plot of chunk unnamed-chunk-37" width="650px" height="650px" /><img src="QPot-figure/unnamed-chunk-37-3.png" title="plot of chunk unnamed-chunk-37" alt="plot of chunk unnamed-chunk-37" width="650px" height="650px" />

QPOT: a disagreement in stability
========================================================
 $$
 \begin{aligned}
  \frac{dx(t)}{dt} &= \alpha x(t)\left(1 - \frac{x(t)}{\beta}\right) - \frac{\delta x^2(t)y(t)}{\kappa + x^2(t)} \\
  \frac{dy(t)}{dt} &= \frac{\gamma x^2(t)y(t)}{\kappa + x^2(t)} - \mu y^2(t)
  \end{aligned}
$$

<img src="QPot-figure/unnamed-chunk-39-1.png" title="plot of chunk unnamed-chunk-39" alt="plot of chunk unnamed-chunk-39" width="700px" height="700px" />
<img src="QPot-figure/unnamed-chunk-40-1.png" title="plot of chunk unnamed-chunk-40" alt="plot of chunk unnamed-chunk-40" width="500px" height="313px" />
<img src="QPot-figure/unnamed-chunk-41-1.png" title="plot of chunk unnamed-chunk-41" alt="plot of chunk unnamed-chunk-41" width="700px" height="700px" />


Thank you
========================================================
- Website: [mutualismecology.com](http://www.mutualismecology.com)
- Email: [life.dispersing@gmail.com](life.dispersing@gmail.com)
- Twitter: [@lifedispersing](https://twitter.com/lifedispersing)
