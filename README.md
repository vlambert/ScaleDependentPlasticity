# Elastoplastic contact with scale-dependent plastic strength
This directory includes a numerical program to solve for the distribution of elastic, plastic and non-contact area along a randomly rough surface at varying scales including a scale-dependent plastic yield stress following the theoretical developments of Persson ([2001](https://doi.org/10.1016/S0039-6028(98)00051-X);[2006](https://doi.org/10.1016/j.surfrep.2006.04.001)).

## Problem Statement
The problem statement considers a rough surface with linear scale $$L$$ that can be examined at different magnifications $$\zeta = L/\lambda$$, where $$\lambda$$ is the shortest wavelength of roughness which is resolved at magnification $$\zeta$$. 

$$ \frac{\partial P(\sigma,\zeta)}{\partial \zeta} = f(\zeta) \frac{\partial^2 P(\sigma,\zeta)}{\partial \sigma^2} $$ 

## Numerical Solution
To solve the governing equation subject to the evolving yield stress boundary condition, we implement an implicit Crank-Nicholson finite difference scheme using a modified Thomas Algorithm.

The interior solution of the diffusion problem is solved wiht a standard finite difference scheme, with the solutions to the governing equation being given at stress increments $j$ and magnification step $i+1$ in the form:

$$    a_{j}^{i+1}P_{j-1}^{i+1}+b_{j}^{i+1}P_{j}^{i+1} + c_{j}^{i+1}P_{j+1}^{i+1} = d_{j}^{i} $$

Following an implicit Crank-Nicolson scheme, the coefficients take the form:

$$    a_{j}^{i+1} = -\frac{\Delta\zeta}{2}\frac{F^{i+1}}{\Delta\sigma^2} $$

$$    b_{j}^{i+1} = 1+ 2\frac{\Delta\zeta}{2}\frac{F^{i+1}}{\Delta\sigma^2} $$

$$    c_{j}^{i+1} =-\frac{\Delta\zeta}{2}\frac{F^{i+1}}{\Delta\sigma^2} $$

$$    d_{j}=\frac{\Delta\zeta}{2}\frac{F^{i}}{\Delta\sigma^2}\bigg( P_{j-1}^{i} + P_{j+1}^{i}\bigg) + \bigg(1-2\frac{\Delta\zeta}{2}\frac{F^{i}}{\Delta\sigma^2}\bigg) P_{j}^{i} $$

which can all be evaluated based on information from the current magnification step $\zeta_i$ and knowledge of the diffusivity $F^{i} = f(\zeta_i)$, which is precribed in our calculations. The system of equations form matrix vector product with a tridiagonal matrix and can be solved in parallel using a modified Thomas Algorithm. This implicit Crank-Nicolson scheme is numerically stable and 2nd-order accurate in space and time within the domain $[0,Y(\zeta))$.

### Discretization scheme
We construct a 1-D grid in stress space spanning the domain $$[0,Y(1)]$$ discretized by N grid points with uniform spacing $\Delta\sigma = Y(1)/(N-1)$. We choose magnification increments of the form:

$$    \Delta\zeta = \alpha \frac{\Delta \sigma^2}{f(\zeta)}  $$

to satisfy the Courant conditon. Higher diffusivities $f(\zeta)$ require smaller magnification increments for adequate integration. The effective diffusivity can vary by several orders of magnitude with varying magnficiation $\zeta$. 

Assuming that the power spectral density follows a power law of form $C(\zeta)^{1D} = C_{0} \zeta^{m}$ where $m = -2H-1$ for self affine roughness, one can show that:
$$  \frac{f'(\zeta)}{f(\zeta)} = (2+m)\zeta^{-1} $$
Thus $f'(\zeta) \le 0$ if $m\le-2$ or $H\ge 0.5$. This feature is convenient for an adaptive magnification integration scheme for roughness distributions with Hurst exponents greater than 0.5 seeing that the diffusivity decreases with increasing magnification, allowing for larger time steps with increasing magnification. Thus at a given time station $i$, we can compute a reasonable magnification step $\Delta\zeta_{i}$ to the next time station $i+1$ using the diffusivity evaluated at the time station i $F_{i} = f(\zeta_{i})$:

$$    \Delta\zeta_{i} = \alpha \frac{\Delta \sigma^2}{F^{i}} $$

We numerically compute the solutions to Equation~\ref{eq:govern} subject to the boundary condition Equation~\ref{eq:gov1} which is evolved by Equation~\ref{eq:gov2} for a suite of stress exponents $n$ and Hurst exponents $H$ to examine the trade-offs between elastic and plastic deformation as a function of scale (Figures~\ref{fig:area} and \ref{fig:force}). The equations are solved numerically using an implicit Crank-Nicholson finite difference scheme that is 2nd-order accurate in stress- and magnification-space.



