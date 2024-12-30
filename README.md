# Elastoplastic contact with scale-dependent plastic strength
This directory includes a numerical program to solve for the distribution of elastic, plastic and non-contact area along a randomly rough surface at varying scales including a scale-dependent plastic yield stress following the theoretical developments of Persson ([2001](https://doi.org/10.1016/S0039-6028(98)00051-X);[2006](https://doi.org/10.1016/j.surfrep.2006.04.001)).

## Problem Statement
The problem statement considers a rough surface with linear scale $$L$$ that can be examined at different magnifications $$\zeta = L/\lambda$$, where $$\lambda$$ is the shortest wavelength of roughness which is resolved at magnification $$\zeta$$. 

$$ \frac{\partial P(\sigma,\zeta)}{\partial \zeta} = f(\zeta) \frac{\partial^2 P(\sigma,\zeta)}{\partial \sigma^2} $$ 

## Numerical Solution
We solve the governing equation subject to the evolving yield stress boundary condition using implicit Crank-Nicholson finite difference scheme. The numerical solution is implement using a modified Thomas Algorithm which provides a parallilzed framework for efficiently solving the system of equations characterizerd by a tri-diagonalize matrix.

The interior solution of the diffusion problem is solved with a standard finite difference scheme, with the solutions to the governing equation given at stress increments $j$ and magnification step $i+1$ in the form:

$$    a_{j}^{i+1}P_{j-1}^{i+1}+b_{j}^{i+1}P_{j}^{i+1} + c_{j}^{i+1}P_{j+1}^{i+1} = d_{j}^{i} $$

Following an implicit Crank-Nicolson scheme, the coefficients take the form:

$$    a_{j}^{i+1} = -\frac{\Delta\zeta}{2}\frac{F^{i+1}}{\Delta\sigma^2} $$

$$    b_{j}^{i+1} = 1+ 2\frac{\Delta\zeta}{2}\frac{F^{i+1}}{\Delta\sigma^2} $$

$$    c_{j}^{i+1} =-\frac{\Delta\zeta}{2}\frac{F^{i+1}}{\Delta\sigma^2} $$

$$    d_{j}=\frac{\Delta\zeta}{2}\frac{F^{i}}{\Delta\sigma^2}\bigg( P_{j-1}^{i} + P_{j+1}^{i}\bigg) + \bigg(1-2\frac{\Delta\zeta}{2}\frac{F^{i}}{\Delta\sigma^2}\bigg) P_{j}^{i} $$

which can all be evaluated based on information from the current magnification step $\zeta_i$ and knowledge of the diffusivity $F^{i+1} = f(\zeta_{i+1})$, which is precribed in our calculations. The system of equations form matrix-vector product with a tridiagonal matrix and can be solved in parallel using a modified Thomas Algorithm. This implicit Crank-Nicolson scheme is numerically stable and 2nd-order accurate in stress and magnficiation-space (space and time) within the domain $[0,Y(\zeta))$.

The solution for the B.C. at $\sigma=0$ takes the form:

$$    b_{1}P_{0}^{i+1} + c_{1}P_{1}^{i+1} = d_{1}^{i}, $$

where $b_{1} = 1$ and $c_{1}=d_{1} = 0$ to ensure $P_{0} = 0$. Similarly for a constant yield stress ($Y'(\zeta) =0$) $P_{ny} = 0$ where $ny$ indicates the node at the maximum end of the discretized stress axis (i.e. $P(Y(\zeta),\zeta))$. The $\sigma = Y(\zeta)$ boundary condition takes the form:

$$    a_{ny}P_{ny-1}^{i+1} + b_{ny}P_{ny}^{i+1} = d_{ny}^{i} $$

with $b_{ny} = 1$ and $a_{ny}=d_{ny} = 0$ for the specific case of a constant yield stress.

### Discretization scheme
We construct an initial 1-D grid in stress space spanning the domain $$[0,Y(1)]$$ discretized by N grid points with uniform spacing $\Delta\sigma = Y(1)/(N-1)$. Magnification increments of the form:

$$    \Delta\zeta = \alpha \frac{\Delta \sigma^2}{f(\zeta)}  $$

are chosen to satisfy the Courant conditon. Higher diffusivities $f(\zeta)$ require smaller magnification increments for adequate integration. The effective diffusivity can vary by several orders of magnitude with varying magnficiation $\zeta$. 

For a self-affine surface where the roughness power spectral density is a power-law function of scale $C(\zeta) = C_{0} \zeta^{m}$ where $C_0$ is a pre-factor and $m = -2H-1$, the change in the effective diffusivity with increasing magnification is given by:

$$  \frac{f'(\zeta)}{f(\zeta)} = (2+m)\zeta^{-1} $$

Thus $f'(\zeta) \le 0$ if $m\le-2$ or for surfaces with Hurst exponents $H\ge 0.5$, as is common for many natural and engineered surfaces. This outcome is convenient for an adaptive magnification integration scheme for roughness distributions with Hurst exponents greater than 0.5 seeing that the effective diffusivity decreases with increasing magnification, allowing for larger time steps with increasing magnification. At a given magnification step $\zeta_i$, we compute a reasonable magnification step $\Delta\zeta_{i}$ to the next magnification step $\zeta_{i+1}$ using the diffusivity evaluated at step $\zeta_i$, $F_{i} = f(\zeta_{i})$:

$$    \Delta\zeta_{i} = \alpha \frac{\Delta \sigma^2}{F^{i}} $$

## Treatment of evolving yield stress boundary
For a scale-dependent yield stress where the yield stress $Y(\zeta)$ is a function of magnification, we need to both adapt the model domain $\sigma \in [0,Y(\zeta)]$ as a function of magnification and solve for the appropriate boundary conditions $P(Y(\zeta),\zeta)$. Conservation of force across scales requires that $Y'(\zeta) \ge 0$ so we only consider problems where the stress domain increases with increasing magnficiation $\zeta$. 

Let $F^{i}$ and $Sp^{i}$ denote the diffusivity $f(\zeta_i)$ and derivative of the yield stress $Y'(\zeta)$ at the magnification step $\zeta_{i}$. Given the solution to $P_{ny}$ at magnification step $i$, i.e. $P_{ny}^{i}$, one can solve for the corresponding plastic area of contact:

$$ A_{pl}^{i} = \frac{F^{i}}{Sp^{i}} P_{ny}^{i}. $$

Similarly, the boundary condition at the new magnification step $P_{ny}^{i+1}$ should be consistent with the new plastic area of contact $A_{pl}^{i+1}$:

$$    P_{ny}^{i+1} = \frac{Sp^{i+1}}{F^{i+1}} A_{pl}^{i+1}, $$

which we approximate as:

$$    P_{ny}^{i+1}= \frac{Sp^{i+1}}{F^{i+1}} \left[ A_{pl}^{i} + \Delta \zeta \frac{\partial A_{pl}}{\partial \zeta}\right]. $$

We compute the increment $\Delta A_{pl} = \Delta\zeta \partial A_{pl}/\partial \zeta$ using the trapezoidal rule:

$$    \Delta A_{pl} = \frac{\Delta\zeta}{2}\bigg(  -Sp^{i+1}P_{ny}^{i+1} - \frac{F^{i+1}}{\Delta \sigma}\bigg[ P_{ny}^{i+1}-P_{ny-1}^{i+1}\bigg] -Sp^{i}P_{ny}^{i} - \frac{F^{i}}{\Delta \sigma}\bigg[ P_{ny}^{i}-P_{ny-1}^{i}\bigg]\bigg). $$

The above considers the uniform grid spacing $\Delta\sigma = \sigma_{j}-\sigma_{j-1}$. As the yield boundary evolves, the spacing between the yield boundary $\sigma_{ny}$ and its neighboring grid point $\sigma_{ny-1}$ may differ from the uniform grid by a small amount $\delta$ as $\Delta\sigma_{ny} = \Delta\sigma + \delta$ depending on the magnification step.

We adapt the 2nd-order finite difference method for the nodes near the yield boundary with uneven spacing by interpolating the distribution $P(\sigma,\zeta)$ with a 2nd-order Lagrange polynomial between points $P_{ny-2}$ and $P_{ny}$. The equation for our discretized diffusion equation at the second-to-last node can then be determined as:

```math
\frac{P_{j}^{i+1}-P_{j}^{i}}{\Delta\zeta} = \bigg[F^{i}\bigg(\frac{P^{i}_{j+1}}{(\Delta\sigma+\delta^{i})(2\Delta\sigma + \delta^{i})} + \frac{P^{i}_{j-1}}{\Delta\sigma(2\Delta\sigma+\delta^{i})}-\frac{P^{i}_{j}}{\Delta\sigma(\Delta\sigma+\delta^{i})}\bigg) \\
+ F^{i+1}\bigg(\frac{P^{i+1}_{j+1}}{(\Delta\sigma+\delta^{i+1})(2\Delta\sigma + \delta^{i+1})} + \frac{P^{i+1}_{j-1}}{\Delta\sigma(2\Delta\sigma+\delta^{i+1})}-\frac{P^{i+1}_{j}}{\Delta\sigma(\Delta\sigma+\delta^{i+1})}\bigg)\bigg]
```

Following an implicit Crank-Nicolson scheme, the new coefficients for the second-to-last node take the form:

$$    a_{ny-1}^{i+1} = -\Delta\zeta\frac{F^{i+1}}{\Delta\sigma (2\Delta\sigma+\delta^{i+1})} $$
$$    b_{ny-1}^{i+1} = 1+ \Delta\zeta F^{i+1}\frac{1}{\Delta\sigma(\Delta\sigma+\delta^{i+1})} $$
$$    c_{ny-1}^{i+1} =-\Delta\zeta\frac{F^{i+1}}{(\Delta\sigma+\delta^{i+1})(2\Delta\sigma+\delta^{i+1})} $$
$$    d_{ny-1}^{i}=\Delta\zeta F^{i} \bigg( \frac{P_{ny-2}^{i}}{\Delta\sigma(2\Delta\sigma + \delta^{i})} + \frac{P_{ny}^{i}}{(\Delta\sigma+\delta^{i})(2\Delta\sigma+\delta^{i})}\bigg) + \bigg(1-\Delta\zeta F^{i} \frac{1}{\Delta\sigma(\Delta\sigma + \delta^{i})}\bigg) P_{ny-1}^{i} $$

Now we consider the condition for the node at the yield boundary:

```math
    \Delta A_{pl} = \frac{\Delta\zeta}{2}\bigg(  -Sp^{i+1}P_{ny}^{i+1} - F^{i+1}\bigg[ \frac{3\Delta\sigma + 2 \delta^{i+1}}{(\Delta\sigma + \delta^{i+1})(2\Delta\sigma+\delta^{i+1})}P_{ny}^{i+1}-\frac{2\Delta\sigma + \delta^{i+1}}{\Delta\sigma(\Delta\sigma+\delta^{i+1})}P_{ny-1}^{i+1} + \frac{\Delta\sigma + \delta^{i+1}}{\Delta\sigma(2\Delta\sigma + \delta^{i+1})}P_{ny-2}^{i+1}\bigg] 
```
```math
   - Sp^{i}P_{ny}^{i} - F^{i}\bigg[ \frac{3\Delta\sigma + 2 \delta^{i}}{(\Delta\sigma + \delta^{i})(2\Delta\sigma+\delta^{i})}P_{ny}^{i}-\frac{2\Delta\sigma + \delta^{i}}{\Delta\sigma(\Delta\sigma+\delta^{i})}P_{ny-1}^{i} + \frac{\Delta\sigma + \delta^{i}}{\Delta\sigma(2\Delta\sigma + \delta^{i})}P_{ny-2}^{i}\bigg]\bigg) 
```

Note that this 2nd-order formulation of the boundary equation compromises the tridiagonal nature of the finite difference operator. We can define a new tridiagonal matrix to solve by scaling and adding the second-to-last equation to the yield boundary equation:

$$ a_{ny-1}P_{ny-2}+ b_{ny-1}P_{ny-1}+ c_{ny-1}P_{ny} = d_{ny-1} $$
$$    z_{ny}P_{ny-2}+ a_{ny}P_{ny-1}+ b_{ny}P_{ny} = d_{ny} $$

$$ \rightarrow (a_{ny}- \frac{z_{ny}}{a_{ny-1}}b_{ny-1})P_{ny-1} + (b_{ny}-\frac{z_{ny}}{a_{ny-1}}c_{ny-1}) P_{ny} = d_{ny} - \frac{z_{ny}}{a_{ny-1}}d_{ny-1} $$

The new coefficients for our final equation become:

$$ \tilde{a}_{ny} =  a_{ny} - \frac{z_{ny}}{a_{ny-1}} b_{ny-1} $$

$$ \tilde{b}_{ny} =  b_{ny}-\frac{z_{ny}}{a_{ny-1}} c_{ny-1}$$

$$    \tilde{d}_{ny} = d_{ny} - \frac{z_{ny}}{a_{ny-1}} d_{ny-1}$$

where:

$$     -\frac{z_{ny}}{a_{ny-1}} = \frac{Sp^{i+1}}{2F^{i+1}} (\Delta\sigma  + \delta^{i+1})$$

Giving:

$$ \tilde{a}_{ny}^{i+1} = - \frac{\Delta\zeta}{2}Sp^{i+1}\frac{1}{(\Delta\sigma+\delta^{i+1})} + \frac{Sp^{i+1}}{2F^{i+1}}(\Delta\sigma  + \delta^{i+1}) $$

$$ \tilde{b}_{ny}^{i+1} =  1+ \frac{\Delta\zeta}{2}\frac{Sp^{i+1}}{F^{i+1}} \bigg[ Sp^{i+1}+F^{i+1}\frac{1}{(\Delta\sigma + \delta^{i+1})}\bigg]$$
```math
   \tilde{d}_{ny}^{i+1} = \frac{Sp^{i+1}}{F^{i+1}} \bigg[ A_{pl}^{i} +   \frac{\Delta\zeta}{2}  \bigg(Sp^{i}P_{ny}^{i} - F^{i}\bigg[ \frac{1}{\Delta\sigma + \delta^{i}}P_{ny}^{i}-\frac{1}{\Delta\sigma} P_{ny-1}^{i} \bigg]\bigg) + \frac{1}{2}(\Delta\sigma  + \delta^{i+1}) P_{ny-1}^{i}\bigg ]
```


