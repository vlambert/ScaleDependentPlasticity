# Elastoplastic contact with scale-dependent plastic strength
This directory includes a numerical program to solve for the distribution of elastic, plastic and non-contact area along a randomly rough surface at varying scales including a scale-dependent plastic yield stress following the theoretical of Persson ([2001](https://doi.org/10.1016/S0039-6028(98)00051-X),[2006](https://doi.org/10.1016/j.surfrep.2006.04.001))

$ \frac{\partial P(\sigma,\zeta)}{\partial \zeta} = f(\zeta) \frac{\partial^2 P(\sigma,\zeta)}{\partial \sigma^2}$ 



We numerically compute the solutions to Equation~\ref{eq:govern} subject to the boundary condition Equation~\ref{eq:gov1} which is evolved by Equation~\ref{eq:gov2} for a suite of stress exponents $n$ and Hurst exponents $H$ to examine the trade-offs between elastic and plastic deformation as a function of scale (Figures~\ref{fig:area} and \ref{fig:force}). The equations are solved numerically using an implicit Crank-Nicholson finite difference scheme that is 2nd-order accurate in stress- and magnification-space.
