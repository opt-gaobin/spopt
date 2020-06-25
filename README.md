# spopt
A Matlab solver for Riemannian optimization on the symplectic Stiefel manifold

## Problems
This solver is to solve the following optimization problem,
+ min f(X), s.t.   X' J2n X = J2p
  
where X is a 2n-by-2p matrix, J2n = [0 In; -In 0], and In is the n-by-n identity matrix.
### Examples
1. the nearest symplectic matrix problem

  min ||X-A||^2_F, s.t.  X' J2n X = J2p.
  
2. minimization of the Brockett cost function

  min trace(X'AXN-2BX'), s.t.  X' J2n X = J2p.
  
3. symplectic eigenvalue problem

  min trace(X'AX), s.t.  X' J2n X = J2p.
  
4. symplectic model order reduction

  min ||M-XX^\dag M||, s.t.  X' J2n X = J2p, where X^\dag = J2p' X' J2n

# References
Bin Gao, Nguyen Thanh Son, P.-A. Absil, Tatjana Stykel
1. [Riemannian optimization on the symplectic Stiefel manifold]()
2. [Euclidean--metric-based Riemannian gradient method on the symplectic Stiefel manifold]()

# Authors
+ [Bin Gao](https://www.gaobin.cc/) (UCLouvain, Belgium)

# Copyright
Copyright (C) 2020, Bin Gao

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/)
