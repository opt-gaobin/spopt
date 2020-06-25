# spopt
A Matlab solver for Riemannian optimization on the symplectic Stiefel manifold

## Problems
This solver is to solve the following optimization problem,
> min f(X),
  s.t.   X' J2n X = J2p
  
where X is a 2n-by-2p matrix, J2n = [0 In; -In 0], and In is the n-by-n identity matrix.
### Examples
+ min ||X-A||^2_F, s.t.  X' J2n X = J2p.

# References
+ Bin Gao, Nguyen Thanh Son, P.-A. Absil, Tatjana Stykel, [Riemannian optimization on the symplectic Stiefel manifold]()

# Authors
+ [Bin Gao](https://www.gaobin.cc/) (UCLouvain, Belgium)

# Copyright
Copyright (C) 2020, Bin Gao

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/)
