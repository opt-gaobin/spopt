# spopt
A Matlab solver for Riemannian **opt**imization on the **s**ym**p**lectic Stiefel manifold

## Problems
This solver is to solve the following optimization problem,

+ min f(X), s.t.   X' J2n X = J2p
  
where X is a 2n-by-2p matrix, J2n = [0 In; -In 0], and In is the n-by-n identity matrix.
## Applications
1. the nearest symplectic matrix problem:

> min ||X-A||^2_F, s.t.  X' J2n X = J2p.

2. the extrinsic mean problem:
  
> min 1/N \sum_{i=1}^{i=N}||X-A_i||^2_F, s.t.  X' J2n X = J2p.
  
3. minimization of the Brockett cost function:

> min trace(X'AXN-2BX'), s.t.  X' J2n X = J2p.
  
4. symplectic eigenvalue computation of spd or spsd:

> min trace(X'AX), s.t.  X' J2n X = J2p.
  
5. symplectic model order reduction:

> min ||M-XX^\dag M||, s.t.  X' J2n X = J2p, where X^\dag = J2p' X' J2n

## References
1. [Riemannian optimization on the symplectic Stiefel manifold](https://doi.org/10.1137/20M1348522), SIAM Journal on Optimization, 31-2 (2021), 1546-1575.
2. [Geometry of the symplectic Stiefel manifold endowed with the Euclidean metric](https://doi.org/10.1007/978-3-030-80209-7_85), Geometric Science of Information: 5th International Conference, GSI 2021, Lecture Notes in Computer Science, 12829 (2021), 789–796.
3. [Computing symplectic eigenpairs of symmetric positive-definite matrices via trace minimization and Riemannian optimization](https://doi.org/10.1137/21M1390621),  SIAM Journal on Matrix Analysis and Applications, 42-4 (2021), 1732–1757.
4. [Optimization on the symplectic Stiefel manifold: SR decomposition-based retraction and applications](https://arxiv.org/abs/2211.09481), arXiv:2211.09481, (2022).

## Authors
+ [P.-A. Absil](https://sites.uclouvain.be/absil/) (UCLouvain, Belgium)
+ [Bin Gao](https://www.gaobin.cc/) (UCLouvain, Belgium)
+ [Nguyen Thanh Son](https://sites.google.com/view/ntson) (Thai Nguyen University of Sciences, Vietnam)
+ [Tatjana Stykel](https://www.uni-augsburg.de/en/fakultaet/mntf/math/prof/numa/team/tatjana-stykel/) (University of Augsburg, Germany)


## Copyright
Copyright (C) 2020, P.-A. Absil, Bin Gao, Nguyen Thanh Son, Tatjana Stykel

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/)
