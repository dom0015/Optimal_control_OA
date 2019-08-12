# Optimal_control_OA
Matlab files for preconditioner testing for optimal control problems for PDEs

# Structure of repository
1. grid - procedures for uniform triangulation on given rectange + setting boundary conditions
2. assemblers - contains procedures for assembly of matrices and rhs (general grids), Assembly_all encapsulates all procedures needed for tested problems (additionaly contains tag for non-uniform grid smoother around the boundary)
3. plotting - ploting procedures for plot grid (with numbering), solution, gradient of solution, difference on Dirichlet boundary.
4. itersolvers - so far empty
5. preconditioners - so far empty

# Entry points 
- scripts starting with "test_" containing selected examples 

# Two examples
1. model problem according to sent pdf
2. problem with exact solution u=x*x+y*y