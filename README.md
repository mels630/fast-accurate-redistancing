# fast-accurate-redistancing
Implementation of a fast, high-ordered redistancing method for uniform finite-difference grids, following this paper:

M. Elsey and S. Esedoglu. Fast and accurate redistancing via directional optimization. SIAM Journal on Scientific Computing, 36:1 (2014), pp. A219-A231.

To-do:
* Implement a RTS to ensure that results aren't affected by any changes
* Comment
* Convert 2D containers to use size_t types rather than signed integers for indexing
* Utilize STL algorithms, std::vector::assign, etc. to simplify code
* Consistency: place const to right of types, use Hungarian notation for variables

