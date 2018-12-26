# algorithm

1. remove duplicate element, and sort by x axis.
2. prepare a large diamond region for initial triangle.
3. insert every point and update triangulation according to delaunay's rule.
4. finally remove the 4 assistant points.

# dependency

Eigen 3 (?)

# problem

* speed may be slow, some optimization is needed.
* triangle orientation could be correct when it's created, rather than repair them after creation.

# reference
[wikipedia](https://en.m.wikipedia.org/wiki/Bowyer-Watson_algorithm)

[link to c++ implemention in wiki](https://github.com/Bl4ckb0ne/delaunay-triangulation)