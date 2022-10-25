- AUTHOR: Kevin Nisterenko
- FILE: README.md
- ASSIGNMENT: Programming Assignment 11 - The Traveling Salesman Problem
- COURSE: CSc 210; Fall 2021

Timing Results for big11.mtx on 5 different runs:

heuristic: cost = 3.4, 0 milliseconds
mine: cost = 3.4, 1 milliseconds
backtrack: cost = 1.4, 8577 milliseconds

heuristic: cost = 3.4, 0 milliseconds
mine: cost = 3.4, 1 milliseconds
backtrack: cost = 1.4, 8216 milliseconds

heuristic: cost = 3.4, 1 milliseconds
mine: cost = 3.4, 0 milliseconds
backtrack: cost = 1.4, 8217 milliseconds

heuristic: cost = 3.4, 0 milliseconds
mine: cost = 3.4, 0 milliseconds
backtrack: cost = 1.4, 8092 milliseconds

heuristic: cost = 3.4, 0 milliseconds
mine: cost = 3.4, 1 milliseconds
backtrack: cost = 1.4, 8027 milliseconds

Looking at the structure of each algorithm:
-Heuristic: it is clear to see that if we iterate over every vertex, and then iterate
further to find the minimum edge and get it, we are doing a quadratic operation to then
find the optimal path in the heuristic form. However, we are not guaranteed the best
possible path in the graph and we take this quadratic time to build a single path.

-Backtracking: this brute force approach is helpful in the sense that we do get the
global optimum path. However, we use the recursive backtracking to get every permutation
with a constraint that we already have a known starting vertex. Thus, we have (n-1)!
possible paths that we must test and compare before we actually get to exit the function
with the best possible path. This factorial approach is way too costly, and it is not
efficient even though it is guaranteed the least costly path/solution.

-My Approach: looking at the recursive backtracking, we know that it is one approach we
can clearly make better. Moreover, we can even use a heuristic-like approach to add
constraints to the traversal of vertices, we traverse to the minimum connected vertex just
the same as we do in the heuristic approach. This constraint reduces the scope of paths we
end up producing, since we are not looking at (n-1)! paths, but rather, looking at paths
where the immediate/local cost between two vertices is the least possible. This way, we
vastly improve its performance and can enjoy a much faster method of search using recursion.
