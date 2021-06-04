var documenterSearchIndex = {"docs":
[{"location":"gettingstarted/#Getting-Started","page":"Getting Started","title":"Getting Started","text":"","category":"section"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"We provide a brief overview of the package.","category":"page"},{"location":"gettingstarted/#Installation","page":"Getting Started","title":"Installation","text":"","category":"section"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"This package can be installed by cloning this repository:","category":"page"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"] add https://github.com/exanauts/ExaTron.jl","category":"page"},{"location":"gettingstarted/#Example","page":"Getting Started","title":"Example","text":"","category":"section"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"The following code snippet shows how to use this pacakge to solve a simple quadratic programming problem of the form","category":"page"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"min  05*(x-1)^2  textst  0 leq x leq 20","category":"page"},{"location":"gettingstarted/","page":"Getting Started","title":"Getting Started","text":"using ExaTron\n\n# callback function to evaluate objective\nqp_eval_f_cb(x) = 0.5*(x[1]-1)^2\n\n# callback function to evaluate the gradient \nfunction qp_eval_grad_f_cb(x, grad_f)\n    grad_f[1] = x[1] - 1\nend\n\n# callback function to evaluate the Hessian\nfunction qp_eval_h_cb(x, mode, rows, cols, obj_factor, lambda, values)\n    if mode == :Structure\n        rows[1] = 1\n        cols[1] = 1\n    else\n        values[1] = 1.0\n    end\nend\n\nx_l = zeros(1)\nx_u = zeros(1)\nx_u[1] = 2.0\nobj = 0.0\n\nprob = ExaTron.createProblem(1, x_l, x_u, 1, qp_eval_f_cb, qp_eval_grad_f_cb, qp_eval_h_cb)\n\nExaTron.solveProblem(prob)","category":"page"},{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"This page has been automatically generated by Documenter.jl to list the functions.","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [ExaTron]","category":"page"},{"location":"api/#ExaTron.dbreakpt-NTuple{5, Any}","page":"API","title":"ExaTron.dbreakpt","text":"Subroutine dbreakpt\n\nThis subroutine computes the number of break-points, and the minimal and maximal break-points of the projection of x + alpha*w on the n-dimensional interval [xl,xu].\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.dcauchy-NTuple{10, Any}","page":"API","title":"ExaTron.dcauchy","text":"Subroutine dcauchy\n\nThis subroutine computes a Cauchy step that satisfies a trust region constraint and a sufficient decrease condition.\n\nThs Cauchy step is computed for the quadratic\n\nq(s) = 0.5*s'*A*s + g'*s\n\nwhere A is a symmetric matrix in compressed row storage, and g is a vector. Given a parameter alpha, the Cauchy step is\n\ns[alpha] = P[x - alpha*g] - x,\n\nwith P the projection onto the n-dimensional interval [xl,xu]. The Cauchy step satisfies the trust region constraint and the sufficient decrease condition\n\n|| s || <= delta,    q(s) <= mu_0*(g'*s),\n\nwhere mu_0 is a constant in (0,1).\n\nMINPACK-2 Project. March 1999. Argonne National Laboratory. Chih-Jen Lin and Jorge J. More'.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.dgpstep-NTuple{7, Any}","page":"API","title":"ExaTron.dgpstep","text":"Subroutine dgpstep\n\nThis subroutine computes the gradient projection step\n\ns = P[x + alpha*w] - x,\n\nwhere P is the projection on the n-dimensional interval [xl,xu].\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.dicf-Tuple{Any, Any, TronSparseMatrixCSC, Any, Any, Any, Any, Any}","page":"API","title":"ExaTron.dicf","text":"Subroutine dicf\n\nGiven a sparse symmetric matrix A in compressed row storage, this subroutine computes an incomplete Cholesky factorization.\n\nImplementation of dicf is based on the Jones-Plassmann code. Arrays indf and list define the data structure. At the beginning of the computation of the j-th column,\n\nFor k < j, indf[k] is the index of A for the first\nnonzero l[i,k] in the k-th column with i >= j.\n\nFor k < j, list[i] is a pointer to a linked list of column\nindices k with i = L.rowval[indf[k]].\n\nFor the computation of the j-th column, the array indr records the row indices. Hence, if nlj is the number of nonzeros in the j-th column, then indr[1],...,indr[nlj] are the row indices. Also, for i > j, indf[i] marks the row indices in the j-th column so that indf[i] = 1 if l[i,j] is not zero.\n\nMINPACK-2 Project. May 1998. Argonne National Laboratory. Chih-Jen Lin and Jorge J. More'.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.dicfs-NTuple{9, Any}","page":"API","title":"ExaTron.dicfs","text":"Subroutine dicfs\n\nGiven a symmetric matrix A in compreessed column storage, this subroutine computes an incomplete Cholesky factor of A + alpha*D, where alpha is a shift and D is the diagonal matrix with entries set to the l2 norms of the columns of A.\n\nMINPACK-2 Project. October 1998. Argonne National Laboratory. Chih-Jen Lin and Jorge J. More'.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.dmid-NTuple{4, Any}","page":"API","title":"ExaTron.dmid","text":"Subroutine dmid\n\nThis subroutine computes the projection of x on the n-dimensional interval [xl,xu].\n\nMINPACK-2 Project. March 1999. Argonne National Laboratory. Chih-Jen Lin and Jorge J. More'.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.dprsrch-NTuple{9, Any}","page":"API","title":"ExaTron.dprsrch","text":"Subroutine dprsrch\n\nThis subroutine uses a projected search to compute a step that satisfies a sufficient decrease condition for the quadratic\n\nq(s) = 0.5*s'*A*s + g'*s,\n\nwhere A is a symmetric matrix in compressed column storage, and g is a vector. Given the parameter alpha, the step is\n\ns[alpha] = P[x + alpha*w] - x,\n\nwhere w is the search direction and P the projection onto the n-dimensional interval [xl,xu]. The final step s = s[alpha] satisfies the sufficient decrease condition\n\nq(s) <= mu_0*(g'*s),\n\nwhere mu_0 is a constant in (0,1).\n\nThe search direction w must be a descent direction for the quadratic q at x such that the quadratic is decreasing in the ray x + alpha*w for 0 <= alpha <= 1.\n\nMINPACK-2 Project. March 1999. Argonne National Laboratory. Chih-Jen Lin and Jorge J. More'.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.dsel2-NTuple{4, Any}","page":"API","title":"ExaTron.dsel2","text":"Subroutine dsel2\n\nGiven an array x, this subroutine permutes the elements of the array keys so that\n\nabs(x(keys(i))) <= abs(x(keys(k))),  1 <= i <= k,   abs(x(keys(k))) <= abs(x(keys(i))),  k <= i <= n.\n\nIn other words, the smallest k elements of x in absolute value are x(keys(i)), i = 1,...,k, and x(keys(k)) is the kth smallest element.\n\nMINPACK-2 Project. March 1998. Argonne National Laboratory. William D. Kastak, Chih-Jen Lin, and Jorge J. More'.\n\nRevised October 1999. Length of x was incorrectly set to n.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.dspcg-NTuple{18, Any}","page":"API","title":"ExaTron.dspcg","text":"Subroutine dspcg\n\nThis subroutine generates a sequence of approximate minimizers for the subproblem\n\nmin { q(x) : xl <= x <= xu }.\n\nThe quadratic is defined by\n\nq(x[0]+s) = 0.5*s'*A*s + g'*s,\n\nwhere x[0] is a base point provided by the user, A is a symmetric matrix in compressed column storage, and g is a vector.\n\nAt each stage we have an approximate minimizer x[k], and generate a direction p[k] by using a preconditioned conjugate gradient method on the subproblem\n\nmin { q(x[k]+p) : || L'*p || <= delta, s(fixed) = 0 },\n\nwhere fixed is the set of variables fixed at x[k], delta is the trust region bound, and L is an incomplete Cholesky factorization of the submatrix\n\nB = A(free:free),\n\nwhere free is the set of free variables at x[k]. Given p[k], the next minimizer x[k+1] is generated by a projected search.\n\nThe starting point for this subroutine is x[1] = x[0] + s, where x[0] is a base point and s is the Cauchy step.\n\nThe subroutine converges when the step s satisfies\n\n|| (g + A*s)[free] || <= rtol*|| g[free] ||\n\nIn this case the final x is an approximate minimizer in the face defined by the free variables.\n\nThe subroutine terminates when the trust region bound does not allow further progress, that is, || L'*p[k] || = delta. In this case the final x satisfies q(x) < q(x[k]).\n\nMINPACK-2 Project. March 1999. Argonne National Laboratory. Chih-Jen Lin and Jorge J. More'.\n\nMarch 2000\n\nClarified documentation of nv variable. Eliminated the nnz = max(nnz,1) statement.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.dssyax-Tuple{TronSparseMatrixCSC, Any, Any}","page":"API","title":"ExaTron.dssyax","text":"Subroutine dssyax\n\nThis subroutine computes the matrix-vector product y = A*x, where A is a symmetric matrix with the strict lower triangular part in compressed column storage.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.dtron-NTuple{23, Any}","page":"API","title":"ExaTron.dtron","text":"Subroutine dtron\n\nThis subroutine implements a trust region Newton method for the solution of large bound-constrained optimization problems\n\nmin { f(x) : xl <= x <= xu }\n\nwhere the Hessian matrix is sparse. The user must evaluate the function, gradient, and the Hessian matrix.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.dtrpcg-NTuple{14, Any}","page":"API","title":"ExaTron.dtrpcg","text":"Subroutine dtrpcg\n\nGiven a sparse symmetric matrix A in compressed column storage, this subroutine uses a preconditioned conjugate gradient method to find an approximate minimizer of the trust region subproblem\n\nmin { q(s) : || L'*s || <= delta }.\n\nwhere q is the quadratic\n\nq(s) = 0.5s'As + g's,\n\nA is a symmetric matrix in compressed column storage, L is a lower triangular matrix in compressed column storage, and g is a vector.\n\nThis subroutine generates the conjugate gradient iterates for the equivalent problem\n\nmin { Q(w) : || w || <= delta },\n\nwhere Q is the quadratic defined by\n\nQ(w) = q(s),        w = L'*s.\n\nTermination occurs if the conjugate gradient iterates leave the trust regoin, a negative curvature direction is generated, or one of the following two convergence tests is satisfied.\n\nConvergence in the original variables:\n\n|| grad q(s) || <= tol\n\nConvergence in the scaled variables:\n\n|| grad Q(w) || <= stol\n\nNote that if w = L's, then Lgrad Q(w) = grad q(s).\n\nMINPACK-2 Project. March 1999. Argonne National Laboratory. Chih-Jen Lin and Jorge J. More'.\n\nAugust 1999\n\nCorrected documentation for l, ldiag, lcolptr, and lrowind.\n\nFebruary 2001\n\nWe now set iters = 0 in the special case g = 0.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.dtrqsol-NTuple{4, Any}","page":"API","title":"ExaTron.dtrqsol","text":"Subroutine dtrqsol\n\nThis subroutine computes the largest (non-negative) solution of the quadratic trust region equation\n\n||x + sigma*p|| = delta.\n\nThe code is only guaranteed to produce a non-negative solution if ||x|| <= delta, and p != 0. If the trust region equation has no solution, sigma = 0.\n\nMINPACK-2 Project. March 1999. Argonne National Laboratory. Chih-Jen Lin and Jorge J. More'.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.ihsort-Tuple{Any, Any}","page":"API","title":"ExaTron.ihsort","text":"Subroutine ihsort\n\nGiven an integer array keys of length n, this subroutine uses a heap sort to sort the keys in increasing order.\n\nThis subroutine is a minor modification of code written by Mark Jones and Paul Plassmann.\n\nMINPACK-2 Project. March 1998. Argonne National Laboratory. Chih-Jen Lin and Jorge J. More'.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.insort-Tuple{Any, Any}","page":"API","title":"ExaTron.insort","text":"Subroutine insort\n\nGiven an integer array keys of length n, this subroutine uses an insertion sort to sort the keys in increasing order.\n\nMINPACK-2 Project. March 1998. Argonne National Laboratory. Chih-Jen Lin and Jorge J. More'.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.tron_daxpy-NTuple{6, Any}","page":"API","title":"ExaTron.tron_daxpy","text":"Subroutine daxpy\n\nThis subroutine computes constant times a vector plus a vector. It uses unrolled loops for increments equal to one. Jack Dongarra, LINPACK, 3/11/78.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.tron_ddot-NTuple{5, Any}","page":"API","title":"ExaTron.tron_ddot","text":"Subroutine ddot\n\nThis subroutine forms the dot product of two vectors. It uses unrolled loops for increments equal to one. Jack Dongarra, LINPACK, 3/11/78.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.tron_dnrm2-Tuple{Any, Any, Any}","page":"API","title":"ExaTron.tron_dnrm2","text":"DNRM2 returns the euclidean norm of a vector via the function name, so that\n\nDNRM2 := sqrt( x'*x )\n\n– This version written on 25-October-1982.    Modified on 14-October-1993 to inline the call to DLASSQ.    Sven Hammarling, Nag Ltd.\n\n\n\n\n\n","category":"method"},{"location":"api/#ExaTron.tron_kernel-Tuple{Int64, Int64, Int64, Int64, Float64, Float64, Bool, CUDA.CuDeviceVector{Float64, A} where A, CUDA.CuDeviceVector{Float64, A} where A, CUDA.CuDeviceVector{Float64, A} where A, CUDA.CuDeviceMatrix{Float64, A} where A, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}","page":"API","title":"ExaTron.tron_kernel","text":"Driver to run TRON on GPU. This should be called from a kernel.\n\n\n\n\n\n","category":"method"},{"location":"#ExaTron.jl-Documentation","page":"Home","title":"ExaTron.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"ExaTron is the Julia implementaion of novel GPU-accelerated algorithm for bound-constrained nonlinear nonconvex optimization problems of the form","category":"page"},{"location":"","page":"Home","title":"Home","text":"min_x  f(x)  textsubject to  l leq x leq u","category":"page"},{"location":"","page":"Home","title":"Home","text":"where x in mathbfR^d is the optimization variable and lu in mathbfR^d cup -inftyinfty^d are respectively lower and upper bounds (allowing negative and positive infinite values). Bound constraints hold componentwise, and the objective function f mathbfR^d rightarrow mathbfR is a generic nonlinear nonconvex function. Bound-constrained problems play an important role as a building block to solve problems with more general constraints such as h(x)=0, where h is a linear or a nonlinear function. The algorithm is a variant of TRON (Lin and Moré, 1999) with the complete Cholesky factorization for preconditioning, which has been carefully designed for solving extremely many small nonlinear problems as GPU batching.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package also provides the implementation of adaptive ADMM for solving large-scale alternating current optimal power flow by using the algorithm on multiple NVIDIA GPUs (or CPUs). ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"admm.md\"]\nDepth = 2","category":"page"},{"location":"#Citing-this-package","page":"Home","title":"Citing this package","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"@misc{ExaTron.jl.0.0.0,\n  author       = {Kim, Youngdae and Pacaud, Fran\\ccois and Kim, Kibaek},\n  title        = {{ExaTron.jl: GPU-capable TRON solver in Julia}},\n  month        = Mar,\n  year         = 2021,\n  version      = {0.0.0},\n  url          = {https://github.com/exanauts/ExaTron.jl}\n}","category":"page"},{"location":"#Acknowledgements","page":"Home","title":"Acknowledgements","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This material is based upon work supported by the U.S. Department of Energy, Office of Science, under contract number DE-AC02-06CH11357.","category":"page"},{"location":"admm/#Distributed-Optimization-of-ACOPF","page":"Distributed Optimization of ACOPF","title":"Distributed Optimization of ACOPF","text":"","category":"section"},{"location":"admm/","page":"Distributed Optimization of ACOPF","title":"Distributed Optimization of ACOPF","text":"This presents the use case of ExaTron.jl for solving large-scale alternating current optimal power flow (ACOPF) problem. In this pacakge, we also provide the implementation of adaptive ADMM for distributed ACOPF introduced by Mhanna et al. (2019). We have implemented the ADMM algorithm fully on GPUs without data transfer to the CPU, where ExaTron.jl is used to solve many small nonlinear nonconvex problems, each of which represents a branch subproblem of the ADMM.","category":"page"},{"location":"admm/#Numerical-Experiment","page":"Distributed Optimization of ACOPF","title":"Numerical Experiment","text":"","category":"section"},{"location":"admm/","page":"Distributed Optimization of ACOPF","title":"Distributed Optimization of ACOPF","text":"All experiments were performed on a compute node of the Summit supercomputer at Oak Ridge National Laboratory using Julia@1.6.0 and CUDA.jl@2.6.1. Note, however, that our implementation is not limited to a single node. Each compute node of the Summit supercomputer has 2 sockets of POWER9 processors having 22 physical cores each, 512 GB of DRAM, and 6 NVIDIA Tesla V100 GPUs evenly distributed to each socket.","category":"page"},{"location":"admm/#ACOPF-Problem-Instances","page":"Distributed Optimization of ACOPF","title":"ACOPF Problem Instances","text":"","category":"section"},{"location":"admm/","page":"Distributed Optimization of ACOPF","title":"Distributed Optimization of ACOPF","text":"The following table presents the data statistics of our test examples from MATPOWER and PGLIB benchmark instances. We note that up to 34,000 nonlinear nonconvex problems are solved by our solver at each ADMM iteration.","category":"page"},{"location":"admm/","page":"Distributed Optimization of ACOPF","title":"Distributed Optimization of ACOPF","text":"Data # Generators # Branches # Buses\n2868rte 600 3,808 2,868\n6515rte 1,389 9,037 6,515\n9241pegase 1,445 16,049 9,241\n13659pegase 4,092 20,467 13,659\n19402goc 971 34,704 19,402","category":"page"},{"location":"admm/#Weak-Scaling:-Performance-on-a-single-GPU","page":"Distributed Optimization of ACOPF","title":"Weak Scaling: Performance on a single GPU","text":"","category":"section"},{"location":"admm/","page":"Distributed Optimization of ACOPF","title":"Distributed Optimization of ACOPF","text":"The following figure depicts the average solution time of ExaTron.jl for different sizes of batches of branch subproblems. The time on the y-axis is the average computation time in milliseconds taken by ExaTron.jl to solve each batch within an ADMM iteration.","category":"page"},{"location":"admm/","page":"Distributed Optimization of ACOPF","title":"Distributed Optimization of ACOPF","text":"(Image: )","category":"page"},{"location":"admm/#Strong-Scaling:-Performance-on-multiple-GPUs","page":"Distributed Optimization of ACOPF","title":"Strong Scaling: Performance on multiple GPUs","text":"","category":"section"},{"location":"admm/","page":"Distributed Optimization of ACOPF","title":"Distributed Optimization of ACOPF","text":"The following figure shows the speedup of ExaTron.jl when we parallelize the computation across different GPUs (up to the 6 GPUs available on a node in the Summit supercomputer). Branch problems are evenly dispatched among 6 MPI processes in the order of branch indices, and the speedup is computed based on the timing of the root process.","category":"page"},{"location":"admm/","page":"Distributed Optimization of ACOPF","title":"Distributed Optimization of ACOPF","text":"(Image: )","category":"page"},{"location":"admm/#Performance-comparison:-6-GPUs-vs.-40-CPUs","page":"Distributed Optimization of ACOPF","title":"Performance comparison: 6 GPUs vs. 40 CPUs","text":"","category":"section"},{"location":"admm/","page":"Distributed Optimization of ACOPF","title":"Distributed Optimization of ACOPF","text":"This experiment was run on a single Summit node with 6 GPUs and 40 CPUs. For the CPU run, we use the MPI library to implement the parallel communication between the CPU processes. In the following figure, the computation time of the CPU implementation shows a linear increase of with respect to the batch size. However, the average computation time increases faster than that of the GPU implementation: the computation time of ExaTron.jl on 6 GPUs is up to 35 times faster than the CPU implementation using 40 cores. Most of the speedup relates to the GPU's massive parallel computation capability.","category":"page"},{"location":"admm/","page":"Distributed Optimization of ACOPF","title":"Distributed Optimization of ACOPF","text":"(Image: )","category":"page"}]
}
