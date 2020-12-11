mutable struct ExaTronProblem{VI, VD}
    n::Cint                 # number of variables
    nnz::Int                # number of Hessian entries
    nnz_a::Int              # number of Hessian entries in the strict lower
    A::TronSparseMatrixCSC{VI, VD}
    B::TronSparseMatrixCSC{VI, VD}
    L::TronSparseMatrixCSC{VI, VD}
    indfree::VI   # a working array of dimension n
    iwa::VI       # a working array of dimension 3*n
    g::VD         # gradient
    xc::VD        # a working array of dimension n
    s::VD         # a working array of dimension n
    wa::VD        # a working array of dimension 7*n

    x::VD
    x_l::VD
    x_u::VD
    rows::VI
    cols::VI
    values::VD

    eval_f::Function
    eval_grad_f::Function
    eval_h::Function

    options::Dict{String, Any}

    # Statistics
    delta::Cdouble
    f::Cdouble
    gnorm_two::Float64
    gnorm_inf::Float64
    nfev::Int
    ngev::Int
    nhev::Int
    minor_iter::Int
    status::Symbol

    ExaTronProblem{VI, VD}() where {VI, VD} = new{VI, VD}()
end

function gpnorm(n, x, x_l, x_u, g)
    two_norm = 0.0
    inf_norm = 0.0

    for i=1:n
        if x_l[i] != x_u[i]
            if x[i] == x_l[i]
                val = (min(g[i], 0.0))^2
            elseif x[i] == x_u[i]
                val = (max(g[i], 0.0))^2
            else
                val = g[i]^2
            end

            two_norm += val
            val = sqrt(val)
            if inf_norm < val
                inf_norm = val
            end
        end
    end

    two_norm = sqrt(two_norm)

    return two_norm, inf_norm
end

function set_default_options!(prob::ExaTronProblem)
    prob.options = Dict{String,Any}()

    prob.options["max_feval"] = 500
    prob.options["max_minor"] = 200
    prob.options["p"] = 5
    prob.options["verbose"] = 0
    prob.options["tol"] = 1e-6
    prob.options["fatol"] = 0
    prob.options["frtol"] = 1e-12
    prob.options["fmin"] = -1e32
    prob.options["cgtol"] = 0.1
    prob.options["tron_code"] = :Julia

    return
end

function getOption(prob::ExaTronProblem, keyword::String)
    if !haskey(prob.options, keyword)
        error("ExaTron does not have option with name $(keyword).")
    end
    return prob.options[keyword]
end

function setOption(prob::ExaTronProblem, keyword::String, value::Symbol)
    if !haskey(prob.options, keyword)
        error("ExaTron does not have option with name $(keyword).")
    end
    prob.options[keyword] = value
end

function setOption(prob::ExaTronProblem, keyword::String, value::Integer)
    if !haskey(prob.options, keyword)
        error("ExaTron does not have option with name $(keyword).")
    end
    prob.options[keyword] = value
    return
end

function setOption(prob::ExaTronProblem, keyword::String, value::Float64)
    if !haskey(prob.options, keyword)
        error("ExaTron does not have option with name $(keyword).")
    end
    prob.options[keyword] = value
    return
end

function instantiate_memory!(tron::ExaTronProblem{VI,VD}, n, nele_hess) where {VI,VD}
    tron.n = convert(Cint, n)
    tron.nnz = nele_hess
    tron.indfree = VI(undef, n)
    tron.iwa = VI(undef, 3*n)

    tron.g = tron_zeros(VD, n)
    tron.xc = tron_zeros(VD, n)
    tron.s = tron_zeros(VD, n)
    tron.wa = tron_zeros(VD, 7*n)

    tron.x = tron_zeros(VD, n)
    tron.x_l = tron_zeros(VD, n)
    tron.x_u = tron_zeros(VD, n)
    tron.rows = VI(undef, nele_hess)
    tron.cols = VI(undef, nele_hess)
    tron.values = tron_zeros(VD, nele_hess)
end

function createProblem(n::Int, x_l::AbstractVector{Float64}, x_u::AbstractVector{Float64},
                       nele_hess::Int, eval_f_cb, eval_grad_f_cb, eval_h_cb)
    @assert n == length(x_l) == length(x_u)
    @assert typeof(x_l) == typeof(x_u)

    VI = Vector{Cint}
    VD = typeof(x_l)

    tron = ExaTronProblem{VI, VD}()
    set_default_options!(tron)
    instantiate_memory!(tron, n, nele_hess)
    copyto!(tron.x_l, 1, x_l, 1, n)
    copyto!(tron.x_u, 1, x_u, 1, n)

    tron.eval_f = eval_f_cb
    tron.eval_grad_f = eval_grad_f_cb
    tron.eval_h = eval_h_cb

    tron.eval_h(tron.x, :Structure, tron.rows, tron.cols, 1.0, Float64[], Float64[])
    # Instantiate sparse matrix
    p = tron.options["p"]
    tron.A = TronSparseMatrixCSC(tron.rows, tron.cols, tron.values, n)
    tron.B = TronSparseMatrixCSC{VI, VD}(n, nele_hess)
    tron.L = TronSparseMatrixCSC{VI, VD}(n, nele_hess + n*p)
    tron.nnz_a = tron.A.nnz
    tron.status = :NotSolved

    return tron
end

function solveProblem(tron::ExaTronProblem)
    task = Vector{UInt8}(undef, 60)
    for (i,s) in enumerate("START")
        task[i] = UInt8(s)
    end

    VI = Vector{Cint}
    VD = typeof(tron.x)

    isave = VI(undef, 3)
    dsave = tron_zeros(VD, 3)
    fatol = tron.options["fatol"]
    frtol = tron.options["frtol"]
    fmin = tron.options["fmin"]
    cgtol = tron.options["cgtol"]
    gtol = tron.options["tol"]
    max_feval = tron.options["max_feval"]
    tcode = tron.options["tron_code"]
    max_minor = tron.options["max_minor"]

    # Project x into its bound. Tron requires x to be feasible.
    tron.x .= max.(tron.x_l, min.(tron.x_u, tron.x))

    itermax = tron.n
    tron.minor_iter = 0
    tron.nfev = tron.ngev = tron.nhev = 0
    tron.status = :NotSolved
    search = true
    while (search)

        # Evaluate the function.

        if Char(task[1]) == 'F' || unsafe_string(pointer(task), 5) == "START"
            tron.f = tron.eval_f(tron.x)
            tron.nfev += 1
            if tron.nfev >= max_feval
                tron.status = :Maximum_Iterations_Exceeded
                search = false
            end
        end

        # Evaluate the gradient and Hessian.

        if (Char(task[1]) == 'G' && Char(task[2]) == 'H') || unsafe_string(pointer(task), 5) == "START"
            tron.eval_grad_f(tron.x, tron.g)
            tron.eval_h(tron.x, :Values, Int[], Int[], 1.0, Float64[], tron.values)
            tron.ngev += 1
            tron.nhev += 1
            tron.minor_iter += 1

            # Copy values in the CSC matrix.
            fill!(tron.A, 0.0)
            @inbounds for i in 1:tron.nnz
                m = tron.A.map[i]
                if m < 0
                    tron.A.diag_vals[-m] += tron.values[i]
                else
                    tron.A.tril_vals[m] += tron.values[i]
                end
            end
        end

        # Initialize the trust region bound.

        if unsafe_string(pointer(task), 5) == "START"
            gnorm0 = norm(tron.g, 2)
            tron.delta = gnorm0
        end

        # Call Tron.

        if search
            if tcode == :Fortran
                delta = Ref{Cdouble}(tron.delta)
                # TODO
                ccall((:dtron_, EXATRON_LIBRARY),
                    Cvoid,
                    (Ref{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                    Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                    Ptr{Cint}, Ptr{Cint}, Ref{Cdouble}, Ref{Cdouble},
                    Ref{Cdouble}, Ref{Cdouble}, Ref{Cint}, Ref{Cdouble},
                    Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint},
                    Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint},
                    Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint},
                    Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}),
                    tron.n, tron.x, tron.x_l, tron.x_u,
                    tron.f, tron.g, tron.A.tril_vals, tron.A.diag_vals,
                    tron.A.colptr, tron.A.rowval, frtol, fatol,
                    fmin, cgtol, itermax, delta,
                    task, tron.B.tril_vals, tron.B.diag_vals, tron.B.colptr,
                    tron.B.rowval, tron.L.tril_vals, tron.L.diag_vals, tron.L.colptr,
                    tron.L.rowval, tron.xc, tron.s, tron.indfree,
                    isave, dsave, tron.wa, tron.iwa)
                tron.delta = delta[]
            else
                tron.delta = ExaTron.dtron(tron.n, tron.x, tron.x_l, tron.x_u,
                                tron.f, tron.g, tron.A,
                                frtol, fatol,
                                fmin, cgtol, itermax, tron.delta,
                                task, tron.B, tron.L,
                                tron.xc, tron.s, tron.indfree,
                                isave, dsave, tron.wa, tron.iwa)
            end
        end

        if unsafe_string(pointer(task), 4) == "NEWX"
            tron.gnorm_two, tron.gnorm_inf = gpnorm(tron.n, tron.x, tron.x_l, tron.x_u, tron.g)
            if tron.gnorm_inf <= gtol
                for (i,s) in enumerate("CONVERGENCE: GTOL TEST SATISFIED")
                    task[i] = UInt8(s)
                end
            end

            if tron.minor_iter >= max_minor
                tron.status = :Maximum_Iterations_Exceeded
                search = false
            end
        end

        if unsafe_string(pointer(task), 4) == "CONV"
            tron.status = :Solve_Succeeded
            search = false
        end
    end

    return tron.status
end