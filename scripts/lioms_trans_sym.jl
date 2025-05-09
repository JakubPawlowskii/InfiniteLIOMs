import PauliStrings as ps
using LinearAlgebra
using Printf
using ProgressMeter
using Base.Threads
using DelimitedFiles
using ArgParse

"""
XXZ(N, J, Δ) -> OperatorTS1D

Build the 1D XXZ Hamiltonian with periodic boundary conditions for N sites,
coupling strength J, and anisotropy Δ.
"""
function XXZ(N, J, Δ)
  H = ps.Operator(N)
  i = 1
  H += J / 2, "S+", i, "S-", mod1(i + 1, N)
  H += J / 2, "S-", i, "S+", mod1(i + 1, N)
  H += J * Δ, "Sz", i, "Sz", mod1(i + 1, N)
  return ps.OperatorTS1D(H, full=false)
end

"""
hs_product(op1, op2) -> UInt128

Compute the Hilbert–Schmidt inner product tr(op1† * op2) normalized by 2^N.
"""
@inline hs_product(op1::ps.OperatorTS1D, op2::ps.OperatorTS1D) = ps.trace_product(ps.dagger(op1), op2) / UInt128(2)^ps.qubitlength(op1)

"""
hs_norm(op) -> Float64

Compute the Hilbert–Schmidt norm of an operator.
"""
@inline hs_norm(op::ps.OperatorTS1D) = sqrt(hs_product(op, op))

"""
hc(o) -> Int

Return the Hermitian conjugate: 1↔3, leaving other labels unchanged.
"""
@inline hc(o) = (o == 1 ? 3 : (o == 3 ? 1 : o))

"""
undigit(olist::Vector{Int}) -> Int

Convert a base-4 digit list to its integer representation.
"""
@inline undigit(olist) = sum([olist[i] * 4^(i - 1) for i in eachindex(olist)])

symbol_map = Dict(0 => "1", 1 => "S+", 2 => "Sz", 3 => "S-")
superscript_map = Dict("S+" => "+", "Sz" => "Z", "S-" => "-", "1" => "1")

"""
map_with_indices(arr) -> Tuple

Convert a list of Pauli strings and indices into the ps.Operator representation tuple.
"""
function map_with_indices(arr::Vector{String})
  result = Any[]
  for (i, val) in enumerate(arr)
    val != "1" && push!(result, val, i)
  end
  return Tuple((1 + 0im, result...))
end

"""
is_odd_under_parity(op_list, time_reversal) -> Bool

Determine if an operator list is odd under parity assuming given time-reversal symmetry.
"""
function is_odd_under_parity(op_list::Vector{Int}, time_reversal::Symbol)
  cnt_z = count(x -> x == 2, op_list)
  sign = (-1)^cnt_z
  if time_reversal == :odd
    return sign == 1
  end
  return sign == -1
end

"""
is_even_under_parity(op_list, time_reversal) -> Bool

Determine if an operator list is even under parity assuming given time-reversal symmetry.
"""
function is_even_under_parity(op_list::Vector{Int}, time_reversal::Symbol)
  return !is_odd_under_parity(op_list, time_reversal)
end

"""
all_M_local_trans_sym(L, M, time_reversal, parity, conserve_Sz) -> (Vector{OperatorTS1D}, Vector)

Generate all M-local operators on L sites that satisfy the given time-reversal,
parity, and Sz conservation constraints.
"""
function all_M_local_trans_sym(L::Int, M::Int, time_reversal::Symbol, parity::Symbol, conserve_Sz::Symbol)
  ops = ps.OperatorTS1D[]
  ops_list = []

  for i in 1:4^M-1

    # convert to base 4
    op_list = digits(i, base=4, pad=M)

    # no identity on first site
    op_list[1] == 0 && continue

    # commutes with Sz check
    Sz_check = sum(op_list .== 1) - sum(op_list .== 3)
    (conserve_Sz == :yes) && Sz_check != 0 && continue
    (conserve_Sz == :no) && Sz_check == 0 && continue

    # hermitian conjugate check
    op_list_hc = hc.(op_list)
    op_hc = undigit(op_list_hc)
    op_hc < i && continue

    # all checks passed, building operator
    symbolic_op = [symbol_map[op] for op in op_list]
    op = ps.Operator(L)
    op += map_with_indices(symbolic_op)

    # consider op odd under time-reversal
    if time_reversal == :odd || time_reversal == :both

      skip_op = false
      if parity == :even
        skip_op = !is_even_under_parity(op_list, :odd)
      elseif parity == :odd
        skip_op = !is_odd_under_parity(op_list, :odd)
      end

    # cannot have operators odd under time-reversal with only Sz and Id
      if sum(in.(op_list, Ref([1, 3]))) == 0
        skip_op = true
      end

      if !skip_op
        push!(ops_list, (symbolic_op, :odd))
        op -= ps.dagger(op)
        push!(ops, ps.OperatorTS1D(op, full=false))
      end
    end

    # consider op even under time-reversal
    if time_reversal == :even || time_reversal == :both
      skip_op = false
      if parity == :even
        skip_op = !is_even_under_parity(op_list, :even)
      elseif parity == :odd
        skip_op = !is_odd_under_parity(op_list, :even)
      end
      if !skip_op
        push!(ops_list, (symbolic_op, :even))
        op += ps.dagger(op)
        push!(ops, ps.OperatorTS1D(op, full=false))
      end
    end
  end

  return ops, ops_list
end

"""
compute_lioms(H, max_supp, time_reversal, parity, conserve_Sz) -> (values, vectors, ops_list)

Compute LIOMs by building the correlation matrix
of commutators with Hamiltonian H up to support max_supp under given symmetries.
"""
function compute_lioms(H::ps.OperatorTS1D, max_supp::Int, time_reversal::Symbol, parity::Symbol, conserve_Sz::Symbol)
  # Generate operators
  println("Generating operators...")
  ops, ops_list = all_M_local_trans_sym(ps.qubitlength(H), max_supp, time_reversal, parity, conserve_Sz)
  println("Generated $(length(ops)) operators")
  println("Computing norms...")
  ops ./= hs_norm.(ops)

  # Precompute commutators
  println("Precomputing commutators...")
  comms = Vector{ps.OperatorTS1D}(undef, length(ops))
  @showprogress for i in eachindex(ops)
    comms[i] = im * ps.com(H, ops[i])
  end

  # Initialize correlation matrix
  corr_mat = zeros(Float64, length(ops), length(ops))

  total = (length(ops) * (length(ops) + 1)) ÷ 2
  p = Progress(total; desc="Computing correlation matrix...", showspeed=true)

  @threads :greedy for i in eachindex(ops)
    comm_i = comms[i]
    for j in i:length(ops)
      comm_j = comms[j]
      corr_mat[i, j] = hs_product(comm_i, comm_j)
      i != j && (corr_mat[j, i] = corr_mat[i, j])
      next!(p)
    end
  end
  finish!(p)

  F = eigen(Symmetric(corr_mat))
  return F.values, F.vectors, ops_list
end

function parse_args()
  s = ArgParseSettings(description="Compute LIOMs for 1D XXZ model in symmetry resolved basis")
  @add_arg_table s begin
    "--delta", "-d"
    help = "Anisotropy parameter Δ (float)"
    arg_type = Float64
    default = 0.5
    "--max-supp", "-M"
    help = "Maximum support M for LIOMs (integer)"
    arg_type = Int
    default = 3
    "--time-reversal", "-T"
    help = "Behavior under time-reversal symmetry (even, odd, both)"
    arg_type = String
    default = "odd"
    "--parity", "-P"
    help = "Behavior under parity (even, odd, both)"
    arg_type = String
    default = "even"
    "--conserve-Sz", "-C"
    help = "Conserve total Sz? (yes, no, both)"
    arg_type = String
    default = "yes"
  end
  return ArgParse.parse_args(s)
end

function main()
  args = parse_args()
  allowed_time_rev = ["even", "odd", "both"]
  !(args["time-reversal"] in allowed_time_rev) && error("Invalid --type: $(args["time-reversal"]). Must be one of $(allowed_time_rev)")
  allowed_parity = ["even", "odd", "both"]
  !(args["parity"] in allowed_parity) && error("Invalid --parity: $(args["parity"]). Must be one of $(allowed_parity)")
  allowed_conserve_Sz = ["yes", "no", "both"]
  !(args["conserve-Sz"] in allowed_conserve_Sz) && error("Invalid --conserve-Sz: $(args["conserve-Sz"]). Must be one of $(allowed_conserve_Sz)")

  J = 1.0
  Delta = args["delta"]
  max_supp = args["max-supp"]
  time_reversal = Symbol(args["time-reversal"])
  conserve_Sz = Symbol(args["conserve-Sz"])
  parity = Symbol(args["parity"])

  L = 2 * max_supp + 1
  (L % 2 == 1) && (L += 1)

  H = XXZ(L, J, Delta)
  tic = time_ns()
  evals, evecs, ops_list = compute_lioms(H, max_supp, time_reversal, parity, conserve_Sz)
  toc = time_ns()

  res_path = "$(@__DIR__)/results/1D_XXZ/symmetry_resolved/"
  if !ispath(res_path)
    mkpath(res_path)
  end
  filename_evals = "$(res_path)/eigenvalues_1D_XXZ_M_$(max_supp)_J_$(J)_d_$(Delta)_T_$(time_reversal)_P_$(parity)_Sz_cons_$(conserve_Sz).txt"
  filename_evecs = "$(res_path)/eigenvectors_1D_XXZ_M_$(max_supp)_J_$(J)_d_$(Delta)_T_$(time_reversal)_P_$(parity)_Sz_cons_$(conserve_Sz).txt"
  filename_operators = "$(res_path)/operators_1D_XXZ_M_$(max_supp)_J_$(J)_d_$(Delta)_T_$(time_reversal)_P_$(parity)_Sz_cons_$(conserve_Sz).txt"
  filename_log = "$(res_path)/log_1D_XXZ_M_$(max_supp)_J_$(J)_d_$(Delta)_T_$(time_reversal)_P_$(parity)_Sz_cons_$(conserve_Sz).txt"
  writedlm(filename_evals, evals)
  writedlm(filename_evecs, evecs)
  writedlm(filename_operators, [[time_rev, [superscript_map[s] for s in op]...] for (op, time_rev) in ops_list])
  open(filename_log, "w") do io
    println(io, "J = ", J)
    println(io, "Δ = ", Delta)
    println(io, "max_supp = ", max_supp)
    println(io, "Behavior under time-reversal symmetry = ", time_reversal)
    println(io, "Behavior under parity = ", parity)
    println(io, "Conserve Sz? = ", conserve_Sz)
    println(io, "Operators basis size = ", length(ops_list))
    println(io, "Number of LIOMs found = ", count(x -> abs(x) < 1e-10, evals))
    println(io, "Elapsed time = ", (toc - tic) / 1e9, " s")
  end

  max_evals = min(10, length(evals))
  println("$(max_evals) smallest eigenvalues:")
  println(evals[1:max_evals])
  println("$(max_evals) eigenvectors corrsponding to smallest eigenvalues:")

  for i in eachindex(ops_list)
    op_str = join([superscript_map[s] for s in ops_list[i][1]], "")
    op_str_hc = join([hc(superscript_map[s]) for s in ops_list[i][1]], "")
    time_rev = ops_list[i][2]
    op_display = op_str
    op_hc_display = op_str_hc

    if time_rev == :odd
      label = "∑_l i("
      sep = " - "
      tail = ")"
      total_width = max_supp * 2 + length(label) + length(sep) + length(tail)
      content = label * lpad(op_display, max_supp) * sep * lpad(op_hc_display, max_supp) * tail
    elseif time_rev == :even
      label = "∑_l "
      sep = " + "
      total_width = max_supp * 2 + length(label) + length(sep)
      content = label * lpad(op_display, max_supp) * sep * lpad(op_hc_display, max_supp)
    else
      label = "∑_l "
      total_width = max_supp + length(label)
      content = label * lpad(op_display, max_supp)
    end

    print(rpad(content, 2 * max_supp + 12)) 

    for j in 1:max_evals
      num_str = @sprintf("% .5f", evecs[i, j])
      num_str = replace(num_str, " -" => "-")
      print(" ", num_str)
    end
    println()
  end
  println("Size of basis for M=$max_supp, restricted to time reversal = $time_reversal, parity = $parity, Sz conservation = $conserve_Sz: ", length(ops_list))
  println("Number of LIOMs found: ", count(x -> abs(x) < 1e-10, evals))
  println("Elapsed time: ", (toc - tic) / 1e9, " s")
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end