import PauliStrings as ps
using LinearAlgebra
using Printf
using ProgressMeter
using Base.Threads
using DelimitedFiles
using ArgParse


"""
XXZ(N, J, Δ) -> OperatorTS1D

Build the 1D XXZ Hamiltonian with periodic boundary conditions on N sites,
coupling strength J, and anisotropy Δ (Pauli basis).
"""
function XXZ(N, J, Δ)
  H = ps.Operator(N)
  i = 1
  H += J, "X", i, "X", mod1(i + 1, N)
  H += J, "Y", i, "Y", mod1(i + 1, N)
  H += J * Δ, "Z", i, "Z", mod1(i + 1, N)

  return ps.OperatorTS1D(H, full=false)
end

"""
hs_product(op1, op2) -> UInt128

Compute the Hilbert–Schmidt inner product tr(op1† * op2) normalized by 2^N.
"""
@inline hs_product(op1::ps.OperatorTS1D, op2::ps.OperatorTS1D) = ps.trace_product(ps.dagger(op1), op2) / UInt128(2)^ps.qubitlength(op1)

"""
hs_norm(op) -> Float64

Compute the Hilbert–Schmidt norm of an operator using hs_product.
"""
@inline hs_norm(op::ps.OperatorTS1D) = sqrt(hs_product(op, op))

"""
undigit(olist::Vector{Int}) -> Int

Convert a base-4 digit list to its integer representation.
"""
@inline undigit(olist) = sum([olist[i] * 4^(i - 1) for i in eachindex(olist)])

symbol_map = Dict(0 => "1", 1 => "X", 2 => "Z", 3 => "Y")
superscript_map = Dict("1" => "1", "X" => "X", "Z" => "Z", "Y" => "Y")

"""
map_with_indices(arr) -> Tuple{Complex,Vararg{Any}}

Build the tuple representation for PauliStrings
"""
function map_with_indices(arr::Vector{String})
  result = Any[]
  for (i, val) in enumerate(arr)
    val != "1" && push!(result, val, i)
  end
  return Tuple((1 + 0im, result...))
end


"""
all_M_local_trans_sym(L, M) -> (Vector{OperatorTS1D}, Vector{Vector{String}})

Generate all M-local Pauli-string operators on L sites (excluding trivial identity).
"""
function all_M_local_trans_sym(L::Int, M::Int)
  ops = ps.OperatorTS1D[]
  ops_list = Vector{String}[]

  for i in 1:4^M-1
    op_list = digits(i, base=4, pad=M)

    # no identity on first site
    op_list[1] == 0 && continue

    # all checks passed, building operator
    symbolic_op = [symbol_map[op] for op in op_list]
    push!(ops_list, symbolic_op)
    op = ps.Operator(L)
    op += map_with_indices(symbolic_op)
    push!(ops, ps.OperatorTS1D(op, full=false))
  end

  return ops, ops_list
end

"""
compute_lioms(H, max_supp) -> (values, vectors, ops_list)

Compute LIOM eigenvalues and eigenvectors by forming the correlation matrix
of commutators [H, op] for all M-local operators up to support max_supp.
"""
function compute_lioms(H::ps.OperatorTS1D, max_supp::Int)
  # Generate operators
  println("Generating operators...")
  ops, ops_list = all_M_local_trans_sym(ps.qubitlength(H), max_supp)
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

  # Set up progress bar
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

"""
parse_args() -> Dict{String,Any}

Parse command-line arguments for the Pauli-basis XXZ LIOM script
(flags: --delta, --max-supp).
"""
function parse_args()
  s = ArgParseSettings(description="Compute LIOMs for 1D XXZ model in Pauli strings basis")
  @add_arg_table s begin
    "--delta", "-d"
    help = "Anisotropy parameter Δ (float)"
    arg_type = Float64
    default = 0.5
    "--max-supp", "-M"
    help = "Maximum support M for LIOMs (integer)"
    arg_type = Int
    default = 3
  end
  return ArgParse.parse_args(s)
end

"""
main()

Entry point: parse CLI args, compute LIOMs, and write results and logs.
"""
function main()
  args = parse_args()

  J = 1.0
  Delta = args["delta"]
  max_supp = args["max-supp"]

  # Parameters
  J = 1.0
  L = 2 * max_supp + 1
  (L % 2 == 1) && (L += 1)
  H = XXZ(L, J, Delta)

  tic = time_ns()
  evals, evecs, ops_list = compute_lioms(H, max_supp)
  toc = time_ns()

  res_path = "$(@__DIR__)/results/1D_XXZ/Pauli_strings"
  if !ispath(res_path)
    mkpath(res_path)
  end
  filename_evals = "$(res_path)/eigenvalues_1D_XXZ_M_$(max_supp)_J_$(J)_d_$(Delta).txt"
  filename_evecs = "$(res_path)/eigenvectors_1D_XXZ_M_$(max_supp)_J_$(J)_d_$(Delta).txt"
  filename_operators = "$(res_path)/operators_1D_XXZ_M_$(max_supp)_J_$(J)_d_$(Delta).txt"
  filename_log = "$(res_path)/log_1D_XXZ_M_$(max_supp)_J_$(J)_d_$(Delta).txt"
  writedlm(filename_evals, evals)
  writedlm(filename_evecs, evecs)
  writedlm(filename_operators, ops_list)
  open(filename_log, "w") do io
    println(io, "J = ", J)
    println(io, "Δ = ", Delta)
    println(io, "max_supp = ", max_supp)
    println(io, "L = ", L)
    println(io, "Operators basis size = ", length(ops_list))
    println(io, "Number of LIOMs found = ", count(x -> abs(x) < 1e-10, evals))
    println(io, "Elapsed time = ", (toc - tic) / 1e9, " s")
  end

  max_evals = min(10, length(evals))
  println("$(max_evals) smallest eigenvalues:")
  println(evals[1:max_evals])
  println("$(max_evals) eigenvectors corrsponding to smallest eigenvalues:")

  for i in eachindex(ops_list)
    op_str = join([superscript_map[s] for s in ops_list[i]], "")
    print("∑_l ", rpad(op_str, max_supp), "   ")

    for j in 1:max_evals
      num_str = @sprintf("% .5f", evecs[i, j])
      num_str = replace(num_str, " -" => "-")
      print(" ", num_str)
    end
    println()
  end
  println("Size of basis for M=$max_supp : $(length(ops_list))")
  println("Number of LIOMs found: ", count(x -> abs(x) < 1e-10, evals))
  println("Elapsed time: ", (toc - tic) / 1e9, " s")

end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end