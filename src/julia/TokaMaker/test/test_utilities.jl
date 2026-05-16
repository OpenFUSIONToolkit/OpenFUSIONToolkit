# Test infrastructure mirroring test_TokaMaker.py:20-146.
# - mp_run: subprocess isolation with timeout (Distributed.jl worker)
# - validate_dict / validate_eqdsk / validate_ifile

using Distributed
using LinearAlgebra: norm
using Test
using TokaMaker: read_eqdsk, read_ifile

const TEST_DIR = @__DIR__
const REPO_TESTS = abspath(joinpath(TEST_DIR, "..", "..", "..", "tests", "physics"))

"""
    mp_run(target_module, target_func, args; timeout=30)

Run `target_module.target_func(args...)` on a fresh worker. Returns `nothing`
on timeout, otherwise the worker's result. Mirrors `mp_run` in
test_TokaMaker.py.
"""
function mp_run(target_module::Module, target_func::Symbol, args::Tuple;
                timeout::Real=30)
    if get(ENV, "OFT_DEBUG_TEST", "") != ""
        timeout *= 4
    end
    pid = addprocs(1; exeflags="--project=$(Base.active_project())")[1]
    try
        @everywhere [pid] using TokaMaker
        @everywhere [pid] cd($TEST_DIR)
        future = remotecall(getfield(target_module, target_func), pid, args...)
        t0 = time()
        while time() - t0 < timeout
            isready(future) && return fetch(future)
            sleep(0.5)
        end
        @info "mp_run timeout reached after $(timeout)s"
        return nothing
    finally
        rmprocs(pid)
    end
end

"""
    validate_dict(results, expected) -> Bool

1% relative tolerance default; 5% for triangularity (`delta`, `deltaU`,
`deltaL`). Per-element comparison for arrays. Returns true if all keys
match within tolerance, false otherwise (with diagnostics printed).
"""
function validate_dict(results, expected::AbstractDict)
    tol = Dict("delta" => 5e-2, "deltaU" => 5e-2, "deltaL" => 5e-2)
    if results === nothing
        @info "FAILED: error in solve"
        return false
    end
    ok = true
    for (key, exp_val) in expected
        result_val = haskey(results, key) ? results[key] :
                     (results isa NamedTuple && hasproperty(results, Symbol(key)) ?
                      getproperty(results, Symbol(key)) : nothing)
        if result_val === nothing
            @info "FAILED: key \"$key\" missing from results"
            ok = false; continue
        end
        rel = get(tol, key, 1e-2)
        if exp_val isa AbstractArray
            for i in eachindex(exp_val)
                exp_val[i] === nothing && continue
                if abs((result_val[i] - exp_val[i]) / exp_val[i]) > rel
                    @info "FAILED: $key[$i] expected $(exp_val[i]) got $(result_val[i])"
                    ok = false
                end
            end
        else
            if abs((result_val - exp_val) / exp_val) > rel
                @info "FAILED: $key expected $exp_val got $result_val"
                ok = false
            end
        end
    end
    return ok
end

function validate_eqdsk(file_test::AbstractString, file_ref::AbstractString)
    test_data = try
        read_eqdsk(file_test)
    catch e
        @info "FAILED: could not read result EQDSK: $e"
        return false
    end
    ref_data = try
        read_eqdsk(file_ref)
    catch e
        @info "FAILED: could not read reference EQDSK: $e"
        return false
    end
    ok = true
    for (key, exp_val) in ref_data
        key == "case" && continue
        result_val = get(test_data, key, nothing)
        if result_val === nothing
            @info "FAILED: key \"$key\" missing"
            ok = false; continue
        end
        if exp_val isa AbstractArray
            err = norm(exp_val .- result_val) / norm(exp_val)
            if err > 1e-2
                @info "FAILED: $key relative norm error $err"
                ok = false
            end
        else
            if abs((result_val - exp_val) / exp_val) > 1e-2
                @info "FAILED: $key expected $exp_val got $result_val"
                ok = false
            end
        end
    end
    return ok
end

function validate_ifile(ifile_test::AbstractString, ifile_ref::AbstractString)
    test_data = try
        read_ifile(ifile_test)
    catch e
        @info "FAILED: could not read result i-file: $e"
        return false
    end
    ref_data = try
        read_ifile(ifile_ref)
    catch e
        @info "FAILED: could not read reference i-file: $e"
        return false
    end
    ok = true
    for (key, exp_val) in ref_data
        result_val = get(test_data, key, nothing)
        if result_val === nothing
            @info "FAILED: key \"$key\" missing"
            ok = false; continue
        end
        if exp_val isa AbstractArray
            err = norm(exp_val .- result_val) / norm(exp_val)
            if err > 1e-2
                @info "FAILED: $key relative norm error $err"
                ok = false
            end
        else
            if result_val != exp_val
                @info "FAILED: $key expected $exp_val got $result_val"
                ok = false
            end
        end
    end
    return ok
end

# Helper for marker-based filtering. Tests can be tagged with strings;
# the test runner consults ENV["OFT_TEST_FILTER"] (e.g. "coverage", "slow").
const TEST_FILTER = get(ENV, "OFT_TEST_FILTER", "")

should_run(tags::AbstractVector{<:AbstractString}) = isempty(TEST_FILTER) || (TEST_FILTER in tags)
