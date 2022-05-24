using StaticArrays
using Distributions, Statistics

function f(q,d,p)
    @inbounds dq1 = q[1]+(p - cdf(d,q[1]))/pdf(d,q[1])
    @SVector [dq1]
end

function solve_system(f, q0, d, p, n) 
    q = copy(q0)
    for i in 1:n-1
        q = f(q,d,p)
    end
    q
end

function myquantile(d::UnivariateDistribution,p::Float64)
    q0 = @SVector [mean(d)]
    res = solve_system(f, q0, d, p, 10)
    res[1]
end

using BenchmarkTools
d = Normal(0,1)
@btime myquantile(d,0.99)
@btime quantile(d,0.99)

d = Gamma(5,1)
@btime myquantile(d,0.90)
@btime quantile(d,0.90)

d = Beta(2, 4)
@btime myquantile(d,0.90)
@btime quantile(d,0.90)

using InteractiveUtils
@code_llvm myquantile(d,0.90)
@code_warntype myquantile(d,0.9)
