function logistic(x, r)
    dx = r * x * (1-x)
    dx
end

function solve_system(f, x0, p, n) 
    x = copy(x0)
    for i in 1:n-1
        x = f(x, p)
    end
    x
end

function solve_system_save(f, x0, p, n) 
    x = Vector{typeof(x0)}(undef,n)
    @inbounds x[1] = x0
    @inbounds for i in 1:n-1
        x[i+1] = f(x[i],p)
    end
    x
end

function solve_system_save!(out, f, x0, p, n) 
    @inbounds out[1] = x0
    @inbounds for i in 1:n-1
        out[i+1] = f(out[i],p)
    end
end

function calc_attractor!(out,f,x0,p,num_attract=150,warmup=400)
    attractor = solve_system(f,x0,p,warmup)
    solve_system_save!(out, f, attractor, p, num_attract)
end

r = 2.9
x0 = 0.25
warmup = 400
n = 150

using BenchmarkTools
@btime solve_system(logistic, x0, r, warmup)
@btime solve_system_save(logistic, x0, r, warmup)

const out = Vector{typeof(x0)}(undef,n)
@btime solve_system_save!(out,logistic, x0, r, n)
@btime calc_attractor!(out,logistic,x0,r)

# Part 2: Bifurcation plot
using StaticArrays
r = 2.9:0.001:4
m = length(r)
result = @SMatrix{m,m}(undef)

using Plots
xs = solve_system_save(logistic, x0, r, 400)
ts = [i for i in 1:length(to_plot)]
plot(ts, xs)
