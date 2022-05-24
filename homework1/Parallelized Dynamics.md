# Problem 1: A Ton of New Facts on Newton
Assumptions:
-	$\frac{dg}{dx}$ is non-singular

## Part 1

*Newton's method update process* $$x_{n+1}=f_1(x_n)=x_n-(\frac{dg}{dx}(x_n))^{-1}g(x_n)$$

**Prove that if $x^∗$ is a steady state of the equation, then $g(x^∗)=0$.**

Suppose $x^*$ is a steady state of Newton's method.
We have 
$$\begin{align}
x^*=f_1(x^*) & \iff x^*=x^*-(\frac{dg}{dx}(x^*))^{-1}g(x^*)\\
& \iff 0 = -(\frac{dg}{dx}(x^*))^{-1}g(x^*)& \\
& \iff 0 = (\frac{dg}{dx}(x^*))(\frac{dg}{dx}(x^*))^{-1} g(x^*) \\
& \iff g(x^*) = 0 \\
\end{align}$$


## Part 2
*Quasi-Newton approximation update process* $$x_{n+1}=f_2(x)=x_n-(\frac{dg}{dx}(x_0))^{-1}g(x_n)$$

**Derive the stability of the Quasi-Newton approximation in the form of a matrix whose eigenvalues need to be constrained.**

First, we compute the Jacobian of $f_2$ 
$$\begin{align}
J(f_2)(x)&=I-(\frac{dg}{dx}(x_0))^{-1}J(g)(x) \\
&=I-(J(g)(x_0))^{-1}J(g)(x) \\
\end{align}$$

If $x^*$ is a steady state, and $J(f_2)(x^*)=I-(J(g)(x_0))^{-1}J(g)(x^*)$ is diagonalizable and has eigenvalues with absolute value less than 1, then the dynamical system $x_{n+1}=f_2(x_n)$ is stable at $x^*$.

**Argue that if $x_0$ is sufficiently close to $x^*$, then the steady state is stable. **


## Part 3
*Relaxed Quasi-Newton approximation update process* $$x_{n+1}=f_{3}(x)=x_n-\alpha(\frac{dg}{dx}(x_0))^{-1}g(x_n)$$
**Argue that for some sufficiently small $\alpha$ that the Quasi-Newton iterations will be stable if the eigenvalues of $((\frac{dg}{dx}(x_0))^{-1}g(x))'$ are all positive for every $x$**

## Part 4
*Fixed-point iteration* $$x_{n+1}=f_{4}(x)=g(x_n)$$

**Fixed-point iteration converges to $x^*=g(x^*)$. What is a small change to the dynamical system that could be done such that $g(x)=0$ is the steady state?**
We can simply subtract $x_n$ or define the dynamical system as $x_{n+1}=x_n-g(x_n)$. The new system converges to $x^{**}$ such that $$x^{**}-g(x^{**})=x^{**} \iff g(x^{**})=0$$

**Change the $(\frac{dg}{dx}(x_0))^{-1}$ term from the Quasi-Newton iteration to get a method equivalent to fixed point iteration.**


# Problem 2
## Part 1
Let $X$ be a random variable, $F$ its CDF and $f$ its PDF. 
The $p$-quantile of $X$ is the number $q$ such that $$F(q) =\mathbb{P}(X\le q)= p$$

Let $p \in [0,1]$  and $g_p(x)=F(x)-p$.
Finding the $p$-quantile is equivalent to solving the root-finding problem $g_p(q)=0$.

*Applying Newton's method*: 
We have $g_p'(x)=F'(x)=f(x)$.
Thus, Newton's method yields the following discretization: $$\begin{align}
q_{n+1}&=q_n-\frac{g_p(q_n)}{g_p'(q_n)} \\
&=q_n-\frac{F(q_n)-p}{f(q_n)} \\
\end{align}$$

## Part 2

# Problem 3
## Part 1
```julia
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
n = 400

using BenchmarkTools
@btime solve_system(logistic, x0, r, n)
@btime solve_system_save(logistic, x0, r, n)

const out = Vector{typeof(x0)}(undef,n)
@btime solve_system_save!(out,logistic, x0, r, n)
@btime calc_attractor!(out,logistic,x0,r)
```
