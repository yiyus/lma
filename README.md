# LMA: Levenberg-Marquardt Algorithm

The `LMA` operator implements the Levenberg-Marquardt algorithm, an iterative procedure widely used for solving non-linear least squares problems or for finding roots of non-linear systems of equations. This implementation is designed to be robust, but at the same time provides sensible defaults to facilitate its usage.

## Introduction

The LMA method interpolates between the more aggressive Gauss-Newton algorithm and the more conservative method of gradient descent.

It is an iterative algorithm. Every iteration, the residual and the Jacobian are calculated for a given set of parameters (the "current guess"). Then a new guess $p-\Delta p$ is calculated, such that:

$$(J^T J + \lambda D)\Delta p = J^T y$$

where $J$ is the Jacobian ($J_{ij}=\partial y_i / \partial p_j$), $D$ is a diagonal matrix such that $D_{kk} = \max(\epsilon(\lambda), (J^T J)_{kk})$ with $\epsilon(\lambda)$ the damping floor (dependent on the damping factor), $y$ is the residual, and $\lambda$ is the damping factor. The damping factor determines how much the next guess approximates the prediction of the Gauss-Newton algorithm (lower damping factors) or the gradient descent methods (higher damping factors), in effect defining a trust-region.

For this new guess, the predicted error reduction is calculated as:

$$\Delta s_p = \frac{1}{2}(\Delta p^T J^T y + \lambda \Delta p^T D \Delta p)$$

If the ratio of the actual error reduction calculated during the next iteration with respect to this prediction is above a defined limit, the current guess is accepted and the damping factor decreased. Else, the new guess is rejected and the damping factor is increased.

The algorithm finishes once one of these conditions is met, returning the last accepted guess:

* Maximum numbers of iterations reached
* Residual below specified tolerance
* Relative change in parameters or residual below specified tolerance
* Maximum damping factor reached

### Normalized damping factor

A particularity of this LMA implementation is the use of a *normalized damping factor*. The damping factor $\lambda$ will oscillate between the specified limits $\lambda_{min}$ and $\lambda_{max}$, starting with an initial value $\lambda_0$. The normalized damping factor is calculated as:

$$\lambda_{norm} = \frac{(\lambda_{max}-\lambda_0)(\lambda-\lambda_{min})}{(\lambda_0-\lambda_{min})(\lambda_{max}-\lambda)}$$

If $\lambda_{norm}$ is 1, the initial damping factor is being used. If it is larger, it means that more damping than indicated was necessary, reaching the maximum $\lambda_{max}$ at infinite; if it is lower, it means that it could be decreased for faster convergence, with the point at 0 corresponding with $\lambda_{min}$.

The usage of normalized damping factors allows to monitor and adjust the effective size of the trust region during successive function calls independently of the actual damping parameters at use.

### Adaptive damping floor for enhance stability

To enhance numerical stability, particularly when dealing with Jacobian matrices $J$ where $J^T J$ may have very small or zero diagonal elements, this LMA implementation employs an adaptive floor for the diagonal elements of the damping scaling matrix $D$. The floor $\epsilon(\lambda)$ is calculated such that it takes the value $\epsilon_0$ (ar arbitrarly small number) when $\lambda$ is $\lambda_0$ or lower than $\lambda_0$ (so, $\lambda_{norm}≤1$), and it approaches 1 as $\lambda$ approaches $\lambda_{max}$.

This adaptive mechanism ensures that while a very low floor is used during optimistic (low damping) phases, a more substantial floor (approaching 1) is automatically applied to weak components when high overall damping is required.

## `LMA` Operator

    R←{X}f LMA Y

`f` is a configuration namespace or a function. Configuration namespaces can define the following configuration parameters (else, the provided defaults are used):

* `toli`: Maximum number of iterations (default `1E6`)
* `tols`: Tolerance for the sum of squared residuals (default `tolr`)
* `tolr`: Tolerance for relative change, either in the solution or the residual (default `⎕CT`)
* `tolg`: Tolerance for the gain ratio to accept or reject a step (default `1E¯2`)
* `dini`: Initial damping factor for `d=1` (default `1E¯2`)
* `dinc`: Increment of damping factor after rejected solution (default `5`)
* `ddec`: Decrement of damping factor after accepted solution (default `÷dinc`)
* `dmax`: Maximum damping factor (default `÷⎕CT`)
* `dmin`: Minimum damping factor (default `÷dmax`)

With the exception of `toli` and `tols`, these parameters should be modified only by expert users.

Additionally, the namespace must contain two functions:

* `Callback`: Callback function (default `⊢`)
* `Eval`: Evaluation function

The evaluation function must return the residual and Jacobian for the set of parameters given as right argument, and an optional left argument `X`.
In case `f` is a function, the default configuration is used, setting the `Eval` function to `f`.
The callback function is called every iteration before checking convergence, with a solution namespace as argument
and discarding the return value.

A solution namespace is a configuration namespace with the additional elements:

* `iter`: Number of iterations
* `sse`: Sum of squared residuals
* `rel`: Relative change metric
* `p0`: Initial guess
* `p`: Current guess

`Y` is a vector. If the depth of the vector is not larger than 1, it is enclosed first.
The first element of `Y` (or `⊂Y` if `1=≡⍵`) contains the initial guess for the parameters.
Additional elements of `Y` must be either configuration namespaces (overwriting the parameters at right with those at left) or a number indicating the value of the initial normalized damping factor.

The returned value `R` is a solution namespace.

#### Note

In addition to being used for the definition of default values, `⎕CT` is also the baseline for adaptive floor damping.


## `LM` Operator

`LM` is a simplified version of `LMA`. Usage:

    R←X f LM g Y

where `f` is the evaluation function, `g` is the callback function (which takes as argument an `iter sse rel dnorm p` vector), `Y` is a two elements vector with the initial guess and normalized damping factor. `X` is a vector defining the parameters `toli tols tolr tolg dini dinc ddec dmax dmin`. The return value `R` will be an `iter sse rel dnorm p` vector.
