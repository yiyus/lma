# LMA: Levenberg-Marquardt Algorithm

The `LMA` operator implements the Levenberg-Marquardt algorithm, an iterative procedure widely used for solving non-linear least squares problems or for finding roots of non-linear systems of equations. This implementation is designed to be robust, but at the same time provides sensible defaults to facilitate its usage.

## Introduction

The LMA method interpolates between the more aggressive Gauss-Newton algorithm and the more conservative method of gradient descent.

It is an iterative algorithm. Every iteration, the residual and the Jacobian are calculated for a given set of parameters $p$. Then a new guess $p-\Delta p$ is calculated, such that:

$$(J^T J + \lambda D)\Delta p = J^T y$$

where $J$ is the Jacobian ($J_{ij}=\partial y_i / \partial p_j$), $D$ is a diagonal matrix such that $D_{kk} = \max(\epsilon(\lambda), (J^T J)_{kk})$ with $\epsilon(\lambda)$ the damping floor (dependent on the damping factor), $y$ is the residual, and $\lambda$ is the damping factor. The damping factor determines how much the next guess approximates the prediction of the Gauss-Newton algorithm (lower damping factors) or the gradient descent methods (higher damping factors), in effect defining a trust-region.

For this new guess, the predicted error reduction is calculated as:

$$\Delta s_p = \frac{1}{2}(\Delta p^T J^T y + \lambda \Delta p^T D \Delta p)$$

If the ratio of the actual error reduction calculated during the next iteration with respect to this prediction is above a defined limit, the current guess is accepted and the damping factor is decreased. Else, the new guess is rejected and the damping factor is increased.

The algorithm finishes once one of these conditions is met, returning the last accepted guess:

* Maximum numbers of iterations reached
* Residual below specified tolerance
* Relative change in parameters or residual below specified tolerance
* Stagnation with maximum damping factor

### Normalized damping factor

A particularity of this LMA implementation is the use of a *normalized damping factor*. The damping factor $\lambda$ will oscillate between the specified limits $\lambda_{min}$ and $\lambda_{max}$, starting with an initial value $\lambda_0$. The normalized damping factor is calculated as:

$$\lambda_{norm} = \frac{(\lambda_{max}-\lambda_0)(\lambda-\lambda_{min})}{(\lambda_0-\lambda_{min})(\lambda_{max}-\lambda)}$$

If $\lambda_{norm}$ is 1, the initial damping factor is being used. If it is larger, it means that more damping than indicated was necessary, reaching the maximum $\lambda_{max}$ at infinite; if it is lower, it means that it could be decreased for faster convergence, with the point at 0 corresponding with $\lambda_{min}$.

The usage of normalized damping factors allows to monitor and adjust the effective size of the trust region during successive function calls independently of the actual damping parameters at use.

### Adaptive damping floor

To enhance numerical stability, particularly when dealing with Jacobian matrices $J$ where $J^T J$ may have very small or zero diagonal elements, this LMA implementation employs an adaptive floor for the diagonal elements of the damping scaling matrix $D$. The floor $\epsilon(\lambda)$ is calculated such that it takes the value $\epsilon_0$ (an arbitrarily small number) when $\lambda$ is $\lambda_0$ or lower than $\lambda_0$ (so, $\lambda_{norm}≤1$), and it approaches 1 as $\lambda$ approaches $\lambda_{max}$.

This adaptive mechanism ensures that while a very low floor is used during optimistic (low damping) phases, a more substantial floor (approaching 1) is automatically applied to weak components when high overall damping is required.

## Usage

### `LMA` Operator

    R←{X}f LMA Y

`LMA` is a monadic operator, which takes a left operand to return a derived ambivalent function. This derived function allows to minimize a residual function with a known Jacobian using the Levenberg-Marquadt algorithm, given an initial set of parameters. Several configuration options are available, with sensible defaults previously defined.

The left operand `f` must be a configuration namespace or a function. Configuration namespaces may define the following configuration parameters:

* `toli`: Maximum number of iterations (default `1E3`)
* `tols`: Tolerance for the sum of squared residuals (default `⎕CT`)
* `tolr`: Tolerance for relative change, either in the solution or the residual (default `⎕CT`)
* `tolg`: Tolerance for the gain ratio to accept or reject a step (default `1E¯2`)
* `dini`: Initial damping factor for `dnorm=1` (default `1E¯2`)
* `dinc`: Increment of damping factor after rejected solution (default `5`)
* `ddec`: Decrement of damping factor after accepted solution (default `÷dinc`)
* `dmax`: Maximum damping factor (default `÷⎕CT`)
* `dmin`: Minimum damping factor (default `÷dmax`)
* `dp`: Perturbation applied to parameters for numerical estimation of the Jacobian (default `⎕CT`)
* `verbose`: If `1`, print `iter ssr rel dnorm p` each iteration (default `0`)

Configuration namespaces may also contain the functions:

* `Callback`: Callback function (default `⊢`)
* `Eval`: Evaluation function

The evaluation function `Eval` must return either the residuals and the Jacobian for the given set of parameters, or only the residuals. Whenever the residual and Jacobian need to be evaluated, the function `Eval` will be called with trial parameters as right argument and left argument `X`, if given (`Eval` will be called monadically if the derived function `f LMA` is called monadically). `Eval` must return either a two elements vector with the residuals in the first element and the Jacobian in the second one, or a vector of residuals, optionally enclosed. If a Jacobian is not returned, a numerical estimation is calculated.

The `Callback` function will be called every iteration before checking convergence, with a solution namespace for the current guess as right argument and `X` as left argument, if given (`Callback` will be called monadically if the derived function `f LMA` is called monadically). Its return value is discarded.

If `f` is a function, the result is equivalent to using as `f` a namespace with an `Eval` fuction `f`
(with default values for the rest of parameters).

`Y` must be a vector.
The first element of `Y`, or `⊂Y` if `1=≡⍵`, contains the initial guess for the parameters.
If the next element of `Y` is a scalar numeric value, it is interpreted as the initial normalized damping factor.
Additional elements of `Y` must be configuration namespaces. The final configuration parameters are obtained
overwriting the parameters in the namespace given as left operand with those given as right argument from right to left. Default values will be used for non-defined parameters and the `Callback` function, but the `Eval` function must be defined by the user either as left operand `f` or as member of a configuration namespace.

The returned value `R` is a solution namespace corresponding to the last accepted solution.
A solution namespace is a configuration namespace with all the parameters used when running the algorithm and the additional elements:

* `iter`: Number of iterations
* `ssr`: Sum of squared residuals
* `rel`: Relative change metric
* `dnorm`: Normalized damping factor
* `p0`: Initial guess
* `p`: Accepted guess

#### Notes

* With the exception of `toli` and `tols`, configuration parameters should be modified only by expert users or in case of convergence problems

* The relative change metric `rel` is the minimum relative change between successive accepted solutions either in the the sum of squared residuals or in the parameters themselves

* In addition to being used for the definition of default values, `⎕CT` is also the baseline for adaptive floor damping


### `LM` Operator

`LM` is a simplified version of `LMA`. It is a dyadic operator from which a dydadic function is derived. Usage:

    R←X f LM g Y

where `f` is a monadic function which returns the residuals and Jacobian given a set of parameters, `g` is a monadic function which takes as argument an `iter ssr rel dnorm p` vector and gets called before every convergence check, `Y` is a two elements vector with the initial guess of parameters and normalized damping factor, and `X` is a vector with the configuration parameters `toli tols tolr tolg dini dinc ddec dmax dmin`. The return value `R` will be an `iter ssr rel dnorm p` vector.
