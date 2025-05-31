# LMA: Levenberg-Marquardt Algorithm

The Levenberg-Marquardt algorithm is an iterative procedure widely used for solving non-linear least squares problems or for finding roots of non-linear systems of equations. This implementation is designed to be robust and offer maximum flexibility, but at the same time provides sensible defaults to facilitate its usage.

## Introduction

The LMA method interpolates between the more aggressive Gauss-Newton algorithm and the more conservative method of gradient descent to find the set of *solution parameters* $p$ that minimize a cost function $F(p)$, defined as:

$$F(p) = \sum_i \rho(y_i(p), c_i)$$

where $y_i(p)$ are the components of the residuals vector $y$, and $\rho$ is a loss function which depends on the residuals and a scaling factor $c_i$. In the standard case ($L_2$ loss function), $F(p)$ is the sum of squared residuals (and $c_i = 1$).

Every iteration of the algorithm, the residuals and the Jacobian are calculated for a given set of parameters $p$. Then, a new guess $p-\Delta p$ is calculated such that:

$$(J^T W J + \lambda D)\Delta p = J^T W y$$

where $J$ is the Jacobian of the residuals ($J_{ij}=\partial y_i / \partial p_j$), $W$ is a weight matrix depending on the choice of loss function (the identity matrix for $L_2$), $D$ is a diagonal matrix such that $D_{kk} = \max(\epsilon(\lambda), (J^T W J)_{kk})$ with $\epsilon(\lambda)$ the damping floor (dependent on the damping factor), $y$ is the residuals vector, and $\lambda$ is the damping factor. The damping factor determines how much the next guess approximates the prediction of the Gauss-Newton algorithm (lower damping factors) or the gradient descent methods (higher damping factors), in effect defining a trust-region.

For this new guess, the predicted error reduction is calculated as:

$$\Delta s_p = \frac{1}{2}(\Delta p^T J^T W y + \lambda \Delta p^T D \Delta p)$$

If the ratio of the actual error reduction calculated during the next iteration with respect to this prediction is above a defined limit, the current guess is accepted and the damping factor is decreased. Else, the new guess is rejected and the damping factor is increased.

The algorithm finishes once one of these conditions is met, returning the last accepted guess:

* Maximum number of iterations reached
* Cost (sum of losses) below specified tolerance
* Relative change in solution parameters or cost below specified tolerance
* Stagnation at maximum damping factor

### Choice of loss function

The standard $L_2$ loss function calculates the loss as the square of the residual. Its corresponding weight matrix is the identity. Moreover, if an individual scaling factor is defined for each residual, this loss function can be used for the solution of *weighted least-squares* problems.

Robust loss functions provide mechanisms to mitigate the effect of outliers in fitting data:

* **Huber** Equivalent to $L_2$ for small residuals and linear for larger residuals
* **Cauchy** Strongly down-weights large outliers
* **L1-Soft** Smooth approximation of $L_1$ loss that behaves like $L_2$ for small residuals
* **Tukey** Redescending M-estimator that completely rejects extreme outliers (but it may lead to convergence issues if scaling is not chosen well)
* **Welsh** Another redescending M-estimator, smoother than Tukey's in its rejection
* **Fair** Less sensitive to large errors than $L_2$, but not redescending
* **Arctan** Limits maximum loss of single residuals

### Normalized damping factor

A particularity of this LMA implementation is the use of a *normalized damping factor*. The damping factor $\lambda$ will oscillate between the specified limits $\lambda_{min}$ and $\lambda_{max}$, starting with an initial value $\lambda_0$. The normalized damping factor is calculated as:

$$\lambda_{norm} = \frac{(\lambda_{max}-\lambda_0)(\lambda-\lambda_{min})}{(\lambda_0-\lambda_{min})(\lambda_{max}-\lambda)}$$

If $\lambda_{norm}$ is 1, the initial damping factor is being used. If it is larger, it means that more damping than indicated was necessary, reaching the maximum $\lambda_{max}$ at infinite; if it is lower, it means that it could be decreased for faster convergence, with the point at 0 corresponding with $\lambda_{min}$.

The usage of normalized damping factors allows to monitor and adjust the effective size of the trust region during successive function calls independently of the actual damping parameters at use.

### Adaptive damping floor

To enhance numerical stability, particularly when dealing with Jacobian matrices $J$ where $J^T W J$ may have very small or zero diagonal elements, this LMA implementation employs an adaptive floor for the diagonal elements of the damping scaling matrix $D$. The floor $\epsilon(\lambda)$ is calculated such that it takes the value $\epsilon_0$ (an arbitrarily small number) when $\lambda$ is $\lambda_0$ or lower than $\lambda_0$ (so, $\lambda_{norm}≤1$), and it approaches 1 as $\lambda$ approaches $\lambda_{max}$.

This adaptive mechanism ensures that while a very low floor is used during optimistic (low damping) phases, a more substantial floor (approaching 1) is automatically applied to weak components when high overall damping is required.

## Usage

### ` Jacobian` Operator

    R←{X}f Jacobian Y
    
`Jacobian` is a monadic operator that takes a monadic function `f` as left operand to return an ambivalent function. This derived function returns an estimation of the Jacobian matrix of `f`, using the method of finite differences. The right argument `Y` is the value at which the Jacobian is calculated, and the optional left argument `X` is the relative perturbation to apply to `Y` in the finite differences method. If `X` is not a given, `⎕CT*÷2` is used.

### `LMA` Operator

    R←{X}f LMA Y

`LMA` is a monadic operator, which takes a left operand to return a derived ambivalent function. This derived function allows to minimize a residual function with a known Jacobian using the Levenberg-Marquardt algorithm, given an initial set of parameters. Several configuration options are available, with sensible defaults previously defined.

The left operand `f` must be a configuration namespace or a function. Configuration namespaces may define the following configuration options:

* `toli`: Maximum number of iterations (default `1E3`)
* `tolc`: Tolerance for the cost (sum of squared residuals or loss values) (default `⎕CT`)
* `tolr`: Tolerance for relative change, either in the solution or the residual (default `⎕CT`)
* `tolg`: Tolerance for the gain ratio to accept or reject a step (default `1E¯2`)
* `dini`: Initial damping factor for `dnorm=1` (default `1E¯2`)
* `dinc`: Increment of damping factor after rejected solution (default `5`)
* `ddec`: Decrement of damping factor after accepted solution (default `÷dinc`)
* `dmax`: Maximum damping factor (default `÷⎕CT`)
* `dmin`: Minimum damping factor (default `÷dmax`)
* `pert`: Relative perturbation applied to parameters for numerical estimation of the Jacobian (default `⎕CT*÷2`)
* `loss`: Choice of loss function: `L2` `Huber` `Cauchy` `L1Soft` `Tukey` `Welsh` `Fair` `Arctan` or dyadic function (default `L2`)
* `scale`: Scale factor passed as left argument to loss function (default for 95% efficiency in robust loss functions)
* `verbose`: If `1`, print `iter cost rel dnorm p` each iteration (default `0`)

Configuration namespaces may also contain the functions:

* `Callback`: Callback function (default `⊢`)
* `Eval`: Evaluation function

The evaluation function `Eval` must return either the residuals and the Jacobian for the given set of solution parameters, or only the residuals. Whenever the residual and Jacobian need to be evaluated, the function `Eval` will be called with trial parameters as right argument and left argument `X`, if given (`Eval` will be called monadically if the derived function `f LMA` is called monadically). `Eval` must return either a two elements vector with the residuals in the first element and the Jacobian in the second one, or a vector of residuals, enclosed if they are not simple scalars. If a Jacobian is not returned, a numerical estimation is calculated evaluating the residual function after applying small perturbations to the parameters (as defined by `pert`).

The function selected by the option `loss` is used to calculate the loss from the residuals and scaling factor. If a function is provided by the user, it must be a dyadic function which returns the loss values and weights when given the residuals as right argument and scaling factor as left argument.

The `Callback` function will be called every iteration before checking convergence, with the current solution namespace as right argument and `X` as left argument, if given (`Callback` will be called monadically if the derived function `f LMA` is called monadically). Its return value is discarded.

If `f` is a function, the result is equivalent to using as `f` a namespace with an `Eval` function `f`
(with default values for the rest of parameters).

`Y` must be a vector.
The first element of `Y`, or `⊂Y` if `1=≡Y`, contains the initial guess for the solution parameters.
If the next element of `Y` is a scalar numeric value, it is interpreted as the initial normalized damping factor.
Additional elements of `Y` must be configuration namespaces. The final configuration parameters are obtained
overwriting the parameters in the namespace given as left operand with those given as right argument from right to left. Default values will be used for non-defined parameters and the `Callback` function, but the `Eval` function must be defined by the user either as left operand `f` or as member of a configuration namespace.

The returned value `R` is a solution namespace corresponding to the last accepted solution.
A solution namespace is a configuration namespace including all the configuration options used to run the algorithm and the additional elements:

* `iter`: Number of iterations
* `cost`: Sum of loss values (squared residuals for L2)
* `rel`: Relative change metric
* `dnorm`: Normalized damping factor
* `p0`: Initial guess
* `p`: Accepted guess

#### Notes

* With the exception of `toli` and `tolc`, configuration parameters should be modified only by expert users or in case of convergence problems

* The relative change metric `rel` is the minimum relative change between successive accepted solutions either in the the cost or in the solution parameters

* In addition to being used for the definition of default values, `⎕CT` is also the baseline for adaptive floor damping

* The perturbation to estimate the Jacobian `pert` and the scaling factor for loss functions `scale` can be either scalar values, or vectors of the same length of respectively the parameters and the residuals

* Loss functions and their respective weights, as well as their default values for the scaling parameter are defined in the namespace `Loss`


### `LM` Operator

`LM` is a simplified version of `LMA`. It is a dyadic operator from which a dyadic function is derived. Usage:

    R←X f LM g Y

where `f` is a monadic evaluation function, `g` is a monadic function which takes as argument an `iter cost rel dnorm p` vector and gets called before every convergence check, `Y` is a two elements vector with the initial guess of parameters and normalized damping factor, and `X` is a vector with the configuration parameters `toli tolc tolr tolg dini dinc ddec dmax dmin`. The function `f` must return either the residuals (as a simple or nested vector), the residuals and the Jacobian, or the residuals, Jacobian, loss values and weights, for a set of parameters. If no Jacobian is provided, it is estimated numerically using `Jacobian` (with no left argument). If no loss values and weights are provided, squared residuals are used.

The return value `R` is an `iter cost rel dnorm p` vector.

### Tutorial

See [Jupyter notebook](lma.ipynb).

