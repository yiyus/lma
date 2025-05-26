# LMA: Levenberg-Marquardt Algorithm

The `LMA` operator implements the Levenberg-Marquardt algorithm, an iterative procedure widely used for solving non-linear least squares problems or for finding roots of non-linear systems of equations. This implementation is designed to be robust, but at the same time provides sensible defaults to facilitate its usage.

## Introduction

The LMA method interpolates between the more aggressive Gauss-Newton algorithm and the more conservative method of gradient descent.

It is an iterative algorithm. Every iteration, the residual and the Jacobian are calculated for a given set of parameters (the "current guess"). Then a new guess $p-\Delta p$ is calculated, such that:

$$(J^T J + \lambda D)\Delta p = J^T y$$

where $J$ is the Jacobian ($J_{ij}=\partial y_i / \partial p_j$), $D$ is a diagonal matrix such that $D_{kk} = \max(\epsilon, (J^T J)_{kk})$ with $\epsilon$ an arbitrarily small number, $y$ is the residual, and $\lambda$ is the damping factor. The damping factor determines how much the next guess approximates the prediction of the Gauss-Newton algorithm (lower damping factors) or the gradient descent methods (higher damping factors), in effect defining a trust-region.

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

## Usage

See [LMA Jupyter notebook](lma.ipynb).
