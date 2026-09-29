# Posterior covariance and correlation: derivation

## 1. Model and likelihood

We observe $M$ independent trajectories. For path $m = 1,\dots,M$ and time step $k = 0,\dots,N-1$, the log-increments are

$$
\Delta Y_k^{(m)} = \Delta D_k(\theta) + \sigma\sqrt{\Delta t}\;\xi_k^{(m)},
\qquad \xi_k^{(m)} \sim \mathcal{N}(0,1)\ \text{i.i.d.},
$$

with $\theta = (\mu,\lambda)$ and $\Delta D_k(\theta) = D^{\mathrm{sc}}(t_{k+1};\theta) - D^{\mathrm{sc}}(t_k;\theta)$.

The log-likelihood is

$$
\ell(\theta) = \log L(\theta)
= -\frac{1}{2\sigma^2\Delta t}\sum_{m=1}^{M}\sum_{k=0}^{N-1}
\Big(\Delta Y_k^{(m)} - \Delta D_k(\theta)\Big)^2 + \text{const}.
$$

With a flat prior, the posterior is $\pi(\theta\mid\mathcal{D}) \propto e^{\ell(\theta)}$.

---

## 2. Gradient and Hessian of the log-likelihood

Define the residuals $r_k^{(m)}(\theta) = \Delta Y_k^{(m)} - \Delta D_k(\theta)$. Then

$$
\nabla_\theta \ell(\theta)
= \frac{1}{\sigma^2\Delta t}\sum_{m,k} r_k^{(m)}(\theta)\,\nabla_\theta \Delta D_k(\theta),
$$

$$
\nabla_\theta^2 \ell(\theta)
= \frac{1}{\sigma^2\Delta t}\sum_{m,k}
\Big( r_k^{(m)}(\theta)\,\nabla_\theta^2 \Delta D_k(\theta)
\;-\; \nabla_\theta \Delta D_k(\theta)\,\nabla_\theta \Delta D_k(\theta)^{\top} \Big).
$$

At the true parameters $r_k^{(m)} = \sigma\sqrt{\Delta t}\,\xi_k^{(m)}$ has zero mean. So the first term vanishes in expectation, and it is small compared with the second when $\sigma$ is small (the Gauss–Newton / Fisher approximation). Only the second term survives:

$$
-\nabla_\theta^2 \ell(\hat\theta) \;\approx\; F
= \frac{M}{\sigma^2\Delta t}\sum_{k=0}^{N-1}
\nabla_\theta \Delta D_k\,\nabla_\theta \Delta D_k^{\top}
= \frac{1}{\sigma^2}\,G ,
$$

where

$$
\boxed{\;
G = \frac{M}{\Delta t}\sum_{k=0}^{N-1}
\nabla_\theta \Delta D_k\,\nabla_\theta \Delta D_k^{\top}
\;}
$$

depends **only on the forward model**, not on the noise. $F$ is the Fisher information matrix.

---

## 3. Laplace (Gaussian) approximation of the posterior

Expanding $\ell$ to second order around the maximizer $\hat\theta$, where $\nabla\ell(\hat\theta) = 0$:

$$
\ell(\theta) \approx \ell(\hat\theta) - \tfrac12 (\theta-\hat\theta)^{\top} F\, (\theta-\hat\theta).
$$

Hence

$$
\pi(\theta\mid\mathcal{D}) \approx \mathcal{N}\!\left(\hat\theta,\; \Sigma_{\mathrm{post}}\right),
\qquad
\boxed{\;\Sigma_{\mathrm{post}} = F^{-1} = \sigma^2\, G^{-1}\;}
$$

---

## 4. Explicit $2\times 2$ formulas

Write

$$
G = \begin{pmatrix} G_{11} & G_{12} \\ G_{12} & G_{22} \end{pmatrix},
\qquad
\det G = G_{11}G_{22} - G_{12}^2 > 0,
$$

$$
G^{-1} = \frac{1}{\det G}\begin{pmatrix} G_{22} & -G_{12} \\ -G_{12} & G_{11} \end{pmatrix}.
$$

Therefore

$$
\operatorname{Var}(\mu) \approx \sigma^2\,\frac{G_{22}}{\det G},
\qquad
\operatorname{Var}(\lambda) \approx \sigma^2\,\frac{G_{11}}{\det G},
$$

$$
\boxed{\;
\operatorname{Cov}(\mu,\lambda) \approx -\,\sigma^2\,\frac{G_{12}}{\det G}
\;}
$$

and the correlation is

$$
\operatorname{corr}(\mu,\lambda)
= \frac{\operatorname{Cov}(\mu,\lambda)}{\sqrt{\operatorname{Var}(\mu)\operatorname{Var}(\lambda)}}
= \frac{-\sigma^2 G_{12}/\det G}{\sigma^2\sqrt{G_{11}G_{22}}/\det G},
$$

$$
\boxed{\;
\operatorname{corr}(\mu,\lambda) \approx -\,\frac{G_{12}}{\sqrt{G_{11}\,G_{22}}}
\;}
$$

### Consequences

| Quantity | Scaling with noise |
|---|---|
| $\operatorname{std}(\mu),\ \operatorname{std}(\lambda)$ | $\propto \sigma$ |
| $\operatorname{Cov}(\mu,\lambda)$ | $\propto \sigma^2$ |
| $\operatorname{corr}(\mu,\lambda)$ | independent of $\sigma$ |

So a covariance that grows with noise (Fig. 6) does **not** imply stronger coupling. Only the correlation measures coupling.

---

## 5. Interpretation: angle between sensitivity vectors

Our contrast is linear in the amplitude, $\rho = a\, g_\lambda(x)$, so $D^{\mathrm{sc}} = a\,A g_\lambda$. Define the sensitivity vectors (over $k$)

$$
u_k = \frac{\partial \Delta D_k}{\partial a} = \Delta\!\left(A g_\lambda\right)_k,
\qquad
a\,v_k = \frac{\partial \Delta D_k}{\partial \lambda} = a\,\Delta\!\left(A\,\partial_\lambda g_\lambda\right)_k .
$$

Then

$$
G_{11} = \tfrac{M}{\Delta t}\,\|u\|^2,
\qquad
G_{22} = \tfrac{M}{\Delta t}\,a^2\|v\|^2,
\qquad
G_{12} = \tfrac{M}{\Delta t}\,a\,\langle u, v\rangle,
$$

and

$$
\operatorname{corr}(a,\lambda) \approx -\frac{\langle u, v\rangle}{\|u\|\,\|v\|} = -\cos\angle(u,v).
$$

- If the two sensitivities point in **similar** directions ($\langle u,v\rangle > 0$), increasing one parameter can be compensated by **decreasing** the other, which gives a negative correlation.
- If they point in **opposite** directions ($\langle u,v\rangle < 0$), the parameters must move **together**, which gives a positive correlation.
- **Orthogonal** sensitivities mean uncorrelated parameters.

For our forward model, the numerical values are

$$
\operatorname{corr}(a,w) \approx +0.35,
\qquad
\operatorname{corr}(a,p) \approx +0.39,
$$

so $\langle u,v\rangle < 0$ in both cases. For the width: widening the bump weakens the reflections, which come from the contrast gradients, so the change in data from $w$ opposes the change from $a$.

---

## 6. Estimating from MCMC samples

From $N_s$ post-burn-in samples $\theta^{(n)}$:

$$
\widehat{\operatorname{Cov}}(\mu,\lambda)
= \frac{1}{N_s-1}\sum_{n=1}^{N_s}\big(\mu^{(n)}-\bar\mu\big)\big(\lambda^{(n)}-\bar\lambda\big),
\qquad
\hat r = \frac{\widehat{\operatorname{Cov}}(\mu,\lambda)}{\hat s_\mu\,\hat s_\lambda}.
$$

MCMC samples are autocorrelated. The effective sample size is

$$
N_{\mathrm{eff}} = \frac{N_s}{\tau},
\qquad
\tau = 1 + 2\sum_{j\ge 1}\rho_j ,
$$

where $\rho_j$ is the lag-$j$ autocorrelation of the chain. The standard error of the correlation estimate is approximately

$$
\operatorname{SE}(\hat r) \approx \frac{1-r^2}{\sqrt{N_{\mathrm{eff}}}} .
$$

With low acceptance rates, $\tau$ is large and $N_{\mathrm{eff}}$ small, which explains the erratic correlations obtained with the original sampler.

---

## 7. Proposal scaling

The Laplace standard deviations $\sigma\sqrt{(G^{-1})_{ii}}$ give a good random-walk step size (optimal scaling for a $d$-dimensional Gaussian target, here $d = 2$):

$$
s_i = \frac{2.4}{\sqrt{d}}\;\sigma\sqrt{(G^{-1})_{ii}}
\;\approx\; 1.7\,\sigma\sqrt{(G^{-1})_{ii}} ,
$$

which gives acceptance rates of about 20–40%.

---

## 8. Validity

The approximation holds when:

1. the prior is flat (or weak) over the region of high posterior probability,
2. the noise is small enough that $\Delta D_k(\theta)$ is approximately **linear** over the posterior width, which breaks down at large $\sigma$ (bias at $\sigma = 0.8$),
3. $\hat\theta$ is close to the true parameters, so $G$ can be evaluated at $\theta_{\mathrm{true}}$,
4. the parameter is resolved by the discretization. In this model the position is snapped to the grid, so it fails for $p$ once $\operatorname{std}(p) < 1/n \approx 3.3\times10^{-4}$.
