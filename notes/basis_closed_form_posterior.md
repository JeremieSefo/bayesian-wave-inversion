# Basis expansion: the posterior in closed form (no MCMC)

In the basis case the posterior of the coefficients is **exactly Gaussian**, with a mean and a covariance given by explicit formulas. Uncertainty is then quantified by drawing **independent** samples from that Gaussian, or even analytically, with no Markov chain. This note derives it step by step and says where each step appears in `regularized_inversion.m`.

---

## 1. Model

**Basis expansion of the contrast** (hat or Gaussian functions $\phi_j$ on $[0, 0.8]$):

$$
\rho(x) = \sum_{j=1}^{N_b} \beta_j\,\phi_j(x)
\qquad\Longleftrightarrow\qquad
\boldsymbol\rho = \Phi\,\boldsymbol\beta,
\quad \Phi \in \mathbb{R}^{n\times N_b}.
$$

**Linear forward map.** With $J = \Delta A$ (row $k$ of $J$ is $A_{k+1,:} - A_{k,:}$), the noise-free data increments are

$$
\Delta D^{\mathrm{sc}} = J\boldsymbol\rho = J\Phi\,\boldsymbol\beta =: B\,\boldsymbol\beta,
\qquad B \in \mathbb{R}^{(N-1)\times N_b}.
$$

**Noise.** For each of the $M$ trajectories the log-increments are

$$
\mathbf y^{(m)} = B\boldsymbol\beta + \sigma\sqrt{\Delta t}\;\boldsymbol\xi^{(m)},
\qquad \boldsymbol\xi^{(m)} \sim \mathcal N(0, I)\ \text{i.i.d.}
$$

The only property that matters: **the data depend linearly on $\boldsymbol\beta$, and the noise is Gaussian.**

---

## 2. Reducing $M$ trajectories to their average

The log-likelihood of all $M$ trajectories is

$$
\log L(\boldsymbol\beta) = -\frac{1}{2\sigma^2\Delta t}\sum_{m=1}^{M} \big\|\mathbf y^{(m)} - B\boldsymbol\beta\big\|^2 + \text{const}.
$$

With the average $\bar{\mathbf y} = \frac1M\sum_m \mathbf y^{(m)}$, the sum splits exactly into

$$
\sum_{m=1}^{M} \big\|\mathbf y^{(m)} - B\boldsymbol\beta\big\|^2
= M\,\big\|\bar{\mathbf y} - B\boldsymbol\beta\big\|^2
+ \underbrace{\sum_{m=1}^{M} \big\|\mathbf y^{(m)} - \bar{\mathbf y}\big\|^2}_{\text{independent of } \boldsymbol\beta}.
$$

So $\bar{\mathbf y}$ is a **sufficient statistic**: nothing is lost by using only the average, which behaves like **one** data vector with $M$ times smaller noise variance:

$$
\bar{\mathbf y} = B\boldsymbol\beta + \mathbf e,
\qquad \mathbf e \sim \mathcal N\big(0,\ \sigma_n^2 I\big),
\qquad \boxed{\sigma_n^2 = \frac{\sigma^2\,\Delta t}{M}}.
$$

*(Code: `ybar{j} = mean(diff(D.' + sigma*W, 1, 2), 1).'` and `sn2 = sigma^2*dt/M`.)*

---

## 3. Gaussian prior

Zero-mean Gaussian prior $\boldsymbol\beta \sim \mathcal N(0, \Gamma)$, written through its **precision** $\Gamma^{-1} = K/\tau^2$:

| prior | penalty on $\rho$ | $K$ |
|---|---|---|
| flat | none | $0$ |
| Tikhonov ($L^2$) | $\|\rho\|_{L^2}^2 / \tau^2$ | $\Phi^\top\Phi\,\Delta x$ |
| gradient | $\|\rho'\|_{L^2}^2 / \tau^2$ | $\Phi'^\top\Phi'\,\Delta x$ |

$$
\log \pi(\boldsymbol\beta) = -\tfrac12\,\boldsymbol\beta^\top \Gamma^{-1}\boldsymbol\beta + \text{const}.
$$

---

## 4. The posterior is exactly Gaussian: completing the square

Bayes' rule, $\log\pi(\boldsymbol\beta\mid\bar{\mathbf y}) = \log L + \log\pi + \text{const}$:

$$
\log\pi(\boldsymbol\beta\mid\bar{\mathbf y})
= -\frac{1}{2\sigma_n^2}\big\|\bar{\mathbf y} - B\boldsymbol\beta\big\|^2
  - \frac12\,\boldsymbol\beta^\top\Gamma^{-1}\boldsymbol\beta + \text{const}.
$$

Expand the square, $\|\bar{\mathbf y} - B\boldsymbol\beta\|^2 = \bar{\mathbf y}^\top\bar{\mathbf y} - 2\boldsymbol\beta^\top B^\top\bar{\mathbf y} + \boldsymbol\beta^\top B^\top B\,\boldsymbol\beta$:

$$
\log\pi(\boldsymbol\beta\mid\bar{\mathbf y})
= -\frac12\,\boldsymbol\beta^\top
\underbrace{\Big(\frac{B^\top B}{\sigma_n^2} + \Gamma^{-1}\Big)}_{=:\,P}
\boldsymbol\beta
+ \boldsymbol\beta^\top \frac{B^\top\bar{\mathbf y}}{\sigma_n^2} + \text{const}.
$$

This is a **quadratic function of $\boldsymbol\beta$, exactly, for every $\boldsymbol\beta$**. Completing the square:

$$
\log\pi(\boldsymbol\beta\mid\bar{\mathbf y}) = -\frac12\,(\boldsymbol\beta - \boldsymbol\mu)^\top P\,(\boldsymbol\beta - \boldsymbol\mu) + \text{const},
$$

which is the log-density of a Gaussian:

$$
\boxed{\;
\boldsymbol\beta \mid \bar{\mathbf y} \ \sim\ \mathcal N\big(\boldsymbol\mu,\ \Sigma\big),
\qquad
\Sigma = P^{-1} = \Big(\frac{B^\top B}{\sigma_n^2} + \Gamma^{-1}\Big)^{-1},
\qquad
\boldsymbol\mu = \Sigma\,\frac{B^\top\bar{\mathbf y}}{\sigma_n^2}.
\;}
$$

*(Code: `Pp = BtB/sn2 + Pk;  mu = Pp \ (B.'*y/sn2);`, then `Sigma = inv(Pp)`.)*

**Why it is exact.** A Gaussian prior times a Gaussian likelihood that is **linear** in the unknowns gives a Gaussian posterior. This is conjugacy: no approximation is involved, unlike the Laplace approximation for $(a, w, p)$. Here the "Laplace approximation" is simply exact.

**Flat prior** ($\Gamma^{-1} = 0$): $\boldsymbol\mu = (B^\top B)^{-1}B^\top\bar{\mathbf y}$ is the **least-squares** solution and $\Sigma = \sigma_n^2 (B^\top B)^{-1}$. It is a proper distribution only if $B$ has full column rank. When $B^\top B$ is badly conditioned (condition number $10^8$–$10^9$ for $N_b = 80$–$120$), $\Sigma$ becomes huge: that is the variance explosion seen at high noise.

**With a prior**, the posterior mean is the **Tikhonov solution**:

$$
\boldsymbol\mu = \arg\min_{\boldsymbol\beta}\ \frac{1}{2\sigma_n^2}\|\bar{\mathbf y} - B\boldsymbol\beta\|^2 + \frac{1}{2\tau^2}\,\boldsymbol\beta^\top K\boldsymbol\beta .
$$

The Bayesian version adds the covariance $\Sigma$, and the regularization parameter is $\lambda = \sigma_n^2/\tau^2$.

---

## 5. Uncertainty quantification without a Markov chain

### Independent samples

Factor $\Sigma = LL^\top$ (Cholesky or eigen-decomposition). Then

$$
\boldsymbol\beta^{(n)} = \boldsymbol\mu + L\,\mathbf z^{(n)},
\qquad \mathbf z^{(n)} \sim \mathcal N(0, I)\ \text{i.i.d.},
$$

are **exact, independent** posterior samples:

- no burn-in, no proposal to tune, no acceptance step;
- no autocorrelation, so $N_{\mathrm{eff}} = N_s$ (2000 draws = 2000 effective samples);
- cost: one $N_b\times N_b$ factorization, then matrix–vector products.

*(Code: `[Vs, Ls] = eig(inv(Pp)); Lc = Vs*sqrt(Ls); betas = mu + Lc*randn(Nb, 2000);`)*

### From coefficients to the medium

Each sample gives a profile $\boldsymbol\rho^{(n)} = \Phi\boldsymbol\beta^{(n)}$ and a wave-speed ratio

$$
\frac{c^{(n)}}{c_0} = \tfrac12\Big(\rho^{(n)} + \sqrt{(\rho^{(n)})^2 + 4}\Big).
$$

The inferred medium is the pointwise mean, and the $95\%$ region is the pointwise $2.5\%$ and $97.5\%$ quantiles, exactly as for the MCMC samples of $(a, w)$.

### Even without sampling

$\boldsymbol\rho = \Phi\boldsymbol\beta$ is linear, so it is Gaussian too, **pointwise in closed form**:

$$
\rho(x_i) \mid \bar{\mathbf y} \ \sim\ \mathcal N\big((\Phi\boldsymbol\mu)_i,\ (\Phi\Sigma\Phi^\top)_{ii}\big).
$$

The map $\rho \mapsto c/c_0$ is **increasing**, so quantiles transform directly:

$$
q_{0.025}\!\left(\tfrac{c}{c_0}\right)(x_i) = \tfrac12\Big(q_\rho + \sqrt{q_\rho^2 + 4}\Big),
\qquad
q_\rho = (\Phi\boldsymbol\mu)_i - 1.96\sqrt{(\Phi\Sigma\Phi^\top)_{ii}},
$$

and similarly with $+1.96$ for the upper bound. The pointwise band therefore needs no samples at all. Samples are only needed for non-pointwise quantities, such as the posterior of the peak value $\max_x c/c_0$.

---

## 6. Choosing $\tau$ from the data (evidence), also in closed form

Integrating $\boldsymbol\beta$ out, the data alone are Gaussian:

$$
\bar{\mathbf y} \mid \tau \ \sim\ \mathcal N\big(0,\ C(\tau)\big),
\qquad C(\tau) = \sigma_n^2 I + \tau^2\, B K^{-1} B^\top .
$$

With the eigen-decomposition $BK^{-1}B^\top = Q\,\mathrm{diag}(d)\,Q^\top$, computed once, and $\mathbf z = Q^\top\bar{\mathbf y}$:

$$
\log p(\bar{\mathbf y}\mid\tau) = -\frac12\sum_i \frac{z_i^2}{\sigma_n^2 + \tau^2 d_i} - \frac12\sum_i \log\big(\sigma_n^2 + \tau^2 d_i\big) + \text{const},
$$

evaluated on a grid of $\tau$ and maximized ("empirical Bayes").

*(Code: `Gk = B / chol(K); [Q, d] = eig(Gk*Gk.', 'vector'); logev(t) = ...`)*

As seen on appendix slide A2, this choice can over-shrink at high noise. A physically chosen $\tau$ gave honest bands.

---

## 7. Where MCMC *was* used: only as a check

For $N_b = 20$ and $\sigma = 0.1$, a random-walk Metropolis–Hastings chain (60 000 iterations, proposal $\mathcal N(0, (2.38^2/N_b)\,\Sigma)$, started at $\boldsymbol\beta = 0$) reproduced the closed-form posterior:

- acceptance $24.8\%$;
- posterior means within $0.08$ standard deviations of $\boldsymbol\mu$;
- standard deviations within $3\%$ of $\sqrt{\operatorname{diag}\Sigma}$.

This confirms the formulas. But the chain needed the exact $\Sigma$ for its proposal, which is only available *because* the posterior is known in closed form. In a linear-Gaussian problem, sampling $\mathcal N(\boldsymbol\mu, \Sigma)$ directly is always better.

---

## 8. When the closed form no longer applies (and MCMC is needed again)

| situation | why it breaks | typical method |
|---|---|---|
| nonlinear parameters ($w$, $p$, several bumps) | $B\boldsymbol\beta$ becomes $D(\boldsymbol\theta)$, not linear | MH / MALA / HMC |
| unknown $\tau$ with a hyperprior | joint $(\boldsymbol\beta, \tau)$ not Gaussian | **Gibbs**: $\boldsymbol\beta \mid \tau$ Gaussian (this note), $\tau \mid \boldsymbol\beta$ closed form |
| unknown noise level $\sigma$ | marginal becomes Student-$t$ | Gibbs, or estimate $\sigma$ first (spread across trajectories) |
| non-Gaussian priors (positivity, total variation) | posterior not quadratic | MCMC, e.g. pCN or proximal samplers |
| Wasserstein misfit | not a Gaussian likelihood | MCMC on the generalized posterior |

---

## 9. Cost comparison

| | bump model $(a, w)$ with MCMC | basis model, closed form |
|---|---|---|
| work | $8000$ iterations, each one forward solve $A\rho$ ($270\times3000$) | one $N_b\times N_b$ solve and factorization |
| samples | correlated, $N_{\mathrm{eff}} \approx 300$–$400$ | independent, $N_{\mathrm{eff}} = N_s$ |
| tuning | proposal scaling, burn-in | none |

The closed form is only possible because of **linearity in the unknowns plus Gaussian noise and prior**. That is also why it does not carry over to the bump parameters $w$ and $p$.
