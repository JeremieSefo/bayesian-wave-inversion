# Bayesian Wave Inversion and Uncertainty Quantification

MATLAB code for inferring the parameters of a localized wave-speed contrast from noisy scattered-wave data with **Markov chain Monte Carlo (MCMC)**, and for **quantifying the uncertainty** of the result.

- **Report:** [`docs/report.pdf`](docs/report.pdf)
- **Slides:** [`docs/slides.pdf`](docs/slides.pdf)

The contrast is a smooth bump

$$\rho(x) = a\,\exp\!\left(-\frac{(x-x_0(p))^4}{2w^2}\right),$$

and we invert either its **amplitude and width** $(a, w)$ or its **amplitude and position** $(a, p)$. The data are $M = 100$ noisy trajectories $X_t = X_0\exp\!\big(D^{\mathrm{sc}}(t;\theta) + \sigma W_t\big)$ of the scattered wave, where $D^{\mathrm{sc}} = A\rho$ is given by a **linear forward operator** $A$ (ROM linearization).

---

## Quick start

**Requirements**

- MATLAB (tested with R2025b) with the **Statistics and Machine Learning Toolbox** (`quantile`, `chi2inv`, `skewness`).
- *Optional:* [matlab2tikz](https://github.com/matlab2tikz/matlab2tikz), install it from the MATLAB Add-On Explorer, for TikZ versions of the figures. Without it, every script still runs and saves PNG figures.

**First run**

```matlab
cd path/to/bayesian-wave-inversion
setup                       % adds src/, data/ and experiments/ to the path (once per session)
make_data                   % checks the data and plots the true medium and the data
WAVE_MCMC                   % first inversion: amplitude and width
```

Figures are written to `figures/` and results to `results/`. Both folders are created automatically and are not tracked by git.

---

## Repository layout

```
bayesian-wave-inversion/
├── setup.m                    run first: sets the MATLAB path
├── data/
│   ├── waveData.mat           forward operator and true medium (see "The data")
│   └── make_data.m            regenerates / checks the data, or builds a new setting
├── src/                       functions used by the experiments
├── experiments/               scripts, numbered in learning order (01 ... 08)
├── notes/                     derivations (Markdown with LaTeX math)
├── docs/                      report.pdf, slides.pdf
├── results/                   (generated) .mat results reused by later scripts
└── figures/                   (generated) report/ and slides/ figures
```

---

## The data

`data/waveData.mat` contains **this project's specific setting**:

| variable | size | meaning |
|---|---|---|
| `TrueOper` | 270 × 3000 | linear forward operator $A$: contrast on the grid → scattered data at 270 times $t \in [0, 1.5]$ |
| `Contrast` | 3000 × 1 | true contrast $\rho(x)$ on $x \in [0, 2]$ ($a = 0.75$, $w = 5\cdot10^{-3}$, $p = 0.14$) |
| `C0`, `C1` | 3000 × 1 | background and true wave speed ($c/c_0$ = `C1./C0`) |
| `TrueScattData_noisy` | 270 × 1 | **noise-free** scattered data $D = A\rho$ (historical name: the noise is added in the scripts) |

**What is primary and what can be regenerated.** Only the operator `TrueOper` is primary data. It was computed with a wave solver outside this repository, using the known internal wave field of the reduced-order-model framework (Borcea et al., *SIAM Review* 66(3), 2024). Everything else is regenerated exactly by `data/make_data.m`:

- the true medium, from the parameters with `MediumGenerator.medium3` (identical to the stored values up to rounding);
- the data $D = A\rho$ (relative difference $3\cdot10^{-6}$ from the stored values).

**Testing other settings.** The experiments use this specific data, but can be modified to test other settings:

- **Another true medium** (other amplitude, width or position): edit the parameters at the top of `make_data.m`, set `write = true`, and run it. The experiments then use the new `waveData.mat`. The true values used for comparison are set at the top of each experiment (`amp_true`, `width_true`, `pos_true`).
- **Another operator** (other sources/receivers, frequencies, background): replace `TrueOper` with your own $N_t \times n$ matrix and regenerate the data with `make_data.m`. The scripts assume the spatial grid (3000 points on $[0, 2]$) and record length $T = 1.5$ used here; change `n`, `Le` and `T` in the scripts if your operator uses another grid.
- **Another noise level or number of trajectories:** `sigNoise` / `noise_levels` and `M` at the top of the scripts.

---

## Learning path

Run the scripts in this order. Each step builds on the previous one. Run times are for a laptop (R2025b).

### 0. Setup and data
| file | what it does | time |
|---|---|---|
| `setup.m` | adds the repository folders to the MATLAB path | < 1 s |
| `data/make_data.m` | regenerates the true medium and data from the operator, compares them with `waveData.mat`, plots them | 2 s |

### 1. First inversion: amplitude and width — `experiments/01_amplitude_width/`
| file | what it does | time |
|---|---|---|
| `WAVE_MCMC.m` | simulates the noisy trajectories, computes the MLE, runs Metropolis–Hastings for $(a, w)$, plots chains, samples and Monte Carlo error decay | 12 s |
| `see_posteriorWave.m` | marginal and joint posterior histograms (run after `WAVE_MCMC`, uses its workspace) | 2 s |
| `ampWidth_uncerntainty.m` | uncertainty region of the medium: every sample → $c/c_0$, pointwise mean and 2.5–97.5% quantiles (after `WAVE_MCMC`) | 1 s |
| `UQ_vs_noise.m` | repeats the inversion for 6 noise levels: posterior mean, 95% credible intervals, std, covariance (after `WAVE_MCMC`) | 17 s |

Uses `src/`: `MediumGenerator.m` (contrast profile), `wave_loglik.m` (log-likelihood of the log-increments), `run_mcmc_wave.m` (sampler).

### 2. Amplitude and position — `experiments/02_amplitude_position/`
| file | what it does | time |
|---|---|---|
| `WAVEpos_MCMC.m` | same as step 1 for $(a, p)$; ends by calling `amp_uncerntainty.m` and `UQ_vs_noisePos.m` | 45 s |
| `amp_uncerntainty.m`, `UQ_vs_noisePos.m` | uncertainty region and noise study for $(a, p)$ | (included above) |

Uses `src/wave_pos_loglik.m` and `src/run_mcmc_wavePos.m`.

### 3. Covariance, correlation and the Laplace approximation — `experiments/03_covariance_correlation/`
| file | what it does | time |
|---|---|---|
| `run_both_samplers.m` | runs `compare_cov_noise.m` twice: original fixed-step sampler and proposal scaled to the posterior | 40 s |
| `compare_cov_noise.m` | both experiments on the same noisy data: posterior covariance and correlation vs noise, compared with the Laplace prediction $\Sigma \approx \sigma^2 G^{-1}$ | (called above) |
| `plot_cov_corr_figures.m` | figures: $\mathrm{Cov}/\sigma^2$ and correlation vs noise | 11 s |

Background: [`notes/covariance_derivation.md`](notes/covariance_derivation.md). Covariance grows like $\sigma^2$, and the correlation does not depend on $\sigma$.

### 4. Why the proposal scaling matters — `experiments/04_proposal_scaling/`
| file | what it does | time |
|---|---|---|
| `plot_chain_comparison.m` | original vs scaled sampler: traces, samples with the Laplace 95% ellipse, autocorrelation, effective sample size (ESS) table | 9 s |
| `make_slide_proposal.m` | the same comparison with accuracy measures (std and correlation vs Laplace) | 2 s |
| `mh_dimension_test.m` | random-walk MH in 2, 20 and 80 dimensions: diagonal vs full-covariance proposal, ESS per iteration | 50 s |

### 5. Two-phase sampler: inversion and uncertainty quantification — `experiments/05_two_phase_sampler/`
| file | what it does | time |
|---|---|---|
| `run_slide_results.m` | pilot phase from a far initial guess (fixed steps), then a production phase with steps from $G$ at the best pilot state; estimates, credible intervals, uncertainty regions for all noise levels | 33 s |
| `make_slide_ciwidth.m` | credible-interval width vs noise on log-log axes, with fitted slopes (≈ 1: width ∝ $\sigma$) | 2 s |

### 6. Shape of the posterior and the Wasserstein misfit — `experiments/06_posterior_shape_wasserstein/`
| file | what it does | time |
|---|---|---|
| `gaussianity_check.m` | is the posterior Gaussian away from the truth? Likelihood profiles along $w$ and $p$, local maxima, skewness and kurtosis of the samples | < 1 s |
| `make_slide_landscape.m` | $L^2$ misfit (local minima from cycle skipping) vs Wasserstein $W_2$ misfit (single well) along the position | 2 s |

### 7. Beyond two parameters — `experiments/07_beyond_two_parameters/`
| file | what it does | time |
|---|---|---|
| `make_slide_outlook.m` | SVD of the operator: how many modes of the contrast rise above the noise (Picard plot) | 2 s |
| `hat_basis_inversion.m` | contrast in a hat-function basis: closed-form Gaussian posterior, checked with MCMC | 3 s |
| `basis_compare.m` | hat vs Gaussian basis: model error, posterior error and coverage | 2 s |
| `regularized_inversion.m` | Bayesian regularization ($L^2$ and gradient priors) on top of the basis; prior strength from the evidence or fixed physically | 10 s |
| `make_appendix_regularized.m` | figure: flat prior vs evidence-chosen vs physical prior | 2 s |

Background: [`notes/basis_closed_form_posterior.md`](notes/basis_closed_form_posterior.md), on why the basis case needs no MCMC.

### 8. Slide figures — `experiments/08_slide_figures/`
| file | what it does | time |
|---|---|---|
| `make_slide_figures.m` | model, noise and variance-vs-noise figures of the slides (needs step 3) | 4 s |
| `make_all_tikz.m` | re-runs the slide-figure scripts and exports every panel to TikZ in `figures/slides/tikz` (needs matlab2tikz) | ~3 min |

---

## Which results feed which scripts

Some scripts reuse saved results, so run their producer first:

| result file (`results/`) | produced by | used by |
|---|---|---|
| `compare_cov_noise_results.mat`, `..._tuned.mat` | `run_both_samplers.m` (step 3) | `plot_cov_corr_figures`, `plot_chain_comparison`, `make_slide_proposal`, `make_slide_ciwidth`, `make_slide_figures` |
| `slide_results.mat` | `run_slide_results.m` (step 5) | `make_slide_ciwidth`, `make_slide_proposal`, `gaussianity_check` |
| `regularized_results.mat` | `regularized_inversion.m` (step 7) | `make_appendix_regularized` |

Running all steps in order takes about **5 minutes** (plus about 3 minutes for the TikZ export).

## Helper functions (`src/`)

| file | purpose |
|---|---|
| `repo_paths.m` | locations of `data/`, `results/`, `figures/`, independent of where the repository is cloned |
| `MediumGenerator.m` | builds the contrast profile and wave-speed model on the grid |
| `wave_loglik.m`, `wave_pos_loglik.m` | log-likelihood of the log-increments for $(a, w)$ and $(a, p)$ |
| `run_mcmc_wave.m`, `run_mcmc_wavePos.m` | fixed-step random-walk Metropolis–Hastings samplers |
| `rel_noise_factor.m`, `noise_axis.m`, `pct_str.m` | convert $\sigma$ to the effective relative noise (≈ 3.4σ with M = 100) and label the axes |
| `saveTikz.m`, `export_tikz_panels.m` | save figures as PNG and (if available) TikZ |

## Reproducibility

All random numbers are seeded (the noisy data use `rng(100 + j)` for noise level `j`), so the scripts reproduce the numbers in the report and slides exactly.
