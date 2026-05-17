#set page(
  paper: "us-letter",
  margin: (x: 1in, y: 1in),
  numbering: "1",
)
#set par(justify: true, leading: 0.65em)
#set text(font: "New Computer Modern", size: 11pt)
#set heading(numbering: "1.")
#set math.equation(numbering: "(1)")

#align(center)[
  #text(size: 18pt, weight: "bold")[
    Hamiltonian Equipartition as a Convergence Diagnostic for Markov Chains
  ]

  #v(0.5em)

  #text(size: 12pt)[A Theoretical Report on the Momentum Checker]

  #v(0.5em)

  #text(size: 10pt, style: "italic")[
    May 2026
  ]
]

#v(1em)

#align(center)[
  #block(width: 90%, [
    *Abstract.*  We describe and justify a convergence diagnostic for Markov
    chain Monte Carlo (MCMC) samplers that exploits the factorization of the
    canonical distribution in phase space.  Given a Markov chain targeting the
    Boltzmann distribution $pi(x) prop e^(-beta V(x))$, we draw positions from
    the chain, attach freshly sampled momenta from the canonical momentum
    distribution, evolve the joint state under Hamiltonian dynamics, and
    examine the empirical variance of the resulting momenta.  By the
    equipartition theorem, this variance must equal $m \/ beta$ at
    equilibrium.  Any persistent deviation is direct evidence that the
    position chain has not converged.  The diagnostic is global, uses a known
    exact reference distribution, and naturally admits a one-parameter family
    of probes indexed by particle mass.
  ])
]

#v(1em)

= Introduction

Markov chain Monte Carlo methods underpin much of modern Bayesian inference
and computational statistical mechanics.  Their utility is asymptotic: in the
limit of infinite samples, ergodic chains produce draws from the desired
target distribution.  In any finite simulation, however, the sample path
carries a bias inherited from its starting state and from the slow modes of
the transition kernel.  Detecting when this bias has decayed to negligible
levels — that is, deciding when a chain has "converged" — is a central
practical problem.

Standard diagnostics fall into two broad families.  The first compares
trajectories of multiple independent chains, the prototype being the
Gelman–Rubin $hat(R)$ statistic [1], which contrasts
between-chain and within-chain variance.  The second analyzes a single chain
through its autocorrelation structure, summarized in the effective sample
size [2].  Both families are agnostic to the form of the
target: they treat the chain as a generic stationary process and ask whether
its empirical moments have stabilized.

The diagnostic developed here is different in kind.  It exploits a structural
property of distributions of the form $pi(x) prop e^(-beta V(x))$: when such
a distribution is augmented with a Gaussian momentum and evolved under
Hamiltonian dynamics, the marginal momentum distribution is preserved
exactly.  Deviations from the expected momentum variance are therefore
unambiguous evidence that the position distribution differs from $pi$.  This
gives a diagnostic that is global (it senses any disagreement, not just in
specific test functions chosen a priori), parametric (the mass introduces a
tunable time scale), and grounded in a known closed-form reference.

The rest of the report develops this in detail.  Section 2 sets up the
canonical distribution and its phase-space factorization.  Section 3 reviews
the relevant facts about Hamiltonian dynamics.  Section 4 specifies the
algorithm.  Section 5 contains the central theoretical claim and its proof.
Section 6 discusses the mass scan as a multi-scale probe.  Sections 7 and 8
compare the diagnostic to standard tools and discuss its limitations.

= The Canonical Distribution and Its Factorization

Let $V: bb(R)^d -> bb(R)$ be a potential energy function and let
$beta > 0$ be an inverse temperature.  The MCMC target is
$
  pi(x) = e^(-beta V(x)) / Z, quad Z = integral_(bb(R)^d) e^(-beta V(x)) dif x,
$
where $Z$ is the normalizing constant.  We work in the augmented phase space
$bb(R)^d times bb(R)^d$, introducing an auxiliary momentum variable
$p in bb(R)^d$ and a mass parameter $m > 0$.  Define the Hamiltonian
$
  H(x, p) = V(x) + (||p||^2)/(2 m).
$
The associated canonical (Gibbs) distribution on phase space is
$
  rho(x, p) prop e^(-beta H(x, p)) = e^(-beta V(x)) dot e^(-beta ||p||^2 / (2 m)).
$ <eq:canonical>

The decisive observation is that the right-hand side of @eq:canonical
factorizes.  Writing
$
  rho(x, p) = pi(x) dot mu(p),
$
we identify
$
  mu(p) prop e^(-beta ||p||^2 / (2 m)) = product_(i = 1)^d e^(-beta p_i^2 / (2 m))
$
as a product of independent Gaussians with mean zero and variance $m \/ beta$
in each component.  In particular, in the canonical ensemble *position and
momentum are statistically independent*, and the marginal momentum
distribution is the same regardless of the potential $V$.

This independence is at the heart of the diagnostic.  The momentum has a
universal, parameter-free marginal distribution at equilibrium; any chain
purporting to sample $pi$ can be tested by checking whether the implied joint
distribution, when propagated forward in time by Hamiltonian dynamics,
respects the universal momentum law.

#pagebreak()

= Hamiltonian Dynamics and Preservation of the Canonical Measure

Given the Hamiltonian $H(x, p)$, the equations of motion are
$
  dot(x) = (partial H)/(partial p) = p / m, quad
  dot(p) = -(partial H)/(partial x) = -nabla V(x).
$ <eq:hamilton>

Let $Phi_t: (x_0, p_0) |-> (x_t, p_t)$ denote the flow generated by
@eq:hamilton over a time interval $t$.  The flow has two properties of
central importance.

#par(first-line-indent: 0pt)[
  *Energy conservation.*  Along any trajectory,
  $
    (dif)/(dif t) H(x(t), p(t)) = nabla V dot dot(x) + p/m dot dot(p)
    = nabla V dot p/m - p/m dot nabla V = 0.
  $
  Thus $H compose Phi_t = H$ for all $t$.
]

#par(first-line-indent: 0pt)[
  *Volume preservation (Liouville's theorem).*  The flow is symplectic, hence
  preserves the $2d$-dimensional Lebesgue measure on phase space:
  $
    integral_A dif x dif p = integral_(Phi_t (A)) dif x dif p
    quad "for every measurable" A subset bb(R)^(2d).
  $
]

Combining these two facts, the canonical density $rho(x, p) prop e^(-beta H(x, p))$
is invariant under $Phi_t$.  Concretely, if $(X, P) tilde rho$ and
$(X', P') = Phi_t(X, P)$, then $(X', P') tilde rho$ as well.  This is the
classical statement that Hamiltonian dynamics preserves the canonical
ensemble.

Two corollaries are immediate.  First, since $rho = pi dot mu$ is the
product of two independent marginals, propagation by $Phi_t$ maps the joint
law $pi times mu$ to itself.  Second, the marginal of $rho$ over the position
coordinate is $pi$, so the position-conditional momentum distribution at any
$t$ is again $mu$.  These will be reused below.

= The Momentum Checker Algorithm

Let $X_1, dots, X_N$ be position samples produced by a Markov chain targeting
$pi$, and let $hat(pi)$ denote the empirical distribution of these samples.
The algorithm constructs a secondary Markov chain $P_0, P_1, dots, P_K$ on
momentum space as follows.

+ *Initialization.*  Draw $P_0 tilde mu = cal(N)(0, (m\/beta) I_d)$.
+ *Iteration.*  For $k = 0, 1, dots, K - 1$:
  + Draw an index $j_k$ uniformly from ${1, dots, N}$.
  + Set $X^((k)) <- X_(j_k)$ (resample position from $hat(pi)$).
  + Compute $(tilde(X), tilde(P)) = Phi_(t_k)(X^((k)), P_k)$ using a
    numerical integrator (here, leapfrog over a Poisson-sampled number of
    steps).
  + Discard $tilde(X)$ and set $P_(k + 1) <- tilde(P)$.

+ *Statistic.*  Compute the empirical variance
  $
    hat(sigma)^2 = 1/(K - 1) sum_(k = 1)^K (P_k - macron(P))^2,
    quad macron(P) = 1/K sum_(k = 1)^K P_k.
  $
  Compare $hat(sigma)^2$ to the theoretical value $m \/ beta$.

In words: at each step we replace the position by a fresh draw from
$hat(pi)$, then let the joint state $(X, P)$ evolve under the Hamiltonian
flow, then forget the position.  The momentum alone is carried forward.  The
diagnostic asks whether the long-run distribution of $P_k$ matches the
canonical momentum law $mu$.

= Why It Detects Non-Convergence

#par(first-line-indent: 0pt)[
  *Claim.*  Suppose $hat(pi) = pi$ exactly.  Then $P_k tilde mu$ for every
  $k$, and consequently $hat(sigma)^2 -> m \/ beta$ almost surely as
  $K -> infinity$.
]

#par(first-line-indent: 0pt)[
  *Proof.*  By induction on $k$.  At $k = 0$ the claim holds by construction.
  Suppose $P_k tilde mu$.  In step (ii), $X^((k))$ is an independent draw
  from $hat(pi) = pi$.  Therefore $(X^((k)), P_k)$ has joint law
  $pi times mu = rho$, the canonical distribution on phase space.  By the
  invariance result of Section 3, $(tilde(X), tilde(P)) = Phi_(t_k)(X^((k)), P_k)$
  also has law $rho$.  The marginal of $rho$ in the momentum coordinate is
  $mu$, so $P_(k + 1) = tilde(P) tilde mu$.  The strong law of large numbers
  then gives $hat(sigma)^2 -> "Var"(mu) = m \/ beta$. $square$
]

This is the positive half of the diagnostic.  The negative half — that
$hat(sigma)^2 != m \/ beta$ when $hat(pi) != pi$ — is the operational claim
on which the diagnostic depends, and we treat it next.

#par(first-line-indent: 0pt)[
  *Why disagreement shows up.*  Suppose now that $hat(pi)$ is some
  distribution $tilde(pi) != pi$.  At step $k$, the joint law of
  $(X^((k)), P_k)$ is $tilde(pi) times mu$, which is *not* invariant under
  $Phi_t$ unless $tilde(pi) prop e^(-beta V)$.  Energy conservation along a
  trajectory gives
  $
    V(tilde(X)) + (||tilde(P)||^2)/(2 m) = V(X^((k))) + (||P_k||^2)/(2 m),
  $
  so the post-evolution kinetic energy is
  $
    (||tilde(P)||^2)/(2 m) = (||P_k||^2)/(2 m) + V(X^((k))) - V(tilde(X)).
  $ <eq:kinetic-update>
  Taking expectations under the joint law $tilde(pi)(X^((k))) times mu(P_k)$,
$
    bb(E)((||tilde(P)||^2)/(2 m)) = bb(E)((||P_k||^2)/(2 m))
      + bb(E)_(tilde(pi))(V(X^((k)))) - bb(E)(V(tilde(X))).
  $
  Under the *true* equilibrium $tilde(pi) = pi$, the two potential terms
  cancel because the law of $tilde(X)$ marginalizes to $pi$.  Under any other
  $tilde(pi)$, the cancellation fails, and the kinetic energy systematically
  drifts.  Iterating, the stationary kinetic energy of the $P_k$ chain
  deviates from its canonical value of $d \/ (2 beta)$.  Equivalently, the
  per-coordinate variance of $P_k$ deviates from $m \/ beta$.
]

The diagnostic is therefore not detecting some subtle distributional feature:
it is detecting a violation of equipartition.  This makes it both
interpretable and sharp.

#pagebreak()

= The Mass Scan as a Multi-Scale Probe

The diagnostic depends on a free parameter, the mass $m$.  Two different
choices probe the same chain on different physical time scales, and a scan
over $m$ provides a richer view than any single value.

The natural velocity scale of a particle at temperature $beta^(-1)$ and mass
$m$ is
$
  v_("rms") = sqrt(chevron.l ||p / m||^2 chevron.r) = sqrt(d / (m beta)),
$
so light particles move quickly and heavy particles move slowly.  Over a
fixed simulation time $tau$, the lightest particle in a scan covers a
distance proportional to $tau \/ sqrt(m)$.

This has two important consequences for the diagnostic:

+ *Light masses probe global structure.*  A small $m$ means the particle
  sweeps through a large region of position space during $tau$.  The momentum
  at the end of the trajectory is sensitive to the full mass of position
  samples encountered along the way, including regions $hat(pi)$ may have
  failed to populate correctly.

+ *Heavy masses probe local structure.*  Large $m$ means the particle barely
  moves during $tau$.  The diagnostic then sees only the local potential
  curvature near $X^((k))$, and is less sensitive to long-range distributional
  errors.

A scan plot of $hat(sigma)(m)$ against $m$ thus has interpretive value: the
expected curve is $sqrt(m \/ beta)$, growing as $sqrt(m)$, and systematic
deviations at small $m$ versus large $m$ have different physical meanings.

A second, more practical concern is that small $m$ couples to numerical
stability of the leapfrog integrator.  The position update is
$x <- x + (epsilon \/ m) p$, so for fixed $epsilon$ the per-step displacement
grows as $m^(-1)$.  In the diagnostic code this is mitigated by scaling
$epsilon$ linearly with $m$, keeping the ratio $epsilon \/ m$ — and hence
the numerical stability of the leapfrog — constant across the mass sweep.

= Relation to Standard Diagnostics

It is worth situating the momentum checker among other tools.

#par(first-line-indent: 0pt)[
  *Gelman–Rubin $hat(R)$.*  This contrasts variances within and between
  multiple chains.  It requires running parallel chains from overdispersed
  starting points, and detects failures of mixing between modes that no
  single chain has crossed.  The momentum checker uses a single chain (the
  position chain may be one HMC, one Metropolis–Hastings, or a concatenation
  of several) and detects failures of the *marginal* distribution to equal
  $pi$.  The two diagnostics are not redundant: $hat(R)$ can pass while the
  momentum check fails (the chains agree on a wrong distribution) and vice
  versa (the chains disagree but each marginally matches $pi$ —
  unusual but possible).
]

#par(first-line-indent: 0pt)[
  *Effective sample size and integrated autocorrelation.*  These estimate
  the variance reduction of empirical averages over independent samples.
  They are silent about whether the underlying mean is correct; a chain
  stuck in one mode of a multimodal target may report high ESS and still
  be wrong.  The momentum checker, in contrast, fails systematically if the
  mode weights are wrong, because the implied potential-energy expectation
  used in @eq:kinetic-update is incorrect.
]

#par(first-line-indent: 0pt)[
  *Test-function diagnostics.*  One can choose specific observables
  $f_1, dots, f_J$ and check whether $bb(E)_(hat(pi))(f_j)$ matches some
  known value.  The momentum checker can be seen as a specific instance of
  this strategy, with $f$ implicitly defined by the Hamiltonian flow.  The
  advantage is that no test function need be selected a priori: the dynamics
  itself integrates the position information into a single scalar (the
  momentum variance) whose target value is known exactly from physics.
]

= Limitations

#par(first-line-indent: 0pt)[
  *Discretization bias.*  The Hamiltonian flow $Phi_t$ is approximated by
  the leapfrog integrator.  Leapfrog is symplectic but not exact: it
  conserves a *modified* Hamiltonian $tilde(H) = H + O(epsilon^2)$, not $H$
  itself.  The invariance argument therefore applies to the modified
  canonical density $exp(-beta tilde(H))$, not $exp(-beta H)$.  The
  difference is $O(epsilon^2)$ for leapfrog and decays as $epsilon -> 0$;
  in practice the diagnostic must use $epsilon$ small enough that this bias
  is below the precision of interest.
]

#par(first-line-indent: 0pt)[
  *Finite-sample noise.*  The empirical variance $hat(sigma)^2$ has its own
  sampling distribution.  For a Gaussian sample of size $K$, the standard
  error of the variance estimator scales as $sqrt(2/K) dot sigma^2$, so
  $K = O(10^3)$ suffices for percent-level precision in benign cases but
  may be inadequate for detecting subtle bias.  Burn-in must also be
  discarded before computing $hat(sigma)^2$, and the choice of cutoff
  affects sensitivity.
]

#par(first-line-indent: 0pt)[
  *Multimodal targets.*  If $hat(pi)$ visits only one of several modes, the
  potential-energy expectation $bb(E)_(hat(pi))(V)$ may still happen to
  match its true value if the mode is symmetric or representative.  The
  diagnostic then misses the failure.  Combining the momentum check with
  multi-chain comparisons (Gelman–Rubin) closes this gap.
]

#par(first-line-indent: 0pt)[
  *Cost.*  Each diagnostic iteration requires a leapfrog trajectory of
  length comparable to a full HMC step, so the diagnostic is roughly as
  expensive as generating one HMC sample per check iteration.  A scan over
  $J$ masses with $K$ iterations each multiplies this cost by $J K$.  For
  typical settings ($J approx 30$, $K approx 10^3$) the overhead is
  comparable to running the original chain itself.
]

#pagebreak()

= Conclusion

The momentum checker is a diagnostic with an unusually transparent
theoretical basis: it tests the equipartition theorem, a closed-form
prediction of statistical mechanics, against a Hamiltonian propagation of
the chain's empirical position distribution.  Convergence of the chain
implies the test passes; failure of the test implies the chain has not
converged.

What makes the diagnostic particularly informative in the MCMC setting is
that the reference distribution against which the chain is compared is
*universal* — the same Gaussian regardless of $V$ — and *exact* — known in
closed form independent of any numerics.  Both properties are rare in
convergence diagnostics, which typically rely either on comparisons among
independent chains (whose mutual agreement is not the same as agreement with
$pi$) or on stationarity of summary statistics over time (which says nothing
about the value to which they have stationarized).  By contrast, the
expected value $m \/ beta$ has a derivation a few lines long, and any
sustained deviation is direct evidence of bias.

The mass sweep adds a second axis of diagnosis: rather than producing a
single number, it produces a curve, and the shape of that curve carries
interpretive information about whether the chain has failed at small or
large scales.  Together with conventional tools — $hat(R)$, autocorrelation,
trace plots — the momentum checker provides a physically motivated,
distribution-aware complement that is especially natural in settings where
HMC is already being used and the Hamiltonian machinery is already
implemented.

#v(2em)

#text(weight: "bold")[References]

#text(size: 9pt)[
  - [1] A. Gelman and D. B. Rubin.  _Inference from iterative simulation
    using multiple sequences._  Statistical Science, 7(4):457–472, 1992.
  - [2] C. J. Geyer.  _Practical Markov chain Monte Carlo._  Statistical
    Science, 7(4):473–483, 1992.
  - [3] R. M. Neal.  _MCMC using Hamiltonian dynamics._  In _Handbook of
    Markov Chain Monte Carlo_, pages 113–162. Chapman & Hall/CRC, 2011.
  - [4] M. Betancourt.  _A conceptual introduction to Hamiltonian Monte
    Carlo._  arXiv:1701.02434, 2017.
  - [5] L. D. Landau and E. M. Lifshitz.  _Statistical Physics, Part 1_,
    3rd edition.  Pergamon Press, 1980.
]
