# Package index

## Package overview

- [`raretrans`](https://atiretoo.github.io/raretrans/reference/raretrans-package.md)
  [`raretrans-package`](https://atiretoo.github.io/raretrans/reference/raretrans-package.md)
  : raretrans: Bayesian Priors for Matrix Population Models

## Main functions

Core functions for filling structural zeros

- [`fill_transitions()`](https://atiretoo.github.io/raretrans/reference/fill_transitions.md)
  : Combine prior and data for transition matrix
- [`fill_fertility()`](https://atiretoo.github.io/raretrans/reference/fill_fertility.md)
  : Combine prior and data for fertility matrix
- [`run_app()`](https://atiretoo.github.io/raretrans/reference/run_app.md)
  : Run the interactive Shiny application for Bayesian transition priors

## Simulation

Simulate posterior matrices

- [`sim_transitions()`](https://atiretoo.github.io/raretrans/reference/sim_transitions.md)
  : Simulate population projection matrices

## Credible intervals

Compute and plot credible intervals

- [`transition_CrI()`](https://atiretoo.github.io/raretrans/reference/transition_CrI.md)
  : Calculate beta credible intervals for all transition matrix entries
- [`plot_transition_CrI()`](https://atiretoo.github.io/raretrans/reference/plot_transition_CrI.md)
  : Plot beta credible intervals for transition matrix entries
- [`plot_transition_density()`](https://atiretoo.github.io/raretrans/reference/plot_transition_density.md)
  : Plot posterior beta density curves for all transition matrix entries

## Helper functions

- [`get_state_vector()`](https://atiretoo.github.io/raretrans/reference/get_state_vector.md)
  : Helper functions for generating different priors
- [`rdirichlet()`](https://atiretoo.github.io/raretrans/reference/rdirichlet.md)
  : Dirichlet distributed random numbers

## Data

- [`L_elto`](https://atiretoo.github.io/raretrans/reference/L_elto.md) :

  Transition and recruitment data for *Lepanthes eltoroensis*
