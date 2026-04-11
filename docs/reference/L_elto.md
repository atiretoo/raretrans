# Transition and recruitment data for *Lepanthes eltoroensis*

Individual-level census data from multiple populations of the epiphytic
orchid *Lepanthes eltoroensis* Stimson, a critically rare species
endemic to Puerto Rico. Individuals were permanently marked and surveyed
at 6-month intervals over 12 census periods (6 calendar years). Each row
represents one individual observed at one census period.

## Usage

``` r
data(L_elto)
```

## Format

A data frame with 5233 rows and 13 variables:

- POPNUM:

  Integer. Population identifier number.

- year:

  Integer. Census period number (each period represents a 6-month
  interval, not a calendar year; periods range from 1 to 12).

- seedlings:

  Numeric. Population-level count of plantulas (seedlings, stage `"p"`)
  observed in the population at this census period. The same value is
  repeated for all individuals in the same population-period
  combination.

- adults:

  Numeric. Population-level count of adults (stage `"a"`) observed in
  the population at this census period. The same value is repeated for
  all individuals in the same population-period combination.

- fertility:

  Numeric. Mean per-individual fertility: average number of seedlings
  produced per adult in the population at this census period.

- IND_NUM:

  Integer. Unique individual plant identifier number within its
  population.

- stage:

  Character. Life history stage of the individual at the current census:
  `"p"` (plantula/seedling), `"j"` (juvenile), or `"a"` (adult). See
  Details for stage definitions.

- next_stage:

  Character. Life history stage of the individual at the following
  6-month census: `"p"`, `"j"`, `"a"`, or `"m"` (muerto/dead). See
  Details for stage definitions.

- first_year:

  Integer. 6-month census period number when the individual was first
  observed (1 = first survey).

- last_year:

  Integer. 6-month census period number when the individual was last
  observed alive.

- recruited:

  Logical. `TRUE` if this individual was newly recruited into the
  population (not present at the previous survey), typically via sexual
  reproduction. All individuals at the first census period (`year == 1`)
  are `TRUE` by definition as there is no prior survey for comparison.

- died:

  Logical. Whether the individual was found dead at this census period.

- lifespan:

  Integer. Total number of 6-month census periods during which the
  individual was observed alive (not in calendar years).

## Source

Raymond L. Tremblay, University of Puerto Rico at Humacao, unpublished
data.

## Details

Life history stages are coded as follows:

- `"p"` - plantula (seedling): individuals without a lepanthiform sheet
  on any of the leaves

- `"j"` - juvenile: individuals with no evidence of present or past
  reproductive effort; bases of inflorescences are persistent

- `"a"` - adult: individuals with active or inactive inflorescences
  present

- `"m"` - muerto (dead): individual was dead at the following census
  (used only in `next_stage`)

## References

Tremblay, R.L. and Hutchings, M.J. (2003). Population dynamics in orchid
conservation: a review of analytical methods based on the rare species
*Lepanthes eltoroensis*. In: Dixon, K.W. et al. (eds.) *Orchid
Conservation*. Natural History Publications (Borneo), pp. 183–204.

Tremblay, R.L., Perez, M-E., Tyre, A.J. and Tenhumberg, B. (2021).
Bayesian estimates for matrix population models of rare plants.
*Ecological Modelling*, 440, 109526.
[doi:10.1016/j.ecolmodel.2021.109526](https://doi.org/10.1016/j.ecolmodel.2021.109526)
