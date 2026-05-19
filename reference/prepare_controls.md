# Prepare centred controls and D\*X interactions

Prepare centred controls and D\*X interactions

## Usage

``` r
prepare_controls(df_sample, d, controls)
```

## Arguments

- df_sample:

  Cross-section data frame

- d:

  Treatment indicator column name

- controls:

  Character vector of control variable names (or NULL)

## Value

List with: include, X_centered, interactions, RHS_varnames
