# Fixed-parameter likelihood validation

1,000 fresh trees with 100--500 tips; traditional GeoSSE and unrestricted ClaSSE references only.

3996 of 4000 paired checks passed the prespecified absolute tolerance of 1e-5.

Maximum finite absolute difference: 4.0714792e-07.

Results by model and survival conditioning:

```
 comparison survival    n finite passed max.absolute.error
     GeoSSE    FALSE 1000    998    998       4.071479e-07
     GeoSSE     TRUE 1000    998    998       2.838553e-07
     ClaSSE    FALSE 1000   1000   1000       5.759853e-08
     ClaSSE     TRUE 1000   1000   1000       5.203583e-08
 median.absolute.error
          9.273435e-10
          9.376890e-10
          7.584049e-10
          7.553353e-10
```

Nonfinite or failed paired evaluations: 4. See results/likelihood-checks.csv for error messages; these are not counted as agreement.

These checks establish numerical equivalence only over the tested cases. They do not validate genuinely changing epoch rates, parameter recovery, or stochastic mapping.
Both epochs use identical parameters. Root treatment is ROOT.OBS, with complete sampling. The main plots show survival-conditioned results; the tables retain both settings.
