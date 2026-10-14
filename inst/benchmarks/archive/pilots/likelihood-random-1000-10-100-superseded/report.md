# Fixed-parameter likelihood validation

6000 of 6000 paired checks passed the prespecified absolute tolerance of 1e-5.

Maximum finite absolute difference: 2.1967244e-06.

Results by model and survival conditioning:

```
       comparison survival    n finite passed max.absolute.error
           GeoSSE    FALSE 1000   1000   1000       2.052343e-06
           GeoSSE     TRUE 1000   1000   1000       2.196724e-06
 GeoSSE-as-ClaSSE    FALSE 1000   1000   1000       2.052343e-06
 GeoSSE-as-ClaSSE     TRUE 1000   1000   1000       2.196724e-06
      Full-ClaSSE    FALSE 1000   1000   1000       4.271334e-08
      Full-ClaSSE     TRUE 1000   1000   1000       4.099141e-08
 median.absolute.error
          8.272423e-10
          8.198811e-10
          8.272423e-10
          8.198811e-10
          6.079333e-10
          6.059508e-10
```

These checks establish numerical equivalence only over the tested cases. They do not validate genuinely changing epoch rates, parameter recovery, or stochastic mapping.
Both epochs use identical parameters. Root treatment is ROOT.OBS, with complete sampling. The main plots show survival-conditioned results; the tables retain both settings.
