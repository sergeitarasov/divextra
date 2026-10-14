# Speciation-change recovery: 50 trees per represented boundary-size group

## Scientific design

Old epoch: s=.1,d=.1,x=.05. Recent epoch: s=.3,d=.1,x=.05. Four estimated rates: recent s, old s, shared d and shared x. Boundary is known, fixed to its generating age (5%,10%,25% of total simulated stem duration before present). Start uniformly among A,B,AB from one stem; complete sampling. Crown reconstructed tree supplied to likelihood.

650 distinct trees: 50 per represented group. All five size groups (50–99,100–249,250–499,500–999,1000–2000 tips) at5% and10%; only250–499,500–999,1000–2000 at25%. The two empty small-tree25% groups are explicitly omitted, not filled using duplicates or relaxed filters. All retained histories reach full duration, have50–2000 extant tips and at least30 reconstructed boundary-crossing lineages. No additional state-count filter.

Reuse independently simulated accepted trees from ../speciation-change-pilot-05-10-30 (only5%/10%) and ../speciation-change-pilot-25. Fill shortages using fresh seeds at targeted durations:5%/500–999 at50,5%/1000–2000 at60,10%/50–99 at20,25%/250–499 at30. Stop each top-up when its eligible pool has50 target trees, with a10000-attempt safety limit. Additional eligible non-target trees remain available. Randomly select50 without replacement per group from the combined pool. Each selected seed is unique. Record source, duration and every new proposal; this is an adaptively allocated, size-selected benchmark, not an unconditional simulation sample.

Use DivExtra make.tree.classe.td via the same local unsplit-living-stem safeguard as earlier pilots. No installed package changes. The2000 living-lineage cap is path-dependent selection: capped histories have unknown final sizes. Save selected full tables plus reconstructed trees. Diagnostics include boundary crossings, old/recent reconstructed nodes and branch lengths, crown age and tip-state counts.

## Inference

650 trees x three starts =1950 primary fits. Each tree is fitted only at its generating boundary; no relabeling trees across boundaries. Classical GeoSSE is not an equivalent model here and is not fitted. ROOT.OBS (observed-root treatment), condition.surv=TRUE, full sampling; tol=1e-10,eps=1e-8. Same log-scale nlminb settings and parameter bounds as the constant-rate benchmark. Three fixed jittered starts not centered on truth. Primary fits and near-zero-boundary diagnostics saved separately.

Roihu array:130 tasks, five trees/15 fits each, four R workers per task, at most six concurrent tasks;8GB,three-hour task walltime. Twenty-minute process timeout per fit including any boundary diagnostic. Completed result files permit resumption; errors/timeouts are recorded rather than silently omitted. No automatic resubmission.

Local validation checks13 representative strata against explicitly expanded six-rate likelihoods, perturbations of each free parameter, manifest uniqueness/counts and crossing identities, and three optimization starts on one small tree. First cluster task is a smoke test; remaining tasks depend on its successful completion.

## Files and execution

Run `Rscript generate.R` locally, then `Rscript validate.R`. Generate refuses to overwrite manifest.csv. In the remote directory /scratch/project_2004996/geosse-benchmark/roihu-speciation-size-50, use submit.sh. Datasets are already generated locally; cluster runs inference only. Avoid duplicate submission; consult preparation/submission-status.md.

manifest.csv identifies all selected trees, source, seed, tip counts, crossings and tasks. design.csv lists the13 included groups. preparation/new-proposals.csv and eligible-pool.csv document top-ups and selection. Source pilots retain their original proposal histories. YAML comments are preserved from the earlier benchmark.

## Analysis after completion

Report completeness, convergence, three-start agreement, rate-specific bias/error and zero-boundary diagnostics. Associate recovery with tip number and crossing counts while accounting for duration and boundary; they are correlated. Treat absent small-tree25% cells as missing by design, not successful recovery. There are no repeated equivalent-model fits here. Conditional selection and correct-model/full-sampling assumptions limit generalization. Empirical replicate ranges are not per-fit95% interval coverage. No interval coverage claim without individual-fit intervals.
