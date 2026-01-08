## Simulation Code Structure (Study 1)

`Sim_ES_fixedALL_S_upper_rand.m` and `EvidenceAccumulate_S_upper.m`  
are the main files used to perform evidence accumulation for the selected model

Afterwards, `plotEmpiricalData_auto.m` can be used to summarise simulated draws
and plot them individually to inspect variance

Finally, `combined_plotEmpiricalData.m` can be run to compare simulated and real
data across all draws.

**Posterior draw files:**
- `garcia_replication_For_paper_21_posterior_draws.csv` → ESaDDM, no *z*
- `garcia_replication_For_paper_16_posterior_draws.csv` → aDDM, no *z*
- `garcia_replication_For_paper_7_posterior_draws.csv` → ESaDDM + *z*, *a* ~ OV
- `garcia_replication_For_paper_6_posterior_draws.csv` → ESaDDM + *z*
- `garcia_replication_For_paper_2_posterior_draws.csv` → aDDM + *z*

---

## Simulation Code Structure (Study 2)

`Sim_ES_fixedALL_S_upper_rand.m` and `EvidenceAccumulate_S_upper.m`  
are the main files used to perform evidence accumulation for the selected model

Afterwards, `plotEmpiricalData_auto.m` can be used to summarise simulated draws
and plot them individually to inspect variance

Finally, `combined_plotEmpiricalData.m` can be run to compare simulated and real
data across all draws.

**Posterior draw files:**
- `OV_replication_For_paper_11_posterior_draws.csv` → ESaDDM, no *z*
- `OV_replication_For_paper_8_posterior_draws.csv` → aDDM, no *z*
- `OV_replication_For_paper_6_posterior_draws.csv` → ESaDDM + *z*, *a* ~ OV
- `OV_replication_For_paper_5_posterior_draws.csv` → ESaDDM + *z*
- `OV_replication_For_paper_3_posterior_draws.csv` → aDDM + *z*
