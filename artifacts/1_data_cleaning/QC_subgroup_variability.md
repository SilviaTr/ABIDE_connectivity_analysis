# QC – Subgroup variability

- Total subjects (with .1D): **884**

- Counts by diagnosis:

| DX_GROUP   |   count |
|:-----------|--------:|
| TDC        |     476 |
| ASD        |     408 |


## Sex × Diagnosis

| DX_GROUP   |   F |   M |
|:-----------|----:|----:|
| ASD        |  50 | 358 |
| TDC        |  88 | 388 |

- Δ male proportion (ASD vs TDC) = **0.062** → OK (threshold 0.15)


## Site × Diagnosis (top 10)

| DX_GROUP   |   NYU |   UM_1 |   USM |   UCLA_1 |   YALE |   PITT |   TRINITY |   MAX_MUN |   KKI |   CALTECH |
|:-----------|------:|-------:|------:|---------:|-------:|-------:|----------:|----------:|------:|----------:|
| ASD        |    73 |     36 |    38 |       28 |     22 |     22 |        21 |        18 |    12 |        19 |
| TDC        |    98 |     46 |    23 |       27 |     26 |     23 |        23 |        24 |    27 |        18 |

- ⚠️ **3** cells < **10** in SITE×DX (imbalance / low counts).


## Age by diagnosis – summaries

| DX_GROUP   |   count |    mean |     std |   median |   min |   max |
|:-----------|--------:|--------:|--------:|---------:|------:|------:|
| ASD        |     408 | 17.6904 | 8.93392 |   15     |  7    |  64   |
| TDC        |     476 | 16.7892 | 7.34826 |   14.755 |  6.47 |  56.2 |


## Age – tests

- Welch t = **1.62**, p = **0.105**

- Mann–Whitney U = **100171.0**, p = **0.41778346368540054**

- KS = **0.04971988795518207**, p = **0.6289978813217064**

- Effect sizes: SMD ≈ **0.11**, Cliff’s δ ≈ **-0.03**


## Age bins × Diagnosis (χ²)

| DX_GROUP   |   ≤12 |   13–18 |   19–30 |   31–50 |   ≥51 |
|:-----------|------:|--------:|--------:|--------:|------:|
| ASD        |   110 |     167 |      92 |      34 |     5 |
| TDC        |   124 |     208 |     117 |      26 |     1 |
