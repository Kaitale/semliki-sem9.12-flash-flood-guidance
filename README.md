# Soil-moisture-dependent flash flood guidance for ungauged escarpment tributaries in an operational HEC-RTS forecasting system: the Semliki Basin, Uganda

**Kakooza, C. and Nafula, M. (2026).** KAI Water & Energy, Kampala, Uganda. Zenodo DOI: *to be assigned on publication*.
Corresponding author: Monica Nafula, monicanafula@kaiwaterenergy.com

The paper extends an operational main-river flood forecasting system for the Semliki River (HEC-HMS and HEC-RAS run under HEC-RTS with four forecast members) to six ungauged escarpment tributaries, and derives flash flood guidance directly from the soil state of the operational model at each issue.

## What is here

| File | Content |
|---|---|
| `Kakooza_et_al_2026_Semliki_flash_flood_guidance.pdf` | the paper (13 pages, 6 figures, 6 tables) |
| `reproduce_paper.ipynb` | recomputes Tables 2–7, the sensitivity of the guidance to the constant rate and the fast response, and the data figures from the data; runs in under a minute, locally or in [Google Colab](https://colab.research.google.com/github/Kaitale/semliki-sem9.12-flash-flood-guidance/blob/main/reproduce_paper.ipynb) |
| `semliki_ffg.py` | flash flood guidance from HEC-HMS states (Deficit-and-Constant loss, Clark unit hydrograph, linear-reservoir baseflow) |
| `data_tributary_points.zip` | catchments; catchment-mean daily rain (bias-corrected CHIRPS 1981–2025, IMERG Late, H60B); guidance inputs (parameters and states, 20 Nov 2026); check against HEC-HMS; calibration test; the issue of 20 Nov 2026 |
| `data_workbook_SEM.9.12.xlsx` | the same data and every derived table in one workbook |

## Run

```
pip install -r requirements.txt
jupyter notebook reproduce_paper.ipynb      # Run all
```

The notebook unzips `data_tributary_points.zip` (or downloads it from this repository) and writes every table and figure to `outputs/`.

## Data policy

All rainfall is from public products (CHIRPS, IMERG, EUMETSAT H SAF H60B). The HEC-HMS model, the staff-gauge records and the district flood register are not distributed here; requests may be addressed to the corresponding author. The notebook needs only the parameters and states carried in `ffg_inputs_20261120.json`.

## Licence

Manuscript, data and figures: CC BY 4.0. Code and notebook code cells: MIT. See `LICENSE`.

## Cite

Kakooza, C., Nafula, M., 2026. Soil-moisture-dependent flash flood guidance for ungauged escarpment tributaries in an operational HEC-RTS forecasting system: the Semliki Basin, Uganda. Zenodo. https://doi.org/10.5281/zenodo.[to be assigned]
