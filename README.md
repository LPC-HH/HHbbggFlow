# HHbbggFlow

This repo provides the centralized location for the HHbbgg analysis framework. To run the analysis call:
```
python HHbbgg_flow/analysis_manager/run_analysis.py
```
or
```python
import HHbbgg_flow as hbg

hbg.analysis_runner(**kwargs)
```

ML Training Statement:
All ML models were trained using identical Train(+Validation) / Test splitting to avoid any potential data leakages.
The models are trained on Monte Carlo simulations, using only the *odd* event numbers. The models are tested on *even* event numbers.

This repo was developed and is maintained by the Fermilab-Purdue-Caltech collaboration.
