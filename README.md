# NIRS-EEG
Codes for paper Enhanced spatio-temporal reconstruction of neuronal activities using joint electroencephalography and functional near-infrared spectroscopy

*PrepFEM2.m*: generates the forward models of EEG and DOT

*PrepFEM2_Dense.m*: same as previous, but generates a denser mesh

*Fwd_NIRS_HD.m*: high-density DOT model, as in Fig. 16

*EEG_nirs_prior.m*: main script to simulate the reconstruction

*search_param.m*: scanning the parameters a and b; to generate data used in Fig. 15

*systematic_testing.m*: for systematically testing the efficacy of the algorithm; to generate data used in Fig. 17

## Requires:
- Fieldtrip
- NIRSFAST
- iso2mesh
