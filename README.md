<img width="426" alt="FANCE2" src="https://github.com/nmgountanis/fanCE/assets/143458646/c59fd882-d674-4ab7-87cd-44fedde1c7a8">

# Flexible ANalytic Chemical Evolution

fanCE is an analytic one-zone chemical evolution code. 
It implements the Weinberg, Andrews, & Freudenburg (2017) analytic solution for a constant, exponential declining, 
or linear-exponential (delayed tau) star formation history. 
It also implements a flexible, two-parameter star formation history $\propto(1-e^{-t/\tau_1})\cdot e^{-t/\tau_2}$, 
which rises linearly on a timescale $\tau_1$ and then falls exponentially on a timescale $\tau_2$ for $\tau_1>\tau_2>0$.

## Installation
The source code for fanCE can be downloaded from GitHub and installed by running

    cd <path_to_installation>
    git clone https://github.com/nmgountanis/fanCE
    cd fanCE
    pip install -e .

## Authors
- Nicole Gountanis (gountanis.1@buckeyemail.osu.edu)
- Nathan Sandford
- David Weinberg
- Aliza Beverage

## License & Attribution
Copyright 2024 Nicole Gountanis and contributors.

The source code is made available under the terms of the MIT license.

If you make use of this code, please cite Gountanis et al. (2024, LINK) for the code and Weinberg et al. (2017, https://ui.adsabs.harvard.edu/abs/2017ApJ...837..183W/abstract) for the analytic solutions. The reference for the default CCSN and SNIa yield choices is Weinberg et al. (2024, https://ui.adsabs.harvard.edu/abs/2023arXiv230905719W/abstract)
