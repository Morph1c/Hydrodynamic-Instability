# Hydrodynamic Instability

This repository contains MATLAB code and resources for studying and simulating hydrodynamic instabilities, specifically centrifugal instabilities in vortex flows. The primary goal of the code is to calculate and analyze the growth rates and frequencies (represented by the real and imaginary parts of the complex eigenvalue `\sigma`) for various types of base flows, including Lamb-Oseen vortices, Q-vortices, and more complex velocity profiles such as the Carton & McWilliams model and Francis turbine profiles.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Available Base Flows](#available-base-flows)
- [File Structure](#file-structure)
- [Contributing](#contributing)
- [License](#license)

## Overview

Hydrodynamic instabilities occur when disturbances in a fluid flow grow over time, leading to the breakdown of the flow into turbulence. This repository focuses on instabilities that arise in rotating and swirling flows, where centrifugal forces can play a significant role. The code solves the governing equations using spectral methods and visualizes the resulting growth rates and frequencies for disturbances with different azimuthal modes (`m`).

## Features

- Computes growth rates and frequencies for various base flows.
- Supports a variety of base flow profiles, including Lamb-Oseen vortex, Q-vortex, Carton & McWilliams model, and Francis turbine velocity profile.
- Allows for parameter adjustment to investigate different physical conditions.
- Visualizes the real and imaginary parts of the eigenvalue `\sigma` as a function of wavenumber (`k`).
- Customizable for different azimuthal modes (`m`).

## Installation

### Prerequisites

- MATLAB (Tested on MATLAB R2023 or later)
- Symbolic Math Toolbox (for symbolic calculations)
- Basic understanding of fluid dynamics and MATLAB scripting

### Steps

1. Clone the repository:
    ```bash
    git clone https://github.com/Morph1c/Hydrodynamic-Instability.git
    cd Hydrodynamic-Instability
    ```

2. Open MATLAB and navigate to the repository directory:
    ```matlab
    cd 'path_to_cloned_repo'
    ```

3. Make sure all necessary dependencies are installed in MATLAB (like the Symbolic Toolbox).

## Usage

1. Open the main script `plot_centrifugal_criteria.m` or other scripts like `calculate_sigma_for_m.m` in MATLAB.

2. Define your base flow by selecting from the available options (Lamb-Oseen vortex, Q-vortex, etc.) within the script.

3. Customize parameters such as the azimuthal mode (`m`), step size (`delta_k`), and wavenumber range (`k_max`).

4. Run the script to compute and visualize the growth rates and frequencies for the selected base flow and parameters.

5. The script generates plots for the real and imaginary parts of the growth rate, `\sigma`, as a function of wavenumber, `k`.

### Example

```matlab
m_values = -2:-1:-12;
calculate_sigma_for_m_criteria(m, baseflow_parameters, n, delta_k, k_max, find_initial_guess);
