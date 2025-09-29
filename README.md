# Reduced Order Modeling with Large Deviation Theory

This repository provides code for generating **Reduced Order Models (ROMs)** for PDE systems with random parameters, with a focus on **accurate prediction of rare tail events**.  

The approach combines:
- **Proper Orthogonal Decomposition (POD)** for reduced basis construction  
- **Large Deviation Theory (LDT)** to identify the *LDT optimizer* (the most likely parameter leading to a tail event)  
- **Shifted sampling**: instead of sampling around the original mean, samples are drawn around the optimizer  
- **Importance sampling weights**: used to ensure samples remain consistent with the original distribution  
- **Tail-focused weights**: further bias the snapshot matrix so that extreme (tail) samples contribute more to the reduced model  

---

## Repository Structure

- `supplement_functions.jl`  
  Core functions implementing POD, weighting, sampling, and ROM construction  

- `external/`  
  Example setup borrowed from [DPFEHM](https://github.com/lanl/DPFEHM):  
  steady-state 2D Darcy flow with a random Gaussian log-permeability field (via KL expansion), including a **caprock with vertical fractures**  

- `brute_force_sampling.jl`  
  Generates **10,000 brute force samples**, used for:  
  - Identifying tail events  
  - Building a reference test set  

- `reduced_model.jl`  
  Runs the **LDT-based sampling** and constructs the Reduced Order Model using weighted POD  

---

## Installation & Setup

This project depends on [DPFEHM](https://github.com/lanl/DPFEHM).  

1. **Clone DPFEHM**  
   ```bash
   git clone https://github.com/lanl/DPFEHM.git
2. **Clone this repository** inside the same parent directory:
   ```bash
   git clone https://github.com/mxli630/Tail_Aware_POD.git
   ```

   The directory structure should look like:
   ```
   parent_folder/
   ├── DPFEHM/
   └── Tail_Aware_POD/
  ```

4. **Install Julia dependencies**
   Open Julia inside this repo and run:
   ```julia
   using Pkg
  Pkg.activate(".")
  Pkg.instantiate()
