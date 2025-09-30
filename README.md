# Reduced Order Modeling with Large Deviation Theory

This repository provides code for generating **Reduced Order Models (ROMs)** for PDE systems with random parameters, with a focus on **accurate prediction of rare tail events**. We use **large-deviation-theory(LDT)-guided sampling** to generate informative parameter realizations (including tail/rare events), then build the ROM from the resulting snapshot set. 

The approach combines:
- **Proper Orthogonal Decomposition (POD)** for reduced basis construction  
- **Large Deviation Theory (LDT)** to identify the *LDT optimizer* (the most likely parameter leading to a tail event)  
- **Shifted sampling**: instead of sampling around the original mean, samples are drawn around the optimizer  
- **Importance sampling weights**: used to ensure samples remain consistent with the original distribution  
- **Tail-focused weights**: further bias the snapshot matrix so that extreme (tail) samples contribute more to the reduced model  

> **Two things every user must provide/understand for a new application:**
> 1) the **Rate Function** (consistent with the parameter prior and scaling), and  
> 2) the **Quantity of Interest (QoI)** (what you measure from a forward run).

---

## Repository Structure

- `supplement_functions.jl`  
  Core functions implementing POD, weighting, sampling, and ROM construction  

- `external/`  
  Contains the setup for a **steady-state 2D Darcy flow** example with a random Gaussian log-permeability field (via KL expansion), including a **caprock with vertical fractures**.
  This example was adapted from [DPFEHM](https://github.com/lanl/DPFEHM), and the necessary file has been copied here so that the repository runs standalone.  


- `brute_force_sampling.jl`  
  Generates **10,000 brute force samples**, used for:  
  - Identifying tail events  
  - Building a reference test set  

- `generate_tail_sample.jl`  
  Runs the **LDT-based sampling**:  
  - Computes the LDT optimizer  
  - Shifts the mean to the optimizer and samples around it  
  - Applies importance and tail weights  
  - Produces weighted snapshots for ROM construction  

- `reduced_model.jl`  
  Constructs the **Reduced Order Model** from the LDT-based samples (e.g., builds the reduced POD basis and solves the reduced problem).  

---

## Installation & Setup

This project depends on [DPFEHM](https://github.com/lanl/DPFEHM).  

1. **Clone this repository**  
   ```bash
   git clone https://github.com/mxli630/Tail_Aware_POD.git
   cd Tail_Aware_POD


2. **Install Julia dependencies**
   Open Julia inside this repo and run:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()

## Usage
1. **Brute Force Sampling**
   Generate 10,000 samples for baseline and test sets:
   ```bash
   julia brute_force_sampling.jl
   ```
   **Outputs**: two ``.jld2`` files will be saved in the current directory. These contain the brute force parameters and the corresponding QoI.

2. **Tail_Aware Sampling**
   Run the LDT-based sampling:
   ```bash
   julia generate_tail_sample.jl
   ```
   **Outputs**: one ``.jld2`` file will be saved in the current directory, containing the LDT-shifted parameters, snapshots (solutions to the PDE), QoI, and IS weights.

3. **Reduced Model Construction**
   Build the reduced order model from the LDT-based samples:
   ```bash
   julia reduced_model.jl
   ```
   **Outputs**: one ``.jld2`` file will be saved in the current directory, containing a series of reduced model related errors.

   

