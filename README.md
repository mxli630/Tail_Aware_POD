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

## 1) Model Summary (What is being modeled?)

- **Physical model (forward code):** porous-media flow (e.g., steady Darcy), solved by your application code (external repo / file).
- **Random medium:** the **log-permeability** field is modeled as a **Gaussian Random Field (GRF)** and **parameterized by a truncated Karhunen–Loève (KL) expansion**.
- **KL parameterization (truncated to `M` modes):**

  $\log K(z) = \mu + \sum_{i=1}^M \sqrt{\lambda_i} x_i \phi_i(z),$
  where
  - $\mu$ is the mean log-permeability,
  - $(\lambda_i, \phi_i)$ are KL eigenpairs of the GRF covariance,
  - $x = (x_1,\dots,x_M)$ are the **independent parameters**.

- **Parameter prior:** by default we take the **KL coefficients $x_i$ to be independent standard normals**,  
  $x \sim \mathcal N(0, I_M).$

- **Quantity of Interest (QoI):** a scalar functional extracted from each forward solve (examples: breakthrough time at a sensor, max hydraulic head over the domain, flux through an outflow boundary, etc.). In our example, we set the QoI as the pressure head at a critical point.

---

## 2) LDT Ingredients You Must Define

### 2.1 Rate Function $I(\cdot)$
The **rate function** encodes the prior on parameters and governs the tail-biased sampling:

- **If the parameters are the standardized KL coefficients $x \sim \mathcal N(0, I)$:**
  $I(x) = \tfrac12 \|x\|_2^2 = \tfrac12 \sum_{i=1}^M x_i^2.$

- **If the parameters are the *scaled* coefficients \(\xi_i = s_i\sqrt{\lambda_i}\,x_i\) with \(\xi \sim \mathcal N(0, \mathrm{diag}(\sigma_1^2,\dots,\sigma_M^2))\):**
  \[
  I(\xi) \;=\; \tfrac12 \sum_{i=1}^M \left(\frac{\xi_i}{\sigma_i}\right)^2.
  \]

> ⚠️ **Consistency requirement:** The exact form of \(I\) you implement **must match** the parameterization used by the forward/application code. In particular, the **mode-wise uncertainty constants** (the \(s_i\) and any additional scalings) must be **identical** between the forward model and your rate function. See §3.

**Where to edit in this repo:**  
- Define/modify the rate function in **`rate_function.jl`** (or the corresponding function used by `generate_tail_sample.jl`).  
- Search for `rate(` or `rate_function(` as an entry point if the file name differs.  
- Add comments near the constants explaining the mapping you assume (standardized \(x\) vs scaled \(\xi\)).

### 2.2 Quantity of Interest (QoI)
Define a function that maps a forward solution (and/or boundary traces) to a **scalar** QoI.

**Where to edit in this repo:**  
- Implement/adjust the QoI in **`qoi.jl`** (or the helper where QoI is currently computed during runs; often referenced by `generate_tail_sample.jl` and `brute_force_sampling.jl`).  
- Search for `qoi(` to locate it.

> Tip: Keep the QoI *pure* (no global state); document its physical meaning and units in a comment block.






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

   

