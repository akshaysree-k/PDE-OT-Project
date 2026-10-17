# CAPSTONE PROJECT REPORT

# From Fields to Circuits:
## Locating a Transmission-Line Fault from a TDR Trace
### Using Telegrapher’s Equations, Laplace Transforms and ℓ1 Linear Programming

**Submitted by – Group 2**

- Akshay R
- Akshaysree Kalarickal
- Benson K Binu
- Hargovind KR
- Gautham A

**Department:** Electronics and Communication Engineering

**Institution:** Saintgits College of Engineering

**Course:** PDEs, Transforms & Optimization Techniques

**Semester:** 3 EC

**Academic Year:** 2025–2026

**Date:** October 2026

---

# Abstract

Time-domain reflectometry (TDR) detects discontinuities in a transmission line by launching a fast electrical step and observing reflected changes at the input. This project studies the inverse problem: estimating the location and approximate spatial profiles of resistance R(x) and capacitance C(x) from a noisy input trace.

The supplied computational setup models a 1 m line with 50 spatial nodes and nominal parameters R₀ = 0.1 Ω/m, L₀ = 250 nH/m, G₀ = 0, and C₀ = 100 pF/m. A synthetic fault is placed near x = 0.5 m, where resistance increases and capacitance decreases over three nodes.

The forward model is based on the telegrapher’s equations and their Laplace-domain representation. A finite-difference time-domain simulation produces the nominal and perturbed TDR responses; a sensitivity matrix then maps parameter perturbations to trace changes.

The inverse step is formulated as an ℓ₁ linear program, with non-negativity constraints and slack variables that express absolute residuals linearly. The supplied result plots show a measured-trace bump near the expected round-trip feature and a reconstructed resistance peak and capacitance dip near the fault.

The reconstruction also contains a secondary resistance artefact, so the plots support approximate fault localization but do not establish exact recovery of every parameter. The report documents the mathematical pipeline, computational setup, results, limitations, and verification checks needed for a reproducible implementation.

**Keywords:** time-domain reflectometry, transmission line, telegrapher’s equations, inverse problem, finite differences, ℓ₁ optimization.

---

# Contents

1. Introduction
2. Related Work
3. Methodology
   - 3.1 Field model and telegrapher’s equations
   - 3.2 Laplace-domain model
   - 3.3 Discretization and sensitivity model
   - 3.4 ℓ₁ linear-program formulation
4. Experimental Setup
5. Results and Interpretation
   - 5.1 Input TDR trace
   - 5.2 Reconstructed spatial profiles and residual
6. Limitations
7. Conclusion and Future Scope
8. Pipeline Reflection
9. References

---

# Chapter 1: Introduction

A transmission line is not always electrically uniform. A damaged cable, a connector discontinuity, or a local change in geometry can alter its distributed resistance and capacitance. When a sufficiently fast voltage step travels along the line, a discontinuity changes the local impedance and generates a reflected wave. TDR uses the time and magnitude of this reflection to diagnose signal-integrity problems.

The problem addressed here is to infer the hidden fault from the input trace rather than merely simulate a known line. The project brief specifies a one-metre line with 50 nodes and a fault near the midpoint. In the fault region, R is multiplied by four and C by 0.6 over three nodes. The measured trace is synthetic and contains small noise. The goal is to estimate the spatial profiles R(x) and C(x) and assess whether the anomaly is localized correctly.

The solution follows a connected mathematical pipeline. Vector calculus and circuit conservation lead to the telegrapher’s equations; eliminating current gives a damped-wave partial differential equation (PDE); a Laplace transform expresses propagation and impedance in the complex-frequency domain; finite differences produce a numerical forward model; and linear programming estimates sparse parameter changes from the trace residual.

This connection is important because optimization is not an unrelated post-processing step: it is derived from the discrepancy between the measured and predicted waveforms.

---

# Chapter 2: Related Work

Classical transmission-line theory describes voltage and current propagation using distributed resistance, inductance, conductance, and capacitance. Pozar presents the transmission-line equations, characteristic impedance, and reflection coefficient used to understand impedance discontinuities [1]. Paul develops a systematic treatment of multiconductor transmission lines and transient behaviour, including the role of distributed circuit parameters [2].

These theories explain the forward response when line properties are known; by themselves, they do not uniquely recover an unknown spatial profile from one noisy waveform.

TDR is a well-established diagnostic method: the delay of a reflected feature indicates distance, while its amplitude depends on the impedance change. In the ideal single-discontinuity case, the voltage reflection coefficient is

$$\Gamma =\frac{Z_{\text{fault}}-Z_{\text{ref}}}{Z_{\text{fault}}+Z_{\text{ref}}}\tag{2.1}$$

where $Z_{\text{ref}}$ is the impedance immediately before the discontinuity and $Z_{\text{fault}}$ is the impedance after it.

In a distributed, lossy line, multiple reflections and parameter coupling make this interpretation less direct.

Inverse problems are often ill-conditioned: different parameter profiles can produce similar observations. Sparse regularization and ℓ₁ objectives are widely used when only a small number of locations are expected to deviate from a nominal model. Boyd and Vandenberghe provide the general convex-optimization framework for expressing absolute-value residual minimization as a linear program [4].

The present project combines these ideas by linearizing the trace response around a nominal line and solving a non-negative ℓ₁ fitting problem.

The remaining gap addressed by this small study is therefore computational rather than a claim that TDR itself is new. The work demonstrates an end-to-end pipeline from physical equations to a reproducible inverse estimate.

Its main caveat is identifiability: a single trace may not cleanly distinguish resistance and capacitance changes, particularly on a coarse grid and in the presence of noise.

---

# Chapter 3: Methodology

## 3.1 Field Model and Telegrapher’s Equations

For a short line segment of length $dx$, voltage and current change because of the distributed series resistance $R(x)$ and inductance $L(x)$, and the shunt conductance $G(x)$ and capacitance $C(x)$.

Applying voltage and current conservation to the segment gives the telegrapher’s equations:

$$
\frac{\partial V(x,t)}{\partial x}=-R(x)I(x,t)-L(x)\frac{\partial I(x,t)}{\partial t}\tag{3.1}
$$

$$
\frac{\partial I(x,t)}{\partial x}=-G(x)V(x,t)-C(x)\frac{\partial V(x,t)}{\partial t}\tag{3.2}
$$

Here:

- $V(x,t)$ is voltage in volts.
- $I(x,t)$ is current in amperes.
- $x$ is distance along the line.
- $t$ is time.

The per-unit-length parameters have units:

- Resistance: Ω/m
- Inductance: H/m
- Conductance: S/m
- Capacitance: F/m

The spatial derivatives represent voltage or current change per unit length; the time-derivative terms represent inductive and capacitive effects.

These equations are the one-dimensional transmission-line form consistent with Maxwell’s field laws and circuit conservation for a line whose cross-sectional fields can be represented by distributed parameters.

For constant parameters, differentiate Equation (3.1) with respect to $x$ and substitute Equation (3.2) and its time derivative to eliminate $I$.

The resulting voltage equation is

$$
\frac{\partial^2 V}{\partial x^2}=LC\frac{\partial^2 V}{\partial t^2}+(RC+LG)\frac{\partial V}{\partial t}+RGV\tag{3.3}
$$

The analogous current equation has the same structure.

Equation (3.3) is a damped-wave, hyperbolic PDE when the line parameters are positive and the distributed coefficients are constant.

The initial state is taken to be zero before the source step is applied. At $x=0$, the source launches the incident step; at $x=\ell$, the load condition determines the outgoing and returning waves.

A physically complete simulation must specify the source impedance and load impedance, because both affect the trace.

## 3.2 Laplace-Domain Model

Taking the Laplace transform with respect to time, assuming zero initial conditions, converts time derivatives into multiplication by the complex frequency $s$.

The transformed telegrapher’s equations are

$$\frac{dV(x,s)}{dx}=-(R+sL)I(x,s)\tag{3.4}$$

$$\frac{dI(x,s)}{dx}=-(G+sC)V(x,s)\tag{3.5}$$

The propagation constant and characteristic impedance follow as

$$\gamma(s)=\sqrt{(R+sL)(G+sC)}\tag{3.6}$$

$$Z_c(s)=\sqrt{\frac{R+sL}{G+sC}}\tag{3.7}$$

The propagation constant $\gamma(s)$ describes attenuation and phase change per unit length; $Z_c(s)$ is the voltage-to-current ratio of a travelling wave.

For a low-loss line, the nominal wave speed is approximately

$$v=\frac{1}{\sqrt{LC}}$$

With

$$L_0=250\text{ nH/m}$$

and

$$C_0=100\text{ pF/m}$$

the wave speed is approximately

$$v\approx 2.0\times 10^8\text{ m/s}$$

This gives a one-way delay of about 5 ns over one metre.

A reflection returning to the input from a fault at 0.5 m is therefore expected after approximately

$$t_{\text{round-trip}}=\frac{2x}{v}=5\text{ ns}$$

Using the nominal values, the characteristic impedance is approximately

$$Z_{c0}\approx \sqrt{\frac{L_0}{C_0}}=50\ \Omega$$

If only capacitance changes to $0.6C_0$, the local low-loss impedance estimate becomes

$$Z_{\text{fault}}\approx\sqrt{\frac{L_0}{0.6C_0}}\approx 64.5\ \Omega$$

Substitution into Equation (2.1) gives

$$\Gamma\approx 0.127$$

This is a useful order-of-magnitude check for the observed TDR bump, not an exact prediction for the full lossy, multi-parameter simulation.

## 3.3 Discretization and Sensitivity Model

The line is divided into

$$N_x=50$$

nodes across

$$\ell=1\text{ m}$$

The grid spacing is

$$\Delta x=\frac{\ell}{N_x-1}=\frac{1}{49}\text{ m}\approx 20.4\text{ mm}\tag{3.8}$$

The forward model uses a staggered-grid leapfrog finite-difference time-domain (FDTD) scheme.

Voltage and current are updated at alternating spatial and temporal locations, which is natural for the coupled first-order Equations (3.1)–(3.2).

For wave speed $v$, a conservative Courant condition is

$$\frac{v\Delta t}{\Delta x}\leq 1\tag{3.9}$$

The project configuration uses

$$\Delta t=0.9\frac{\Delta x}{v}\approx 0.092\text{ ns}$$

which is below this limit.

Approximately 435 steps cover the stated simulation window of roughly 40 ns.

Let $y_{\text{meas}}$ be the sampled measured TDR trace and $y_0$ the trace generated by the nominal line. Their difference is the data vector

$$b=y_{\text{meas}}-y_0\tag{3.10}$$

Let $p$ collect the discretized parameter perturbations.

A sensitivity matrix $A$ is assembled so that its $j$-th column represents the approximate change in the TDR samples caused by a perturbation at node $j$.

For small perturbations, the nonlinear forward response is approximated locally by

$$Ap\approx b\tag{3.11}$$

This is a linearized inverse model, not a guarantee that the complete transmission-line equations are linear in the unknown parameters.

Because a fault is expected to occupy only a few nodes, a sparse residual-fitting objective is appropriate.

## 3.4 ℓ₁ Linear-Program Formulation

The inverse estimate minimizes the sum of absolute data residuals:

$$\min_{p\geq 0}\|Ap-b\|_1=\min_{p\geq 0}\sum_{k=1}^{M}|(Ap-b)_k|\tag{3.12}$$

where $M$ is the number of sampled time points and $p\geq 0$ enforces non-negativity for the chosen perturbation representation.

Absolute values are not linear, so introduce a non-negative slack variable $\varepsilon_k$ for each sample.

The condition

$$|(Ap-b)_k|\leq\varepsilon_k$$

is equivalent to two linear inequalities:

$$Ap-\varepsilon\leq b\tag{3.13}$$

$$-Ap-\varepsilon\leq-b\tag{3.14}$$

with constraints

$$\varepsilon\geq 0,\qquad p\geq 0\tag{3.15}$$

The resulting linear program is

$$\min_{p,\varepsilon}\mathbf{1}^{T}\varepsilon\tag{3.16}$$

subject to the inequalities above.

The objective $\mathbf{1}^{T}\varepsilon$ is the total absolute residual. It is linear, as are the constraints, so a standard LP solver such as SciPy’s `linprog` with the HiGHS backend can solve it.

Since resistance is expected to rise while capacitance falls, the implementation must encode signed changes carefully, for example by representing positive and negative perturbations with separate non-negative variables.

Non-negativity of the final physical profiles $R(x)$ and $C(x)$ should be checked after reconstruction.

---

# Chapter 4: Experimental Setup

The supplied project presentation specifies Python as the implementation language.

The following libraries and tools are used:

- **NumPy:** Arrays and numerical calculations.
- **SciPy:** Linear-program solver and sparse/numerical routines.
- **SymPy:** Symbolic checks.
- **Matplotlib:** Plots.
- **Pytest:** Automated checks.
- **Git:** Reproducible code history.

The proposed repository separates configuration, transmission-line theory, the forward solver, sensitivity construction, inversion, metrics, and plotting.

The baseline configuration is summarized in Table 4.1.

The provided figures are output plots from the project workflow. The raw trace samples, complete sensitivity matrix, solver log, and numerical metrics file were not included with the supplied materials; consequently, this report interprets the visible plots and the stated setup without inventing additional precision or claiming that unprovided tests were independently rerun.

### Table 4.1: Transmission-Line Configuration

| Parameter | Value |
|---|---|
| Line length, $\ell$ | 1 m |
| Spatial nodes, $N_x$ | 50 |
| Nominal resistance, $R_0$ | 0.1 Ω/m |
| Nominal inductance, $L_0$ | 250 nH/m |
| Nominal conductance, $G_0$ | 0 S/m |
| Nominal capacitance, $C_0$ | 100 pF/m |
| Fault location | Approximately 0.5 m |
| Fault width | 3 nodes |
| Fault changes | $R\times4,\ C\times0.6$ |
| Nominal impedance | Approximately 50 Ω |
| Nominal wave speed | Approximately $2\times10^8$ m/s |

---

# Chapter 5: Results and Interpretation

## 5.1 Input TDR Trace

Figure 5.1 overlays the nominal-line response and the noisy measured trace.

Both rise rapidly from zero to a plateau near 0.5 V.

The most visible difference is a positive bump around 5–5.5 ns followed by a small undershoot.

This timing is consistent with the predicted round-trip delay from a fault near the midpoint of a one-metre line with wave speed close to $2\times10^8$ m/s.

A smaller feature appears close to the end of the displayed window; without the source/load configuration and raw samples, it should not be attributed uniquely to another fault.

**Figure 5.1:** Provided synthetic TDR trace: nominal response and noisy measured response. The fault-related deviation is visible near the expected midpoint round-trip delay.

> [Original PDF Figure 5.1: Insert the nominal and noisy measured TDR trace plot here.]

## 5.2 Reconstructed Spatial Profiles and Residual

Figure 5.2 shows the estimated resistance and capacitance profiles together with the reference profiles, as well as the waveform residual.

The dashed reference resistance profile has a localized increase around $x=0.5$ m, while the reference capacitance falls in the same region.

The LP reconstruction places its largest resistance spike and its capacitance dip near this location.

This is the main qualitative success: both reconstructed profiles identify the approximate fault region rather than spreading the anomaly uniformly along the line.

The reconstruction is not exact.

The plotted resistance estimate overshoots the reference peak and contains a secondary spike around

$$x\approx0.77\text{ m}$$

that is absent from the reference.

This secondary feature is consistent with an inversion artefact caused by noise, limited resolution, parameter coupling, or an imperfect linearized sensitivity model.

The capacitance dip is more spatially concentrated and appears close to the reference location, but the coarse grid limits how precisely the fault edges can be estimated.

The residual panel oscillates around zero with most visible excursions on the order of a few hundredths of the plotted millivolt scale.

A residual centred around zero suggests that the fitted trace captures the broad measured waveform without a large constant bias.

However, a small residual alone does not prove that the recovered $R(x)$ and $C(x)$ are uniquely correct: distinct parameter profiles may produce similar traces.

Numerical residual norms and a noise-sweep table should be exported directly from the solver before quantitative error claims are made.

**Figure 5.2:** Provided inversion results: reconstructed $R(x)$ and $C(x)$ against the reference profiles, and the residual $V_{\text{meas}}-V_{\text{fit}}$.

> [Original PDF Figure 5.2: Insert the reconstructed resistance profile, capacitance profile, reference profiles, and residual plot here.]

### Table 5.1: Results Supported by the Supplied Plots and Configuration

| Measure | Interpretation |
|---|---|
| Expected fault position | Approximately 0.5 m; three-node region. |
| Expected reflection timing | Approximately 5 ns round-trip delay from the midpoint. |
| TDR observation | A visible positive deviation occurs around 5–5.5 ns, consistent with the expected timing. |
| Resistance reconstruction | Main spike is near the reference fault but overshoots; a secondary non-reference spike is visible. |
| Capacitance reconstruction | A localized dip appears near the reference fault region. |
| Residual | Oscillates around zero; a numeric norm or RMS value is not available in the supplied materials. |
| Quantitative accuracy | Cannot be certified from screenshots alone; raw arrays and solver outputs are required. |

The project presentation lists further verification targets:

1. No internal reflection for a uniform line.
2. Agreement between the finite-difference impedance/delay and theory.
3. Fault recovery within two grid cells.
4. Stability under a noise sweep.
5. Non-negative reconstructed profiles.
6. A successful LP solver status.

These are appropriate acceptance checks, but the presentation describes them as checks to perform and does not include machine-readable logs for each one.

They should therefore be reported as verification targets unless the team can attach the corresponding test output.

---

# Chapter 6: Limitations

The inverse problem is underdetermined in practice because resistance and capacitance changes can affect the same trace features.

A single input trace may not separate their effects reliably, especially when the parameter perturbations are estimated from a first-order sensitivity model.

The 50-node spatial grid has a spacing of about 20.4 mm, which limits localization and edge resolution.

Noise and multiple reflections can create correlated columns in the sensitivity matrix, leading the LP to place some anomaly energy at a false location, as suggested by the secondary resistance spike.

The displayed results are screenshots rather than raw datasets. This prevents independent computation of:

- Peak-location error.
- Relative profile error.
- Residual RMS.
- Confidence intervals.

The report also cannot independently verify source and load boundary settings, solver convergence status, the exact noise realization, or the claimed runtime.

A stronger evaluation would include multiple noise seeds, a sweep of regularization or weighting choices, a comparison with a nonlinear forward-model refinement, and a table of metrics calculated from the saved arrays.

---

# Chapter 7: Conclusion and Future Scope

This project links electromagnetic transmission-line modelling to numerical inversion.

The telegrapher’s equations describe the forward dynamics; the Laplace-domain propagation constant and characteristic impedance explain the timing and scale of reflections; finite differences provide the simulated response; and an ℓ₁ linear program estimates localized parameter perturbations from the measured-minus-nominal trace.

In the supplied figures, the TDR feature occurs near the expected midpoint round-trip delay, and the reconstructed resistance and capacitance profiles show a spike and dip near the reference fault.

These results support approximate localization, but the resistance overshoot and secondary spike show that the profile recovery is not fully accurate.

## Future Scope

Future work should save and publish:

- The raw measured, nominal, and fitted traces.
- The sensitivity matrix.
- The LP solver status.
- Numerical metrics such as fault-location error, profile error, and residual RMS.

Testing multiple noise realizations and using multiple source waveforms or measurements at both ends of the line could improve identifiability.

A final nonlinear refinement around the LP estimate may also reduce bias introduced by the local linearization.

---

# Pipeline Reflection

The project uses vector calculus and physical conservation to derive the field model, PDE methods to classify the propagation equation, Laplace transforms to obtain the frequency-domain line quantities, finite differences and linear algebra to build the forward/sensitivity system, and optimization to estimate the fault.

Python joins these stages into a testable pipeline.

The main lesson is that each stage constrains the next: a physically inconsistent forward model cannot be repaired simply by using a more powerful optimizer.

### Table 7.1: Mathematical Pipeline and Implementation Mapping

| Stage | Course Concept | Equation or Implementation |
|---|---|---|
| Field model | Vector calculus / electromagnetics | Telegrapher equations (3.1)–(3.2) |
| PDE | Partial differential equations | Damped-wave equation (3.3) |
| Transform | Laplace transforms | $\gamma(s)$ and $Z_c(s)$ in (3.6)–(3.7) |
| Discretization | Numerical methods / linear algebra | Staggered-grid FDTD and sensitivity system $Ap\approx b$ |
| Optimization | Linear programming | Slack-variable formulation (3.16) |
| Verification | Python computational lab | Trace/profile plots, CFL and non-negativity checks |

---

# References

## Bibliography

[1] D. M. Pozar, *Microwave Engineering*, 4th ed. Hoboken, NJ: John Wiley & Sons, 2012.

[2] C. R. Paul, *Analysis of Multiconductor Transmission Lines*, 2nd ed. Hoboken, NJ: John Wiley & Sons, 2008.

[3] W. H. Hayt, J. A. Buck, and M. J. Kemmerly, *Engineering Electromagnetics*, 8th ed. New York: McGraw-Hill, 2012.

[4] S. Boyd and L. Vandenberghe, *Convex Optimization*. Cambridge, UK: Cambridge University Press, 2004.

[5] R. Tibshirani, “Regression Shrinkage and Selection via the Lasso,” *Journal of the Royal Statistical Society: Series B*, vol. 58, no. 1, pp. 267–288, 1996.

[6] P. Virtanen et al., “SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python,” *Nature Methods*, vol. 17, pp. 261–272, 2020, doi: 10.1038/s41592-019-0686-2.
