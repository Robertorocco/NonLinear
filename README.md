# Applying Nonlinear Control to Goodwin Oscillators for Circadian Rhythm Tracking

This repository contains the research paper and supplementary materials for the project **“Applying Nonlinear Control to Goodwin Oscillators for Circadian Rhythm Tracking in Neurospora.”** The study explores the regulation of biological oscillatory systems using advanced control theory.

---

## 📖 Overview

This project investigates the application of various control strategies to synchronize the natural oscillations of a biological system with a desired reference signal, a process known as **entrainment**. We use the Goodwin oscillator, a well‑established mathematical model, to represent the circadian rhythms of the fungus *Neurospora*.

The core of this work is a comparative analysis of one linear and two nonlinear control techniques to regulate the oscillator's behavior. The performance of each controller is rigorously evaluated based on its ability to:

- **Track a sinusoidal reference:** Make the system’s output follow a desired periodic pattern.
- **Reject external disturbances:** Maintain performance despite unexpected external inputs.
- **Handle parametric uncertainties:** Remain stable and effective even when the model’s parameters are not precisely known.

---

## 🔬 System Model: The Goodwin Oscillator

We focus on a three‑dimensional Goodwin oscillator model, which describes the concentration of:

- **mRNA** (`x`)
- **Protein** (`y`)
- **Transcriptional inhibitor** (`z`)

The dynamics are governed by:

\[
\begin{aligned}
\dot{x} &= \frac{1 + z^n}{\alpha} - x + \alpha\,u,\\
\dot{y} &= x - y,\\
\dot{z} &= y - z.
\end{aligned}
\]


where:

- \(u\) is the external control input (representing light exposure).
- \(\alpha\) is a system parameter.
- \(n\) is the Hill coefficient, crucial for generating stable oscillations.

---

## ⚙️ Control Strategies Implemented

Three distinct control architectures were designed, simulated, and compared:

1. **Proportional–Integral (PI) Control**  
   A linear strategy based on a linearized version of the Goodwin model. Effective only within a small region around the equilibrium.

2. **Input–Output Feedback Linearization (FBL)**  
   A nonlinear technique that algebraically transforms the system to a linear input–output relationship, allowing a linear controller to achieve precise tracking.

3. **Sliding Mode Control (SMC)**  
   A robust nonlinear method that forces the system’s trajectories onto a predefined sliding surface, delivering excellent disturbance rejection.

---

## 📊 Key Findings

The performance of the FBL and SMC controllers was systematically compared across three scenarios:

| Scenario       | Controller | Input Norm | Settling Time (s) | Steady‑State Error |
|----------------|------------|------------|------------------|--------------------|
| **Tracking**   | FBL        | 0.0064     | 0.27             | 0                  |
|                | SMC        | 0.0064     | 1.30             | 0                  |
| **Disturbance**| FBL        | 0.0064     | 0.27             | 0.0025             |
|                | SMC        | 0.0064     | 1.51             | 0.002              |
| **Robustness** | FBL        | 0.0092     | —                | 0.25               |
|                | SMC        | 0.06       | —                | —                  |

**Conclusion:**  
The Feedback Linearization (FBL) controller generally demonstrated superior performance, offering faster convergence and better robustness to parameter variations. While Sliding Mode Control (SMC) can perfectly reject matched disturbances using a signum function, it was more sensitive to parameter changes and exhibited a slower response in the tested scenarios.

---

## 📜 Read the Full Paper

For a detailed explanation of the mathematical models, control design, and simulation results, please read the full paper:

➡️ [View Paper.pdf](./Paper.pdf)
