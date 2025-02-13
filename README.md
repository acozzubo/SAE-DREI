# SAE-DREI: Small Area Estimation for Disaggregated Race and Geography Research

## Overview
This repository contains scripts for a project analyzing best use of **Small Area Estimation (SAE) methods for sociodemographic research**. The repository is structured to facilitate the implementation of SAE methodologies in addressing research questions that involve race and geographic cross-classifications.

### Scripts repository
**This repository does not include data or proprietary research findings.** It provides only the scripts used for analysis

## Purpose of the Code
The scripts in this repository were developed to implement large **design-based simulations** using public datasets that evaluate the application of SAE techniques in demographic related research questions. Specifically, the simulations focus on:
- Disaggregating race cross-classified with geography.
- Using 2021 ACS microdata as an artificial "known population" to study the impact of different sampling designs.
- Evaluating homeownership estimates by race and state in the United States.
- Comparing different SAE modeling choices and their effectiveness through Mean Squared Error (MSE) analysis.

## Summary of Simulation Findings
The project undertook a **design-based simulation** where repeated samples were drawn under **complex sampling designs** to mimic survey uncertainty. The goal was to assess the impact of different modeling choices for SAE.

### Sampling Designs
Five distinct sampling designs were explored, each with 500 replicated samples:
1. **Equal Sample Size for All Races (eq_1):**
   - Uniform sample sizes across racial groups and states.
   - Total sample size: 30,269 observations (1% of ACS population).
2. **Equal Sample Size for All Races with GSS Size (eq_gss):**
   - Similar to eq_1 but aligned with General Social Survey (GSS) sample size.
   - Total sample size: 4,080 observations.
3. **Probability Proportional to Race Size with 1% of ACS Size (pr_01):**
   - Sample size based on probability proportional to race size.
   - Total sample size: 30,943 observations.
4. **Probability Proportional to Race Size with 0.5% of ACS Size (pr_05):**
   - Similar to pr_01 but with a reduced total sample size.
   - Total sample size: 15,487 observations.
5. **Probability Proportional to Race Size with GSS Size (pr_gss):**
   - Probability proportional to race size, with total sample size approximating the GSS.
   - Total sample size: 4,080 observations.

**Note:** Small differences in sample sizes may occur due to proportional representation adjustments.

### Modeling Considerations
- Models with and without fixed effects and covariates were tested.
- Different methods of incorporating covariates were explored.
- Separate models were fitted for racial groups to analyze their impact.
- The effect of modeling choices on **small area proportion estimates** and **inter-group differences** was evaluated.

## Repository Structure
```
├── scripts/        # SAE implementation scripts
├── README.md       # Project documentation
```

## How to Use
### Prerequisites
- **R (>=4.0)** or **Python (>=3.8)** depending on the script requirements.
- Required packages (install as needed):
  ```r
  install.packages(c("survey", "sae", "lme4", "ggplot2"))
  ```

### Running the Scripts
1. Clone the repository:
   ```bash
   git clone https://github.com/acozzubo/SAE-DREI.git
   ```
2. Navigate to the scripts folder:
   ```bash
   cd SAE-DREI/scripts
   ```
3. Run the analysis in the order of the script:

## Additional Information
For a more detailed discussion of the **methodology, simulation steps, and findings**, please refer to the full research report.

## License and Usage
- The scripts provided are intended **solely for research and methodological purposes**.

## Contact
For inquiries or collaboration requests, please contact the project maintainers through GitHub Issues or via professional channels.

## License
This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

---
🚀 **For further insights into SAE applications, explore the repository and engage with the research community!**

