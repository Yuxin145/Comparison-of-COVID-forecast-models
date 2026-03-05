# COVID-19 Epidemic Forecasting: Mechanistic vs. Foundation Models

This repository contains the complete codebase and evaluation pipeline for a comparative study of short-term COVID-19 forecasting models. The project evaluates traditional mechanistic and statistical models alongside modern zero-shot foundation models. The evaluation spans across both continuous incidence trajectories and discrete peak-time prediction tasks. 

Performance is assessed using general probabilistic metrics, including CRPS, calibration, and sharpness. Furthermore, the models are evaluated using decision-theoretic measures, such as the Cost-Loss framework and Murphy diagrams, to quantify actionable forecast value for public health planning.

## Full Research Report
You can read the comprehensive methodology, experimental design, and findings in the full project report. Please navigate to the [`docs/Yuxin_urops_report.pdf`](./docs/Yuxin urops report.pdf) file to view the detailed analysis.

## Repository Architecture

To accommodate distinct framework dependencies and modeling paradigms, this repository is divided into two operational pillars.

* **`baseline_r/`**: This directory contains the statistical time series and mechanistic model pipeline. It handles hyperparameter tuning, forecast generation, and general probabilistic evaluation for models including EpiNow2, EpiEstim, SARIMA, Gaussian Process, and Prophet.
* **`foundation_python/`**: This directory contains foundation model pipeline. It handles data preprocessing, zero-shot incidence and peak-time forecasting for Chronos-2 and DeepAR, and the decision-oriented evaluation framework.

## Environment Setup

Due to specific framework requirements, the models operate across distinct environments. 

* **R Baseline Environment**: This environment requires an installation of R and RStudio. Key required packages include `EpiNow2` and `forecast`. 
* **Chronos-2 & Evaluation Environment**: 

  ```bash
  conda env create -f foundation_python/environment_chronos.yml
  conda activate covid_chronos
  ```

* **DeepAR Environment**: 

  ```bash
  conda env create -f foundation_python/environment_deepar.yml
  conda activate covid_deepar
  ```

## To Run

### Statistical Baselines (R)
* **Prepare the Data**: Download the dataset and preprocess the raw incidence data:

  ```bash
  Rscript baseline_r/data/selected_covid_data_us.r
  ```
* **Experiments**: Run `experiments/Experiment_models.r` and `experiments/experiments.rmd` to find best hyperparameters. 
* **Generate the Forecasts**: Generate incidence forecasts using the optimized model parameters:

  ```bash
  Rscript baseline_r/run_prediction.r
  ```
* **Convert the Formats**: Run `baseline_r/data_convertion.rmd` to convert the output `.rds` files into `.npy` format for evaluation.

### Foundation Models (Python)
* **Prepare the Data**: Run `foundation_python/data_prep.ipynb` to process the raw COVID-19 dataset.
* **Generate DeepAR Forecasts**: Navigate to `foundation_python/models/DeepAR/`. Preprocess training and forecasting dataset:

  ```bash
  python preprocess_elect.py
  python preprocess_covid.py
  ```
  Train the model:
  
  ```bash
  python train.py --sampling
  ```
  
  Generate sequential incidence forecast and adjust data structure:

  ```bash
  python evaluation.py
  python split_forecast.py
  ```

* **Generate Chronos-2 Forecasts**: Navigate to `foundation_python/models/Chronos2/`. Run the notebooks `chronos2_incidence_forecast.ipynb` and `chronos2_peak_time_forecast.ipynb` to generate standard incidence forecasts and engineered near-peak zero-shot forecasts:

* **Convert the Formats**: Run `baseline_r/data_convertion.rmd` to convert the output `.npy` files into `.rds` format for evaluation.

### Evaluations
* **Evaluate General Probabilistic Metrics**: Run `baseline_r/evaluation.rmd` to conduct general probabilistic evaluation.
* **Transform the Targets**: Run `foundation_python/generate_peaks.ipynb` to convert the incidence forecasts into peak-probability forecasts.
* **Evaluate Task-specific Performance**: Run `foundation_python/evaluation.ipynb` to compute the Cost-Loss elementary scores, construct Murphy diagrams, and evaluate the Brier scores.