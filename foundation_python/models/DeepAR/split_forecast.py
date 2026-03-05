import numpy as np
import pandas as pd

import os
import matplotlib
import matplotlib.pyplot as plt

num_series   = 6
window_size  = 70
horizon      = 14
slide_size   = 7
time_len     = len(dates)
states = ["az", "ca", "il", "md", "nj", "ny"]
model_name = "base_model"


forecasts_wht = np.load(f"experiments/{model_name}/forecast_sample_sequential.npy")            # shape: (N_windows_total, 14, 200) already
incidence = pd.read_csv("../../data/2_processed/covid_incidence.csv", index_col= 0, parse_dates=True)
dates = incidence.index       # length = time_len

windows_per_series = (time_len - (window_size - slide_size)) // slide_size
assert windows_per_series * num_series == forecasts_wht.shape[0], (
    f"Expected {windows_per_series*num_series} windows, got {forecasts_wht.shape[0]}"
)

offset = window_size - horizon   # = 56 (context length)

for s_idx, state in enumerate(states):
    start = s_idx * windows_per_series
    end   = start + windows_per_series

    F_state = forecasts_wht[start:end, :, :]   # shape: (56, 14, 200) for your 458-day example
    np.save(f"../../data/3_incidence predictions/deepar/{state}_forecasts.npy", F_state)
    print(state, F_state.shape)
