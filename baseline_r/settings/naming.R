# models
models <- c("epiestim", "epinow2", "epidemia", "arima", "prophet", "gp")
models_nonqt <-c("epiestim", "epinow2", "epidemia", "arima", "prophet", "gp", "deepar")
models_all <- c("epiestim", "epinow2", "epidemia", "arima", "prophet", "gp", "deepar", "chronos2")
quantile_models <- c("chronos2")

model_names <- c("EpiEstim", "EpiNow2", "Epidemia", "SARIMA", "Prophet", "GP")
model_names_nonqt <- c("EpiEstim", "EpiNow2", "Epidemia", "SARIMA", "Prophet", "GP", "DeepAR")
model_names_all <- c("EpiEstim", "EpiNow2", "Epidemia", "SARIMA", "Prophet", "GP", "DeepAR","Chronos2")
quantile_model_names <- c("Chronos2")

names(model_names) <- models
names(model_names_nonqt) <- models_nonqt
names(model_names_all) <- models_all
names(quantile_model_names) <- quantile_models

alpha_levels <- c(0.05, 0.1, 0.15, 0.20, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)

# states
states <- c("az","ca","il","md","nj","ny")
state_names <- c("Arizona", "California", "Illinois", 
                 "Maryland", "New Jersey", "New York")
names(state_names) <- states

# phases
phase_names <- c("exponential growth", "subexponential growth", "plateau", 
                 "subexponential decline", "exponential decline")
phase_names_abrv <- c("Exponential\ngrowth", "Subexponential\ngrowth", "Plateau", 
                      "Subexponential\ndecline", "Exponential\ndecline")
phase_names_abrv2 <- c("Exponential growth", "Subexponential growth", "Plateau", 
                       "Subexponential decline", "Exponential decline")
names(phase_names_abrv) <- phase_names
names(phase_names_abrv2) <- phase_names