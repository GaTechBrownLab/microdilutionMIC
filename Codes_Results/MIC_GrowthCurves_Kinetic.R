# Kinetic MIC + growth-curve analysis for evolved/ancestor lineages
# Author: Canan Karakoc
# Each plate: one continuous-read OD600 CSV (Time x 96 wells) + Design CSV
# Outputs:
#   - MIC at user-chosen sampling timepoints (e.g. 16h, 24h)
#   - Growth-curve plots (all wells; growth-control wells overlaid)
#   - Life-history traits (lag, mu_max, K, AUC) fit on zero-antibiotic wells
#   - Relative metrics (evolved / ancestor) for both MIC and life-history traits

# ---- Dependencies -----------------------------------------------------------
# install.packages("gcplyr")  # one-time
suppressPackageStartupMessages({
  library(tidyverse)
  library(gcplyr)
  library(cowplot)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

project_dir <- "/Users/canankarakoc/Documents/GitHub/microdilutionMIC"

experiment_id <- "Exp_Kinetic"

# Timepoints (in hours) at which to snapshot OD for MIC determination.
# Pipeline picks the kinetic read closest to each requested hour.
mic_timepoints_h <- c(16, 24, 48)

# OD cutoff used to call MIC (matches MIC_Evolved_Lineages.R)
od_cutoff <- 0.05

# Optional truncation: drop timepoints beyond this many hours before fitting.
# Set to NA to keep the full curve. Inspect the all-wells plot, then decide.
truncate_hours <- NA

# Smoothing window (number of timepoints) for spline-based growth metrics.
# Plate1 reads ~every 30 min, so window = 5 ≈ 2.5 h moving average.
smooth_window_n <- 5

# Window used to smooth the per-capita derivative itself.
deriv_window_n <- 7

# Minimum smoothed OD considered "real growth" — derivative is masked below
# this so noise spikes near zero (OD floor / blank artefacts) don't inflate
# mu_max or distort lag_time.
mu_min_OD <- 0.02

# Parametric growth model: "gompertz" (asymmetric, default) or "logistic"
# (symmetric). Both parameterized in K, mu_max, lag (Zwietering 1990).
growth_model <- "gompertz"

# Trim each well at its peak OD (+ buffer) before fitting. Prevents late
# post-stationary decline / contamination spikes from biasing K downward.
trim_post_peak     <- TRUE
post_peak_buffer_h <- 4

# Hard cap on fit window (hours from t = 0). Set to NA to disable. Useful
# when curves drift post-stationary and you want a uniform cutoff across
# wells (try e.g. 36 / 48 / 60 if K still looks off after trimming).
fit_truncate_hours <- NA

# Concentration treated as "growth control" for life-history fits
control_conc <- 0

data_dir <- file.path(project_dir, "Data_Maps", experiment_id)
out_dir  <- file.path(project_dir, "Codes_Results", "Results", experiment_id)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Auto-detect plates
design_files <- sort(Sys.glob(file.path(data_dir, "Plate*", "Design_Plate*.csv")))

plates <- lapply(design_files, function(design_path) {
  plate_dir  <- dirname(design_path)
  plate_name <- basename(plate_dir)
  csv_files  <- list.files(plate_dir, pattern = "\\.csv$", full.names = TRUE)
  od_file    <- csv_files[!grepl("Design_", basename(csv_files))]
  if (length(od_file) != 1)
    stop("Expected exactly one OD CSV in ", plate_dir, " (found ", length(od_file), ")")
  list(name = plate_name, design = design_path, od = od_file)
})

cat("Found", length(plates), "plates in", data_dir, "\n")
for (p in plates) cat("  ", p$name, "\n")

# =============================================================================
# HELPERS
# =============================================================================

# Convert "HH:MM:SS" to numeric hours
parse_time_hours <- function(x) {
  parts <- strsplit(as.character(x), ":", fixed = TRUE)
  vapply(parts, function(p) {
    p <- as.numeric(p)
    p[1] + p[2] / 60 + p[3] / 3600
  }, numeric(1))
}

read_kinetic <- function(file_path) {
  # Read everything as character so trailing reader-summary rows
  # (e.g. "Max V", "Lagtime") don't break column types.
  df <- read_csv(file_path, show_col_types = FALSE,
                 col_types = cols(.default = "c"))
  names(df)[1] <- "Time"

  # Keep only true kinetic rows: Time matches HH:MM:SS
  df <- df %>% filter(grepl("^\\d+:\\d+:\\d+$", Time))

  df %>%
    mutate(Time_h = parse_time_hours(Time)) %>%
    select(Time_h, everything(), -Time) %>%
    pivot_longer(cols = -Time_h, names_to = "Well", values_to = "OD") %>%
    mutate(OD = suppressWarnings(as.numeric(OD)))
}

read_design <- function(file_path) {
  d <- read_csv(file_path, show_col_types = FALSE)
  if ("AB" %in% names(d) && !"Antibiotic" %in% names(d)) d <- rename(d, Antibiotic = AB)
  if (!"Lineage" %in% names(d)) d$Lineage <- "evolved"
  d
}

# =============================================================================
# LOAD + BLANK-CORRECT EACH PLATE
# =============================================================================

process_plate <- function(plate_info) {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("Processing:", plate_info$name, "\n")

  design <- read_design(plate_info$design)
  kin    <- read_kinetic(plate_info$od)

  if (!is.na(truncate_hours)) {
    kin <- kin %>% filter(Time_h <= truncate_hours)
    cat("  Truncated to", truncate_hours, "h\n")
  }

  # Per-timepoint blank mean (handles drift better than a global blank)
  blank_wells <- design %>% filter(Type == "blank") %>% pull(Well)
  blank_per_t <- kin %>%
    filter(Well %in% blank_wells) %>%
    group_by(Time_h) %>%
    summarize(blank_mean = mean(OD, na.rm = TRUE), .groups = "drop")

  joined <- kin %>%
    left_join(blank_per_t, by = "Time_h") %>%
    left_join(design,      by = "Well") %>%
    mutate(
      OD_blanked = OD - blank_mean,
      Plate      = plate_info$name
    )

  cat("  Wells:", n_distinct(joined$Well),
      " | Reads:", n_distinct(joined$Time_h),
      " | Span:", round(max(joined$Time_h, na.rm = TRUE), 1), "h\n")

  joined
}

all_data <- map_dfr(plates, process_plate)

# =============================================================================
# REPLICATE NUMBERING (within plate x sample x lineage x concentration)
# =============================================================================

# Number replicate wells once (consistent across timepoints) using the order
# of distinct wells per group.
well_replicates <- all_data %>%
  filter(Type == "sample") %>%
  distinct(Plate, Sample, Antibiotic, Lineage, Concentration, Well) %>%
  group_by(Plate, Sample, Antibiotic, Lineage, Concentration) %>%
  arrange(Well, .by_group = TRUE) %>%
  mutate(replicate = row_number()) %>%
  ungroup()

all_data <- all_data %>%
  left_join(well_replicates,
            by = c("Plate", "Sample", "Antibiotic", "Lineage",
                   "Concentration", "Well"))

# =============================================================================
# 1. MIC AT SAMPLED TIMEPOINTS
# =============================================================================

# For each requested timepoint, pick the kinetic read closest to it
# **per plate** — different plates have slightly different read schedules,
# so a single global "nearest" timepoint matches only one plate and silently
# drops data from the others.
snap_at <- function(df, target_h) {
  df %>%
    group_by(Plate) %>%
    group_modify(~ {
      reads   <- sort(unique(.x$Time_h))
      nearest <- reads[which.min(abs(reads - target_h))]
      .x %>%
        filter(Time_h == nearest) %>%
        mutate(timepoint_label  = paste0(round(target_h), "h"),
               timepoint_actual = nearest)
    }) %>%
    ungroup()
}

snapshots <- map_dfr(mic_timepoints_h, ~ snap_at(all_data, .x))

# Diagnostic: max OD reached in growth-control wells at each MIC timepoint.
# If this is barely above the cutoff, the cutoff is too high and MICs read low.
control_snap_check <- snapshots %>%
  filter(Type == "sample", Concentration == control_conc) %>%
  group_by(Lineage, Sample, Antibiotic, timepoint_label) %>%
  summarize(
    n_wells          = n(),
    control_OD_mean  = mean(OD_blanked, na.rm = TRUE),
    control_OD_max   = max(OD_blanked,  na.rm = TRUE),
    .groups = "drop"
  )

cat("\n", strrep("=", 60), "\n", sep = "")
cat("CUTOFF DIAGNOSTIC: max OD in growth controls vs cutoff =", od_cutoff, "\n",
    strrep("=", 60), "\n\n", sep = "")
print(control_snap_check, n = Inf)

# MIC = lowest concentration at which blanked OD < cutoff (per replicate)
treatment_snap <- snapshots %>%
  filter(Type == "sample", Concentration > 0)

# Keep every replicate. If no concentration drops below the cutoff,
# record mic = NA and flag censored = TRUE (resistant beyond max tested conc).
mic_per_rep <- treatment_snap %>%
  group_by(Plate, Sample, Antibiotic, Lineage, replicate, timepoint_label) %>%
  arrange(Concentration, .by_group = TRUE) %>%
  summarize(
    max_conc_tested = max(Concentration, na.rm = TRUE),
    mic = {
      below <- which(OD_blanked < od_cutoff)
      if (length(below) > 0) Concentration[min(below)] else NA_real_
    },
    censored = is.na(mic),
    .groups = "drop"
  ) %>%
  rename(timepoint = timepoint_label)

# Diagnostic: how many replicates had no MIC reached at each timepoint?
censoring <- mic_per_rep %>%
  group_by(Sample, Antibiotic, Lineage, timepoint) %>%
  summarize(
    n_total    = n(),
    n_censored = sum(censored),
    .groups = "drop"
  )
cat("\nMIC censoring (replicates where no conc reached OD < ", od_cutoff, "):\n", sep = "")
print(censoring, n = Inf)

mic_summary <- mic_per_rep %>%
  group_by(Sample, Antibiotic, Lineage, timepoint) %>%
  summarize(
    n_replicates       = n(),
    n_censored         = sum(censored),
    max_conc_tested    = max(max_conc_tested, na.rm = TRUE),
    mic_min            = suppressWarnings(min(mic,    na.rm = TRUE)),
    mic_max            = suppressWarnings(max(mic,    na.rm = TRUE)),
    mic_median         = median(mic, na.rm = TRUE),
    mic_mean           = mean(mic,   na.rm = TRUE),
    mic_se             = sd(mic, na.rm = TRUE) / sqrt(sum(!is.na(mic))),
    mic_geometric_mean = exp(mean(log(mic), na.rm = TRUE)),
    mic_geo_se_lower   = exp(mean(log(mic), na.rm = TRUE) -
                             sd(log(mic), na.rm = TRUE) / sqrt(sum(!is.na(mic)))),
    mic_geo_se_upper   = exp(mean(log(mic), na.rm = TRUE) +
                             sd(log(mic), na.rm = TRUE) / sqrt(sum(!is.na(mic)))),
    .groups = "drop"
  ) %>%
  # If every replicate was censored, the geometric mean is NaN — replace with
  # ">max_conc_tested" by setting the value just above max for plotting.
  mutate(across(c(mic_geometric_mean, mic_geo_se_lower, mic_geo_se_upper,
                  mic_mean, mic_min, mic_max, mic_median),
                ~ ifelse(is.finite(.x), .x, NA_real_))) %>%
  arrange(Antibiotic, Sample, timepoint)

# Relative MIC (evolved / ancestor), matched by Sample + Antibiotic + timepoint
relative_mic <- NULL
if (all(c("ancestor", "evolved") %in% mic_summary$Lineage)) {
  ancestor_mic <- mic_summary %>%
    filter(Lineage == "ancestor") %>%
    select(Sample, Antibiotic, timepoint,
           ancestor_mic = mic_geometric_mean)

  relative_mic <- mic_summary %>%
    filter(Lineage == "evolved") %>%
    left_join(ancestor_mic, by = c("Sample", "Antibiotic", "timepoint")) %>%
    mutate(relative_mic = mic_geometric_mean / ancestor_mic)
}

# =============================================================================
# 2. GROWTH CURVES — PLOTS
# =============================================================================

sample_data <- all_data %>% filter(Type == "sample")

# 2a. All wells, faceted by Lineage x Sample, colored by Concentration
p_all <- sample_data %>%
  ggplot(aes(x = Time_h, y = OD_blanked,
             group = Well,
             color = log10(pmax(Concentration, 1e-3)))) +
  geom_line(alpha = 0.7, linewidth = 0.4) +
  scale_color_viridis_c(name = "log10[conc]\n(0 → -3)") +
  facet_grid(Lineage ~ Sample + Antibiotic) +
  labs(title = "All wells: kinetic OD600 (blanked)",
       x = "Time (h)", y = "OD600 (blanked)") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 8))

# 2b. Growth controls only (Concentration == 0): per-replicate curves
p_controls <- sample_data %>%
  filter(Concentration == control_conc) %>%
  ggplot(aes(x = Time_h, y = OD_blanked,
             color = Lineage, group = interaction(Well, Plate))) +
  geom_line(alpha = 0.8, linewidth = 0.6) +
  facet_wrap(~Sample + Antibiotic) +
  labs(title = paste0("Growth controls (Concentration = ", control_conc, "): per-well curves"),
       x = "Time (h)", y = "OD600 (blanked)") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "growth_curves_all_wells.pdf"),
       p_all, width = 14, height = 8)
ggsave(file.path(out_dir, "growth_curves_controls.pdf"),
       p_controls, width = 10, height = 6)

# 2c. Per-plate "as-laid-out" view: each well in its row x col position with
# OD cutoff line and MIC timepoint markers. Lets you visually verify whether
# the MIC call at each timepoint matches what the growth curve actually does.
plot_plate_layout <- function(plate_data, plate_name) {
  pd <- plate_data %>%
    mutate(
      row_letter = factor(gsub("[0-9]+", "", Well), levels = LETTERS[1:8]),
      col_num    = factor(as.integer(gsub("[A-H]", "", Well)),
                          levels = 1:12)
    )
  annot <- pd %>%
    distinct(Well, row_letter, col_num, Sample, Concentration, Lineage, Type) %>%
    mutate(
      label = case_when(
        Type == "blank"  ~ "blank",
        Type == "empty"  ~ "empty",
        Type == "sample" & Concentration == 0 ~
          paste0(Lineage, "\nctrl"),
        Type == "sample" ~
          paste0(Lineage, "\n", Concentration),
        TRUE ~ ""
      ),
      well_color = case_when(
        Type == "blank"  ~ "blank",
        Type == "empty"  ~ "empty",
        Type == "sample" & Concentration == 0 ~ "control",
        Type == "sample" ~ as.character(Lineage),
        TRUE ~ "other"
      )
    )
  pd <- pd %>% left_join(annot %>% select(Well, well_color), by = "Well")

  ggplot(pd, aes(x = Time_h, y = OD_blanked, group = Well)) +
    geom_vline(xintercept = mic_timepoints_h, linetype = "dotted",
               color = "gray50", linewidth = 0.25) +
    geom_hline(yintercept = od_cutoff, linetype = "dashed",
               color = "firebrick", linewidth = 0.3) +
    geom_line(aes(color = well_color), linewidth = 0.45) +
    geom_text(data = annot, aes(label = label),
              x = Inf, y = Inf, hjust = 1.02, vjust = 1.15,
              size = 1.7, color = "gray20", inherit.aes = FALSE) +
    scale_color_manual(values = c(
      ancestor = "#1f77b4", evolved = "#d62728",
      control  = "darkgreen",
      blank    = "gray70",  empty = "gray85", other = "black"
    ), name = NULL) +
    facet_grid(row_letter ~ col_num, switch = "y") +
    labs(title = paste("Plate layout — kinetic curves:", plate_name),
         subtitle = paste("Dashed red = od_cutoff =", od_cutoff,
                          "  |  dotted gray = MIC timepoints (",
                          paste0(mic_timepoints_h, "h", collapse = ", "), ")"),
         x = "Time (h)", y = "OD600 (blanked)") +
    theme_bw() +
    theme(strip.text   = element_text(size = 7),
          axis.text    = element_text(size = 5),
          panel.spacing = unit(0.15, "lines"),
          legend.position = "bottom")
}

for (plate_name in unique(all_data$Plate)) {
  pd <- all_data %>% filter(Plate == plate_name)
  p_layout <- plot_plate_layout(pd, plate_name)
  ggsave(file.path(out_dir, paste0("plate_layout_", plate_name, ".pdf")),
         p_layout, width = 16, height = 10)
}

# 2c. Dose-response at each MIC timepoint with cutoff line — diagnostic for od_cutoff
dose_resp <- snapshots %>%
  filter(Type == "sample") %>%
  group_by(Sample, Antibiotic, Lineage, timepoint_label, Concentration) %>%
  summarize(
    mean_OD = mean(OD_blanked, na.rm = TRUE),
    se_OD   = sd(OD_blanked,   na.rm = TRUE) / sqrt(sum(!is.na(OD_blanked))),
    .groups = "drop"
  ) %>%
  mutate(timepoint_label = factor(timepoint_label,
                                  levels = paste0(round(sort(mic_timepoints_h)), "h")))

p_dr <- ggplot(dose_resp,
               aes(x = pmax(Concentration, min(Concentration[Concentration > 0]) / 3),
                   y = mean_OD, color = Lineage)) +
  geom_line(alpha = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_OD - se_OD, ymax = mean_OD + se_OD),
                width = 0.08) +
  geom_hline(yintercept = od_cutoff, linetype = "dashed", color = "gray40") +
  scale_x_log10() +
  facet_grid(Antibiotic + Sample ~ timepoint_label, scales = "free_y") +
  labs(title = "Dose-response at MIC timepoints (mean ± SE)",
       subtitle = paste("Dashed line = od_cutoff =", od_cutoff,
                        "  |  zero-conc plotted at 1/3 of lowest non-zero"),
       x = "Concentration (µg/mL, log)", y = "OD600 (blanked)") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "dose_response_at_timepoints.pdf"),
       p_dr, width = 11, height = 8)

# =============================================================================
# 3. LIFE-HISTORY TRAITS — gcplyr on growth-control wells
# =============================================================================

controls <- sample_data %>%
  filter(Concentration == control_conc) %>%
  arrange(Plate, Well, Time_h)

# Smooth raw OD (used only for K starting value, AUC, and plotting)
controls <- controls %>%
  group_by(Plate, Well) %>%
  mutate(
    OD_smooth = smooth_data(
      x = Time_h, y = OD_blanked,
      sm_method = "moving-average",
      window_width_n = smooth_window_n
    )
  ) %>%
  ungroup()

# ----- Sigmoidal fits per well -----
# Both parameterized in (K, mu, lag) per Zwietering 1990:
#   Gompertz: y(t) = K * exp(-exp(mu*e/K * (lag - t) + 1))
#   Logistic: y(t) = K / (1 + exp(4*mu/K * (lag - t) + 2))
gompertz_fn <- function(t, K, mu, lag) {
  K * exp(-exp((mu * exp(1) / K) * (lag - t) + 1))
}
logistic_fn <- function(t, K, mu, lag) {
  K / (1 + exp(4 * mu / K * (lag - t) + 2))
}

fit_growth_well <- function(times, od, model = "gompertz",
                            trim = TRUE, buffer_h = 4,
                            hard_cap_h = NA_real_) {
  fn <- switch(model, gompertz = gompertz_fn, logistic = logistic_fn,
               stop("Unknown growth_model: ", model))

  ok <- is.finite(times) & is.finite(od)
  t <- times[ok]; y <- od[ok]
  na_out <- list(K = NA_real_, mu = NA_real_, lag = NA_real_,
                 pred = rep(NA_real_, length(times)),
                 rmse = NA_real_, converged = FALSE)
  if (length(t) < 10) return(na_out)

  # Hard cap (uniform across wells) overrides peak-trim if specified
  if (!is.na(hard_cap_h)) {
    keep <- t <= hard_cap_h
    if (sum(keep) >= 10) { t <- t[keep]; y <- y[keep] }
  } else if (trim) {
    # Trim post-peak: drop points beyond peak smoothed OD + buffer_h
    half <- 3L
    n_y  <- length(y)
    y_sm <- vapply(seq_len(n_y), function(i) {
      lo <- max(1L, i - half); hi <- min(n_y, i + half)
      mean(y[lo:hi], na.rm = TRUE)
    }, numeric(1))
    i_peak   <- which.max(y_sm)
    t_cutoff <- t[i_peak] + buffer_h
    keep     <- t <= t_cutoff
    if (sum(keep) >= 10) { t <- t[keep]; y <- y[keep] }
  }

  # Smoothed y on the (trimmed) data; K_start = max smoothed OD.
  # This is the highest plateau actually observed — avoids the optimizer
  # pushing K above the real asymptote.
  half <- 3L
  n_y  <- length(y)
  y_sm <- vapply(seq_len(n_y), function(i) {
    lo <- max(1L, i - half); hi <- min(n_y, i + half)
    mean(y[lo:hi], na.rm = TRUE)
  }, numeric(1))
  K_start <- max(y_sm, na.rm = TRUE)
  if (!is.finite(K_start) || K_start <= 0.01) return(na_out)

  # Robust mu start: median of the top derivatives in the actual growth phase
  log_y <- log(pmax(y, 1e-4))
  dy_dt <- diff(log_y) / diff(t)
  in_growth <- y[-1] > 0.05 * K_start
  if (any(in_growth)) {
    top <- sort(dy_dt[in_growth], decreasing = TRUE)
    mu_start <- median(head(top, max(3, floor(length(top) / 4))), na.rm = TRUE)
  } else {
    mu_start <- 0.3
  }
  if (!is.finite(mu_start) || mu_start <= 0) mu_start <- 0.3

  # Multi-start lag candidates (in hours)
  first_above <- function(frac) {
    i <- which(y >= frac * K_start)[1]
    if (is.na(i)) NA_real_ else t[i]
  }
  lag_candidates <- unique(pmax(c(
    0,
    first_above(0.05),
    first_above(0.10),
    first_above(0.25),
    max(t) / 8,
    max(t) / 4
  ), 0))
  lag_candidates <- lag_candidates[is.finite(lag_candidates)]

  best <- NULL
  for (lag_start in lag_candidates) {
    fit <- try(
      nls(y ~ fn(t, K, mu, lag),
          start = list(K = K_start, mu = mu_start, lag = lag_start),
          lower = list(K = 1e-3,            mu = 1e-4, lag = 0),
          upper = list(K = 1.05 * K_start,  mu = 5,    lag = max(t)),
          algorithm = "port",
          control = nls.control(maxiter = 300, warnOnly = TRUE)),
      silent = TRUE
    )
    if (!inherits(fit, "try-error")) {
      ssr <- sum(residuals(fit)^2)
      if (is.null(best) || ssr < best$ssr) best <- list(fit = fit, ssr = ssr)
    }
  }

  if (is.null(best)) return(na_out)

  cf <- coef(best$fit)
  # Predictions over all original times (for plotting overlay)
  pred <- fn(times, cf["K"], cf["mu"], cf["lag"])
  # RMSE on the **fit window** only — pred at the trimmed times vs trimmed y
  pred_fit <- fn(t, cf["K"], cf["mu"], cf["lag"])
  rmse <- sqrt(mean((pred_fit - y)^2, na.rm = TRUE))
  list(K = unname(cf["K"]), mu = unname(cf["mu"]), lag = unname(cf["lag"]),
       pred = pred, rmse = rmse, converged = TRUE)
}

# Apply per well — keep both per-timepoint predictions and per-well params
gompertz_fits <- controls %>%
  group_by(Plate, Sample, Antibiotic, Lineage, replicate, Well) %>%
  group_modify(~ {
    f <- fit_growth_well(.x$Time_h, .x$OD_blanked,
                         model = growth_model,
                         trim  = trim_post_peak,
                         buffer_h   = post_peak_buffer_h,
                         hard_cap_h = fit_truncate_hours)
    tibble(
      Time_h    = .x$Time_h,
      OD_pred   = f$pred,
      K_fit     = f$K,
      mu_fit    = f$mu,
      lag_fit   = f$lag,
      rmse_fit  = f$rmse,
      converged = f$converged
    )
  }) %>%
  ungroup()

# Per-well metrics from the parametric fit; AUC stays numerical (model-free)
life_history <- gompertz_fits %>%
  group_by(Plate, Sample, Antibiotic, Lineage, replicate, Well) %>%
  summarize(
    K_max_OD  = first(K_fit),
    mu_max    = first(mu_fit),
    lag_time  = first(lag_fit),
    rmse      = first(rmse_fit),
    converged = first(converged),
    .groups = "drop"
  ) %>%
  left_join(
    controls %>%
      group_by(Plate, Sample, Antibiotic, Lineage, replicate, Well) %>%
      summarize(auc = auc(x = Time_h, y = OD_smooth), .groups = "drop"),
    by = c("Plate","Sample","Antibiotic","Lineage","replicate","Well")
  )

cat("\n", strrep("=", 60), "\n", sep = "")
cat("GROWTH-MODEL FIT (", growth_model, ") — convergence + RMSE\n",
    strrep("=", 60), "\n", sep = "")
print(life_history %>%
        group_by(Sample, Lineage, converged) %>%
        summarize(n = n(),
                  rmse_med = median(rmse, na.rm = TRUE),
                  rmse_max = max(rmse,    na.rm = TRUE),
                  .groups = "drop") %>%
        arrange(Sample, Lineage, desc(converged)),
      n = Inf)
cat("\nWorst 10 fits by RMSE (potentially struggling):\n")
print(life_history %>%
        filter(converged) %>%
        arrange(desc(rmse)) %>%
        select(Plate, Sample, Lineage, replicate, Well,
               K_max_OD, mu_max, lag_time, rmse) %>%
        head(10))

# Summary across replicates
life_history_summary <- life_history %>%
  group_by(Sample, Antibiotic, Lineage) %>%
  summarize(
    n_replicates = n(),
    across(c(K_max_OD, mu_max, lag_time, auc),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd   = ~sd(.x,   na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

# Relative life-history traits (evolved / ancestor)
relative_lh <- NULL
if (all(c("ancestor", "evolved") %in% life_history_summary$Lineage)) {
  anc_lh <- life_history_summary %>%
    filter(Lineage == "ancestor") %>%
    select(Sample, Antibiotic,
           K_anc   = K_max_OD_mean,
           mu_anc  = mu_max_mean,
           lag_anc = lag_time_mean,
           auc_anc = auc_mean)

  relative_lh <- life_history_summary %>%
    filter(Lineage == "evolved") %>%
    left_join(anc_lh, by = c("Sample", "Antibiotic")) %>%
    mutate(
      rel_K   = K_max_OD_mean / K_anc,
      rel_mu  = mu_max_mean   / mu_anc,
      rel_lag = lag_time_mean / lag_anc,
      rel_auc = auc_mean      / auc_anc
    ) %>%
    select(Sample, Antibiotic, n_replicates,
           K_max_OD_mean, K_anc, rel_K,
           mu_max_mean,   mu_anc, rel_mu,
           lag_time_mean, lag_anc, rel_lag,
           auc_mean,      auc_anc, rel_auc)
}

# =============================================================================
# 4. PLOTS — life-history & MIC summaries
# =============================================================================

# Life-history: mean ± SE per Sample x Antibiotic x Lineage
lh_plot_data <- life_history %>%
  pivot_longer(c(K_max_OD, mu_max, lag_time, auc),
               names_to = "trait", values_to = "value") %>%
  group_by(Sample, Antibiotic, Lineage, trait) %>%
  summarize(
    mean_val = mean(value, na.rm = TRUE),
    se_val   = sd(value,   na.rm = TRUE) / sqrt(sum(!is.na(value))),
    .groups  = "drop"
  )

p_lh <- ggplot(lh_plot_data,
               aes(x = Lineage, y = mean_val, color = Lineage)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                width = 0.15) +
  facet_grid(trait ~ Sample + Antibiotic, scales = "free_y") +
  labs(title = "Life-history traits (growth-control wells; mean ± SE)",
       x = NULL, y = "Value") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.y = element_text(size = 8))

# MIC: timepoints across columns (left-to-right in time), antibiotics down rows
mic_summary <- mic_summary %>%
  mutate(timepoint = factor(timepoint,
                            levels = paste0(round(sort(mic_timepoints_h)), "h")))

# Build a plotting frame that shows fully-resistant groups (every replicate
# censored) as ">max_conc_tested" with a triangle marker, instead of dropping.
mic_plot_data <- mic_summary %>%
  mutate(
    fully_censored = (n_censored == n_replicates),
    plot_y    = ifelse(fully_censored, max_conc_tested, mic_geometric_mean),
    plot_low  = ifelse(fully_censored, NA_real_,         mic_geo_se_lower),
    plot_high = ifelse(fully_censored, NA_real_,         mic_geo_se_upper)
  )

p_mic <- ggplot(mic_plot_data,
                aes(x = Sample, y = plot_y, color = Lineage)) +
  geom_point(aes(shape = fully_censored),
             position = position_dodge(width = 0.5), size = 2.5) +
  geom_errorbar(aes(ymin = plot_low, ymax = plot_high),
                position = position_dodge(width = 0.5), width = 0.2,
                na.rm = TRUE) +
  scale_y_log10() +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 2),
                     labels = c(`FALSE` = "MIC reached",
                                `TRUE`  = ">max conc (censored)"),
                     name = NULL) +
  facet_grid(Antibiotic ~ timepoint, scales = "free_x") +
  labs(title = "MIC at sampled timepoints (geometric mean ± SE)",
       x = "Sample", y = "MIC (µg/mL)") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "life_history_traits.pdf"),
       p_lh, width = 12, height = 8)
ggsave(file.path(out_dir, "MIC_kinetic_timepoints.pdf"),
       p_mic, width = 10, height = 6)

# Per-well fit overlay: raw OD points + Gompertz fit curve + K and lag markers.
fit_annot <- life_history %>%
  select(Plate, Sample, Lineage, replicate, Well,
         K_max_OD, mu_max, lag_time, converged)

raw_with_fit <- controls %>%
  left_join(gompertz_fits %>%
              select(Plate, Sample, Antibiotic, Lineage, replicate, Well,
                     Time_h, OD_pred),
            by = c("Plate","Sample","Antibiotic","Lineage","replicate","Well","Time_h"))

p_fits <- ggplot(raw_with_fit, aes(x = Time_h, group = Well)) +
  geom_point(aes(y = OD_blanked, color = Lineage),
             size = 0.4, alpha = 0.4) +
  geom_line(aes(y = OD_pred, color = Lineage),
            linewidth = 0.8, na.rm = TRUE) +
  geom_hline(data = fit_annot,
             aes(yintercept = K_max_OD),
             linetype = "dashed", color = "gray30", linewidth = 0.3) +
  geom_vline(data = fit_annot,
             aes(xintercept = lag_time),
             linetype = "dotted", color = "firebrick", linewidth = 0.4) +
  facet_wrap(~ Sample + Lineage + Well, ncol = 6) +
  labs(title = "Modified Gompertz fits per well (controls only)",
       subtitle = paste("Points = raw blanked OD  |  line = Gompertz fit",
                        "\nDashed gray = K (asymptote)  |  dotted red = lag (lambda)"),
       x = "Time (h)", y = "OD600 (blanked)") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 6))

ggsave(file.path(out_dir, "growth_fits_per_well.pdf"),
       p_fits, width = 14, height = 10)

# =============================================================================
# 5. PRINT + EXPORT
# =============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("MIC SUMMARY\n", strrep("=", 60), "\n\n", sep = "")
print(mic_summary, n = Inf)

if (!is.null(relative_mic)) {
  cat("\nRELATIVE MIC (evolved / ancestor)\n", strrep("-", 40), "\n", sep = "")
  print(relative_mic %>%
          select(Sample, Antibiotic, timepoint,
                 mic_geometric_mean, ancestor_mic, relative_mic),
        n = Inf)
}

cat("\n", strrep("=", 60), "\n", sep = "")
cat("LIFE-HISTORY TRAIT SUMMARY (growth controls)\n",
    strrep("=", 60), "\n\n", sep = "")
print(life_history_summary, n = Inf)

if (!is.null(relative_lh)) {
  cat("\nRELATIVE LIFE-HISTORY TRAITS (evolved / ancestor)\n",
      strrep("-", 40), "\n", sep = "")
  print(relative_lh, n = Inf)
}

# CSV exports
write_csv(all_data,             file.path(out_dir, "kinetic_all_wells_long.csv"))
write_csv(snapshots,            file.path(out_dir, "kinetic_snapshots.csv"))
write_csv(mic_per_rep,          file.path(out_dir, "MIC_per_replicate.csv"))
write_csv(mic_summary,          file.path(out_dir, "MIC_summary.csv"))
write_csv(life_history,         file.path(out_dir, "life_history_per_well.csv"))
write_csv(life_history_summary, file.path(out_dir, "life_history_summary.csv"))
if (!is.null(relative_mic)) write_csv(relative_mic, file.path(out_dir, "MIC_relative.csv"))
if (!is.null(relative_lh))  write_csv(relative_lh,  file.path(out_dir, "life_history_relative.csv"))

cat("\nDone. Outputs in:", out_dir, "\n")
