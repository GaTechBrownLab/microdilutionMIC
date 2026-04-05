# MIC calculation for evolved lineages (M1-3, T1-3) and no antibiotic treatment (NA1)
# Author: Canan Karakoc
# Each plate has: Design CSV, T0h, T16h, T24h OD600 readings
# MIC calculated at both 16h and 24h as double-check

library(readxl)
library(tidyverse)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Change to yours!
project_dir <- "/Users/canankarakoc/Library/CloudStorage/Dropbox-GaTech/Canan Karakoc/CK_PS_LabData/02_Neisseria_ExpEvo/MIC_calculations"

experiment_id <- "Exp1"  # Change for each repeat (e.g., "Exp2", "Exp3")

# Input data lives under Data_Maps/<experiment_id>/Plate1, Plate2, ...
data_dir <- file.path(project_dir, "Data_Maps", experiment_id)

# Auto-detect plates: finds all Design_Plate*.csv files in the experiment folder
design_files <- sort(Sys.glob(file.path(data_dir, "Plate*", "Design_Plate*.csv")))

plates <- lapply(design_files, function(design_path) {
  plate_dir <- dirname(design_path)
  plate_name <- basename(plate_dir)

  # Find OD files by timepoint pattern
  all_xlsx <- sort(list.files(plate_dir, pattern = "\\.xlsx$", full.names = TRUE))
  t0_file  <- all_xlsx[grepl("T0h",  all_xlsx)]
  t16_file <- all_xlsx[grepl("T16h", all_xlsx)]
  t24_file <- all_xlsx[grepl("T24h", all_xlsx)]

  # T0 may not exist for every plate
  if (length(t0_file) == 0) t0_file <- NA_character_

  list(
    name = plate_name,
    experiment = experiment_id,
    design = design_path,
    t0  = t0_file,
    t16 = t16_file,
    t24 = t24_file
  )
})

cat("Found", length(plates), "plates in", data_dir, "\n")
for (p in plates) cat("  ", p$name, "\n")

od_cutoff <- 0.1  # MIC = lowest concentration where blanked OD falls below this

out_dir <- file.path(project_dir, "Codes_Results", "Results", experiment_id)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# HELPER: read plate-format Excel to long format
# =============================================================================

plate_to_long <- function(file_path, value_name) {
  df <- suppressMessages(read_excel(file_path, col_names = FALSE))
  df <- df[, 1:min(13, ncol(df))]

  row_letters <- LETTERS[1:8]
  plate_data <- df[2:9, 2:13]
  colnames(plate_data) <- as.character(1:12)
  plate_data$row <- row_letters

  plate_data %>%
    pivot_longer(cols = -row, names_to = "col", values_to = value_name) %>%
    mutate(
      col = as.integer(col),
      !!value_name := as.numeric(.data[[value_name]]),
      Well = paste0(row, col)
    )
}

# =============================================================================
# HELPER: read design CSV
# =============================================================================

read_design <- function(file_path) {
  design <- read_csv(file_path, show_col_types = FALSE)
  # Harmonize column name if older files still use "AB"
  if ("AB" %in% names(design) & !"Antibiotic" %in% names(design)) {
    design <- rename(design, Antibiotic = AB)
  }
  # Add Lineage column if not present (default to "evolved")
  if (!"Lineage" %in% names(design)) {
    design$Lineage <- "evolved"
  }
  design
}

# =============================================================================
# PROCESS ONE PLATE: merge design + OD, calculate growth & inhibition
# =============================================================================

process_plate <- function(plate_info) {
  plate_name <- plate_info$name
  experiment <- plate_info$experiment
  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat("Processing:", experiment, "-", plate_name, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  # Read design
  design <- read_design(plate_info$design)

  # Check if T0 file exists
  has_t0 <- !is.na(plate_info$t0)
  if (!has_t0) {
    cat("  WARNING: No T0 file found for", plate_name,
        "- using blank correction only (no T0 subtraction)\n")
  }

  # Read OD data
  od_t16 <- plate_to_long(plate_info$t16, "OD_t16")
  od_t24 <- plate_to_long(plate_info$t24, "OD_t24")

  # Merge design with OD readings
  dat <- design %>%
    left_join(od_t16 %>% select(Well, OD_t16), by = "Well") %>%
    left_join(od_t24 %>% select(Well, OD_t24), by = "Well") %>%
    mutate(Plate = plate_name, Experiment = experiment)

  if (has_t0) {
    od_t0 <- plate_to_long(plate_info$t0, "OD_t0")
    dat <- dat %>% left_join(od_t0 %>% select(Well, OD_t0), by = "Well")
  } else {
    dat <- dat %>% mutate(OD_t0 = NA_real_)
  }

  # Blank correction (average of blank wells per timepoint)
  blanks <- dat %>% filter(Type == "blank")
  blank_t16 <- mean(blanks$OD_t16, na.rm = TRUE)
  blank_t24 <- mean(blanks$OD_t24, na.rm = TRUE)

  if (has_t0) {
    blank_t0 <- mean(blanks$OD_t0, na.rm = TRUE)
    cat("  Blank means: T0 =", round(blank_t0, 4),
        ", T16h =", round(blank_t16, 4),
        ", T24h =", round(blank_t24, 4), "\n")
    dat <- dat %>%
      mutate(
        OD_t0_blanked  = OD_t0  - blank_t0,
        OD_t16_blanked = OD_t16 - blank_t16,
        OD_t24_blanked = OD_t24 - blank_t24,
        # T0-corrected growth
        growth_16h = OD_t16_blanked - OD_t0_blanked,
        growth_24h = OD_t24_blanked - OD_t0_blanked
      )
  } else {
    cat("  Blank means: T16h =", round(blank_t16, 4),
        ", T24h =", round(blank_t24, 4), "\n")
    dat <- dat %>%
      mutate(
        OD_t0_blanked  = NA_real_,
        OD_t16_blanked = OD_t16 - blank_t16,
        OD_t24_blanked = OD_t24 - blank_t24,
        growth_16h = OD_t16_blanked,
        growth_24h = OD_t24_blanked
      )
  }

  # Classify well types
  dat <- dat %>%
    mutate(
      well_type = case_when(
        Type == "blank" ~ "blank",
        Type == "empty" ~ "empty",
        Type == "sample" & (is.na(Concentration) | Concentration == 0) ~ "growth_control",
        Type == "sample" & Concentration > 0 ~ "treatment",
        TRUE ~ "other"
      )
    )

  return(dat)
}

# =============================================================================
# PROCESS ALL PLATES
# =============================================================================

all_data <- map_dfr(plates, process_plate)

# =============================================================================
# MIC DETERMINATION
# =============================================================================

# Treatment wells only
treatment <- all_data %>% filter(well_type == "treatment")

# Assign replicate number within each sample x plate
# (rows with same Sample & Concentration on the same plate are replicates)
treatment <- treatment %>%
  group_by(Experiment, Plate, Sample, Antibiotic, Lineage, Concentration) %>%
  mutate(replicate = row_number()) %>%
  ungroup()

# MIC per replicate: lowest concentration where blanked OD < cutoff
calc_mic <- function(df, od_col, timepoint_label) {
  df %>%
    group_by(Experiment, Plate, Sample, Antibiotic, Lineage, replicate) %>%
    arrange(Concentration) %>%
    filter(.data[[od_col]] < od_cutoff) %>%
    slice_min(Concentration, n = 1) %>%
    ungroup() %>%
    transmute(Experiment, Plate, Sample, Antibiotic, Lineage, replicate,
              mic = Concentration,
              timepoint = timepoint_label)
}

mic_16h <- calc_mic(treatment, "OD_t16_blanked", "16h")
mic_24h <- calc_mic(treatment, "OD_t24_blanked", "24h")
mic_all <- bind_rows(mic_16h, mic_24h)

# Summary: average across technical replicates per sample
mic_summary <- mic_all %>%
  group_by(Experiment, Sample, Antibiotic, Lineage, timepoint) %>%
  summarize(
    n_replicates = n(),
    mic_min = min(mic),
    mic_max = max(mic),
    mic_median = median(mic),
    mic_geometric_mean = exp(mean(log(mic))),
    .groups = "drop"
  ) %>%
  arrange(Experiment, Antibiotic, Sample, timepoint)

# =============================================================================
# RELATIVE MIC (evolved / ancestor)
# =============================================================================

# Calculate relative MIC if both ancestor and evolved data are present
has_ancestor <- any(mic_summary$Lineage == "ancestor")
has_evolved  <- any(mic_summary$Lineage == "evolved")

if (has_ancestor & has_evolved) {
  # Each sample has its own ancestor clone — match by Sample
  ancestor_mic <- mic_summary %>%
    filter(Lineage == "ancestor") %>%
    select(Experiment, Sample, Antibiotic, timepoint,
           ancestor_mic = mic_geometric_mean)

  relative_mic <- mic_summary %>%
    filter(Lineage == "evolved") %>%
    left_join(ancestor_mic, by = c("Experiment", "Sample", "Antibiotic", "timepoint")) %>%
    mutate(relative_mic = mic_geometric_mean / ancestor_mic)

  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat("RELATIVE MIC (evolved / ancestor)\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  print(relative_mic %>%
    select(Experiment, Sample, Antibiotic, timepoint,
           mic_geometric_mean, ancestor_mic, relative_mic),
    n = Inf)
} else {
  relative_mic <- NULL
  cat("\nNote: No ancestor lineage found - skipping relative MIC calculation.\n")
  cat("Add a 'Lineage' column to design files with 'ancestor' and 'evolved' levels\n")
  cat("to enable relative MIC calculation.\n")
}

# =============================================================================
# PRINT RESULTS
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("MIC PER REPLICATE\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")
print(mic_all %>% arrange(Experiment, Antibiotic, Sample, timepoint, replicate), n = Inf)

cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("MIC SUMMARY (averaged across technical replicates)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")
print(mic_summary, n = Inf)

# =============================================================================
# VISUALIZATION
# =============================================================================

# 1. Dose-response curves per sample (24h) with OD cutoff line
p1 <- treatment %>%
  group_by(Sample, Antibiotic, Concentration) %>%
  summarize(
    mean_OD = mean(OD_t24_blanked, na.rm = TRUE),
    se_OD = sd(OD_t24_blanked, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = Concentration, y = mean_OD, color = Sample)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_OD - se_OD, ymax = mean_OD + se_OD), width = 0.1) +
  geom_hline(yintercept = od_cutoff, linetype = "dashed", color = "gray50") +
  scale_x_log10() +
  facet_wrap(~Antibiotic, scales = "free_x") +
  labs(title = "Dose-Response (24h, blanked OD)",
       x = "Concentration (µg/mL)", y = "OD600",
       caption = paste("Dashed line = OD cutoff", od_cutoff)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

# 2. Dose-response curves (16h) for comparison
p2 <- treatment %>%
  group_by(Sample, Antibiotic, Concentration) %>%
  summarize(
    mean_OD = mean(OD_t16_blanked, na.rm = TRUE),
    se_OD = sd(OD_t16_blanked, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = Concentration, y = mean_OD, color = Sample)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_OD - se_OD, ymax = mean_OD + se_OD), width = 0.1) +
  geom_hline(yintercept = od_cutoff, linetype = "dashed", color = "gray50") +
  scale_x_log10() +
  facet_wrap(~Antibiotic, scales = "free_x") +
  labs(title = "Dose-Response (16h, blanked OD)",
       x = "Concentration (µg/mL)", y = "OD600",
       caption = paste("Dashed line = OD cutoff", od_cutoff)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

# 3. MIC comparison bar plot (16h vs 24h side by side)
p3 <- mic_summary %>%
  ggplot(aes(x = Sample, y = mic_geometric_mean, fill = timepoint)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(
    aes(ymin = mic_min, ymax = mic_max),
    position = position_dodge(width = 0.7), width = 0.2
  ) +
  scale_y_log10() +
  facet_wrap(~Antibiotic, scales = "free") +
  labs(title = "MIC Comparison (geometric mean ± range)",
       x = "Sample", y = "MIC (µg/mL)", fill = "Timepoint") +
  theme_bw() +
  theme(legend.position = "bottom")

# 4. Plate heatmaps (24h, blanked)
p4 <- all_data %>%
  filter(Type == "sample") %>%
  ggplot(aes(x = factor(as.integer(gsub("[A-H]", "", Well))),
             y = fct_rev(gsub("[0-9]+", "", Well)),
             fill = OD_t24_blanked)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(OD_t24_blanked, 2)), size = 1.8) +
  scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue",
                       midpoint = 0.2, name = "OD600") +
  facet_wrap(~Plate, ncol = 2) +
  labs(title = "Plate Heatmaps (24h, blanked)", x = "Column", y = "Row") +
  theme_minimal() +
  coord_equal()

# Save plots
ggsave(file.path(out_dir, "MIC_dose_response_24h.pdf"), p1, width = 10, height = 5)
ggsave(file.path(out_dir, "MIC_dose_response_16h.pdf"), p2, width = 10, height = 5)
ggsave(file.path(out_dir, "MIC_comparison_16h_vs_24h.pdf"), p3, width = 8, height = 5)
ggsave(file.path(out_dir, "MIC_plate_heatmaps.pdf"), p4, width = 12, height = 10)

# =============================================================================
# EXPORT
# =============================================================================

write_csv(all_data, file.path(out_dir, "MIC_all_well_data.csv"))
write_csv(mic_all, file.path(out_dir, "MIC_per_replicate.csv"))
write_csv(mic_summary, file.path(out_dir, "MIC_summary.csv"))
if (!is.null(relative_mic)) {
  write_csv(relative_mic, file.path(out_dir, "MIC_relative.csv"))
}

