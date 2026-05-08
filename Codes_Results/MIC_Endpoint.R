# Endpoint MIC analysis (T0/T16/T24 single reads), all experiments merged
# Author: Canan Karakoc
#
# Auto-discovers every Data_Maps/Exp*/Plate*/Design_*.csv (excludes Exp_Kinetic),
# stacks the T0/T16/T24 OD plates into one long table, blank-corrects per
# (plate, timepoint), and computes MIC per replicate with censoring.
# Plot styles mirror MIC_GrowthCurves_Kinetic.R so endpoint and kinetic results
# can be compared directly.

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

project_dir <- "/Users/canankarakoc/Documents/GitHub/microdilutionMIC"

# Timepoints (hours) read from filenames. Must match the "T<N>h" tag in xlsx names.
mic_timepoints_h <- c(0, 16, 24)

# OD cutoff for the MIC call (matches kinetic script default; tune if needed)
od_cutoff <- 0.1

# Concentration treated as growth control
control_conc <- 0

# When a filename's antibiotic token disagrees with the design's Antibiotic
# column: "filename" trusts the filename (overrides design), "design" trusts
# the design (skips the file), "stop" aborts. Default is "filename" because
# experimentally the filename usually reflects what was actually plated.
antibiotic_conflict_policy <- "filename"

data_root <- file.path(project_dir, "Data_Maps")
out_dir   <- file.path(project_dir, "Codes_Results", "Results", "Endpoint")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Discover plates: every Design_*.csv under Data_Maps/Exp*/Plate*/, except Exp_Kinetic.
design_files <- Sys.glob(file.path(data_root, "Exp*", "Plate*", "Design*.csv"))
design_files <- design_files[!grepl("Exp_Kinetic", design_files, fixed = TRUE)]
design_files <- sort(design_files)

if (length(design_files) == 0) stop("No design files found under ", data_root)

# Append-only log of every discrepancy detected while loading designs / files.
# process_plate writes into this; the script prints it before exiting.
.discrepancies <- list()

# Parse "T<N>h" out of an xlsx filename (returns NA if not found)
parse_timepoint <- function(fname) {
  m <- regmatches(fname, regexpr("T\\d+h", fname, ignore.case = TRUE))
  if (length(m) == 0) return(NA_integer_)
  as.integer(gsub("[A-Za-z]", "", m))
}

# Pull a likely antibiotic token from filename: Mero | MRP | Tob | TOB.
# Returned uppercased + canonicalised (Mero -> MRP, since both = meropenem).
parse_antibiotic <- function(fname) {
  m <- regmatches(fname, regexpr("Mero|MRP|Tob|TOB", fname, ignore.case = TRUE))
  if (length(m) == 0) return(NA_character_)
  ab <- toupper(m)
  if (ab == "MERO") ab <- "MRP"
  if (ab == "TOB")  ab <- "TOB"
  ab
}

plates <- lapply(design_files, function(design_path) {
  plate_dir  <- dirname(design_path)
  plate_name <- basename(plate_dir)
  experiment <- basename(dirname(plate_dir))

  xlsx <- list.files(plate_dir, pattern = "\\.xlsx$", full.names = TRUE)
  xlsx <- xlsx[!grepl("^~\\$", basename(xlsx))]  # drop Excel lock files

  files <- tibble(
    path       = xlsx,
    timepoint  = vapply(xlsx, parse_timepoint,  integer(1)),
    antibiotic = vapply(xlsx, parse_antibiotic, character(1))
  )

  list(
    experiment = experiment,
    plate      = plate_name,
    design     = design_path,
    files      = files
  )
})

cat("Found", length(plates), "plate folders:\n")
for (p in plates) {
  cat("  ", p$experiment, "/", p$plate,
      "  (", nrow(p$files), " xlsx)\n", sep = "")
}

# =============================================================================
# HELPERS
# =============================================================================

# Read an 8x12 plate-format Excel into a long Well/OD tibble.
plate_to_long <- function(file_path) {
  df <- suppressMessages(read_excel(file_path, col_names = FALSE))
  df <- df[, 1:min(13, ncol(df))]

  plate_data <- df[2:9, 2:13]
  colnames(plate_data) <- as.character(1:12)
  plate_data$row <- LETTERS[1:8]

  plate_data %>%
    pivot_longer(cols = -row, names_to = "col", values_to = "OD") %>%
    mutate(
      col  = as.integer(col),
      OD   = suppressWarnings(as.numeric(OD)),
      Well = paste0(row, col)
    ) %>%
    dplyr::select(Well, OD)
}

# Read + harmonise a design CSV. Handles, in order:
#   - BOM in headers (read_csv handles this)
#   - trailing/leading whitespace on header names (e.g. "Lineage ")
#   - alternate column names: AB -> Antibiotic, Linneage -> Lineage
#   - missing Lineage column (default "evolved")
#   - whitespace + typos in Type values ("emtpy ", "empty ")
#   - whitespace + spelling variation in Lineage ("ancestral" -> "ancestor")
#   - meropenem antibiotic variants (Mero -> MRP)
# Returns the cleaned design + a tibble of discrepancies found, attached as an
# attribute so the caller can roll them up into the startup audit report.
read_design <- function(file_path) {
  d <- read_csv(file_path, show_col_types = FALSE)
  names(d) <- trimws(names(d))

  notes <- character()

  if ("AB" %in% names(d) && !"Antibiotic" %in% names(d)) {
    d <- rename(d, Antibiotic = AB)
    notes <- c(notes, "header 'AB' renamed to 'Antibiotic'")
  }
  if ("Linneage" %in% names(d) && !"Lineage" %in% names(d)) {
    d <- rename(d, Lineage = Linneage)
    notes <- c(notes, "header 'Linneage' renamed to 'Lineage'")
  }
  if (!"Lineage" %in% names(d)) {
    d$Lineage <- "evolved"
    notes <- c(notes, "no Lineage column: defaulted every row to 'evolved'")
  }

  type_raw <- as.character(d$Type)
  if (any(grepl("emtpy", type_raw, fixed = TRUE), na.rm = TRUE))
    notes <- c(notes, "Type value 'emtpy' normalised to 'empty'")
  if (any(type_raw != trimws(type_raw), na.rm = TRUE))
    notes <- c(notes, "trailing whitespace in Type values trimmed")

  lin_raw <- as.character(d$Lineage)
  if (any(grepl("^ancestral\\s*$", lin_raw), na.rm = TRUE))
    notes <- c(notes, "Lineage value 'ancestral' normalised to 'ancestor'")
  if (any(lin_raw != trimws(lin_raw), na.rm = TRUE))
    notes <- c(notes, "trailing whitespace in Lineage values trimmed")

  d <- d %>%
    mutate(
      Type       = trimws(as.character(Type)),
      Type       = if_else(Type == "emtpy", "empty", Type),
      Lineage    = trimws(as.character(Lineage)),
      Lineage    = if_else(Lineage == "ancestral", "ancestor", Lineage),
      Lineage    = if_else(is.na(Lineage) | Lineage == "NA" | Lineage == "",
                           "evolved", Lineage),
      Antibiotic = toupper(trimws(as.character(Antibiotic))),
      Antibiotic = if_else(Antibiotic == "MERO", "MRP", Antibiotic),
      Sample     = trimws(as.character(Sample))
    )

  attr(d, "notes") <- notes
  d
}

# =============================================================================
# LOAD: stack every (experiment, plate, antibiotic, timepoint) into long form
# =============================================================================

process_plate <- function(plate_info) {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("Processing:", plate_info$experiment, "/", plate_info$plate, "\n", sep = "")

  design <- read_design(plate_info$design)
  for (n in attr(design, "notes")) {
    .discrepancies[[length(.discrepancies) + 1L]] <<- list(
      experiment = plate_info$experiment,
      plate      = plate_info$plate,
      issue      = n
    )
  }

  design_abs <- unique(na.omit(design$Antibiotic))

  files <- plate_info$files %>%
    filter(!is.na(timepoint), timepoint %in% mic_timepoints_h)

  # Match xlsx -> design antibiotic. Three cases:
  #   1. filename token matches a design antibiotic -> use it.
  #   2. filename has no token AND design has exactly one antibiotic -> use design's.
  #   3. filename token disagrees with design (e.g. file says MRP, design says
  #      TOB) -> antibiotic_conflict_policy decides.
  resolve_ab <- function(file_ab) {
    if (!is.na(file_ab) && file_ab %in% design_abs) return(file_ab)
    if (is.na(file_ab) && length(design_abs) == 1)  return(design_abs)
    if (!is.na(file_ab) && length(design_abs) >= 1) {
      # Real conflict
      msg <- sprintf("filename antibiotic '%s' not in design (%s)",
                     file_ab, paste(design_abs, collapse = ","))
      .discrepancies[[length(.discrepancies) + 1L]] <<- list(
        experiment = plate_info$experiment,
        plate      = plate_info$plate,
        issue      = msg
      )
      switch(antibiotic_conflict_policy,
             filename = file_ab,
             design   = NA_character_,
             stop     = stop(plate_info$experiment, "/", plate_info$plate,
                             ": ", msg))
    } else {
      NA_character_
    }
  }
  files <- files %>%
    mutate(antibiotic_match = vapply(antibiotic, resolve_ab, character(1)))

  skipped <- files %>% filter(is.na(antibiotic_match))
  if (nrow(skipped) > 0) {
    cat("  Skipping (no design antibiotic match):\n")
    for (p in skipped$path) cat("    -", basename(p), "\n")
  }

  files <- files %>% filter(!is.na(antibiotic_match))
  if (nrow(files) == 0) {
    cat("  No usable xlsx files; skipping plate.\n")
    return(NULL)
  }

  # Build the joined frame per file. For each file we pick the design rows
  # whose Antibiotic matches; if none match (because antibiotic_conflict_policy
  # = "filename" overrode), we use the full design and stamp the file's
  # antibiotic onto every row. The latter only fires on a real conflict that
  # was already logged above.
  # For each file: keep design rows whose Antibiotic matches the file's
  # antibiotic, PLUS blank/empty wells (which carry NA Antibiotic and apply
  # to every read). Then stamp the file's antibiotic onto blank/empty rows
  # so blank correction can group by (Plate, Antibiotic, Time_h).
  joined <- pmap_dfr(
    list(files$path, files$timepoint, files$antibiotic_match),
    function(path, t_h, ab) {
      reads <- plate_to_long(path) %>% mutate(Time_h = t_h)
      design_sub <- design %>%
        filter((!is.na(Antibiotic) & Antibiotic == ab) |
               Type %in% c("blank", "empty"))
      if (!any(design_sub$Type == "sample" & design_sub$Antibiotic == ab,
               na.rm = TRUE)) {
        # Conflict-policy = "filename" path: design has no rows for this AB.
        # Fall back to the entire design with Antibiotic overridden.
        design_sub <- design %>% mutate(Antibiotic = ab)
      } else {
        design_sub <- design_sub %>% mutate(Antibiotic = ab)
      }
      reads %>% inner_join(design_sub, by = "Well")
    }
  ) %>%
    mutate(
      Experiment = plate_info$experiment,
      Plate      = plate_info$plate
    )

  # Per-timepoint blank correction (per Plate x Antibiotic x Timepoint)
  blanks <- joined %>%
    filter(Type == "blank") %>%
    group_by(Plate, Antibiotic, Time_h) %>%
    summarize(blank_mean = mean(OD, na.rm = TRUE), .groups = "drop")

  joined <- joined %>%
    left_join(blanks, by = c("Plate", "Antibiotic", "Time_h")) %>%
    mutate(OD_blanked = OD - blank_mean)

  cat("  Wells:", n_distinct(joined$Well),
      " | Timepoints:", paste(sort(unique(joined$Time_h)), collapse = ","),
      " | Antibiotics:", paste(sort(unique(joined$Antibiotic)), collapse = ","), "\n")

  joined
}

all_data <- map_dfr(plates, process_plate)
if (nrow(all_data) == 0) stop("No data loaded.")

# Classify each sample into a cohort for grouped plots:
#   M-lines  = M1, M2, M3        (evolved under MRP)
#   T-lines  = T1, T2, T3        (evolved under TOB)
#   PA       = PA01_1, PA01_2    (PA01 reference / ancestor strains)
#   N-lines  = NA1, NA2, NA3     (evolved without antibiotic)
classify_sample <- function(s) {
  s <- as.character(s)
  case_when(
    grepl("^PA01",   s) ~ "PA strains",
    grepl("^NA\\d",  s) ~ "N-lines",
    grepl("^M\\d",   s) ~ "M-lines",
    grepl("^T\\d",   s) ~ "T-lines",
    TRUE                ~ NA_character_
  )
}

sample_group_levels <- c("M-lines", "T-lines", "N-lines", "PA strains")

all_data <- all_data %>%
  mutate(sample_group = factor(classify_sample(Sample),
                               levels = sample_group_levels))

# =============================================================================
# REPLICATE NUMBERING (within experiment x plate x sample x lineage x conc)
# =============================================================================

well_replicates <- all_data %>%
  filter(Type == "sample") %>%
  distinct(Experiment, Plate, Antibiotic, Sample, Lineage, Concentration, Well) %>%
  group_by(Experiment, Plate, Antibiotic, Sample, Lineage, Concentration) %>%
  arrange(Well, .by_group = TRUE) %>%
  mutate(replicate = row_number()) %>%
  ungroup()

all_data <- all_data %>%
  left_join(
    well_replicates,
    by = c("Experiment", "Plate", "Antibiotic", "Sample", "Lineage",
           "Concentration", "Well")
  )

# =============================================================================
# MIC PER REPLICATE (handles censoring: every conc still > cutoff -> NA + flag)
# =============================================================================

treatment <- all_data %>%
  filter(Type == "sample", Concentration > 0)

mic_per_rep <- treatment %>%
  mutate(timepoint = paste0(Time_h, "h")) %>%
  group_by(Experiment, Plate, Antibiotic, Sample, sample_group, Lineage,
           replicate, timepoint) %>%
  arrange(Concentration, .by_group = TRUE) %>%
  summarize(
    max_conc_tested = max(Concentration, na.rm = TRUE),
    mic = {
      below <- which(OD_blanked < od_cutoff)
      if (length(below) > 0) Concentration[min(below)] else NA_real_
    },
    censored = is.na(mic),
    .groups = "drop"
  )

# Censoring diagnostic
censoring <- mic_per_rep %>%
  group_by(Experiment, Sample, Antibiotic, Lineage, timepoint) %>%
  summarize(
    n_total    = n(),
    n_censored = sum(censored),
    .groups = "drop"
  )

cat("\n", strrep("=", 60), "\n", sep = "")
cat("MIC censoring (replicates where no conc reached OD < ", od_cutoff, ")\n",
    strrep("=", 60), "\n\n", sep = "")
print(censoring, n = Inf)

mic_summary <- mic_per_rep %>%
  group_by(Experiment, Antibiotic, Sample, sample_group, Lineage, timepoint) %>%
  summarize(
    n_replicates       = n(),
    n_censored         = sum(censored),
    max_conc_tested    = max(max_conc_tested, na.rm = TRUE),
    mic_min            = suppressWarnings(min(mic, na.rm = TRUE)),
    mic_max            = suppressWarnings(max(mic, na.rm = TRUE)),
    mic_median         = median(mic, na.rm = TRUE),
    mic_mean           = mean(mic,   na.rm = TRUE),
    mic_geometric_mean = exp(mean(log(mic), na.rm = TRUE)),
    mic_geo_se_lower   = exp(mean(log(mic), na.rm = TRUE) -
                             sd(log(mic), na.rm = TRUE) / sqrt(sum(!is.na(mic)))),
    mic_geo_se_upper   = exp(mean(log(mic), na.rm = TRUE) +
                             sd(log(mic), na.rm = TRUE) / sqrt(sum(!is.na(mic)))),
    .groups = "drop"
  ) %>%
  mutate(across(c(mic_geometric_mean, mic_geo_se_lower, mic_geo_se_upper,
                  mic_mean, mic_min, mic_max, mic_median),
                ~ ifelse(is.finite(.x), .x, NA_real_))) %>%
  arrange(Experiment, Antibiotic, Sample, Lineage, timepoint)

# Drop the T0 row from MIC tables (T0 is a baseline read, not a kill timepoint)
mic_per_rep <- mic_per_rep %>% filter(timepoint != "0h")
mic_summary <- mic_summary %>% filter(timepoint != "0h")

# =============================================================================
# RELATIVE MIC (evolved / ancestor), matched within antibiotic + timepoint
# =============================================================================

relative_mic <- NULL
lineage_levels <- unique(mic_summary$Lineage)
if (any(grepl("ancestor|ancestral", lineage_levels, ignore.case = TRUE)) &&
    any(grepl("evolved",            lineage_levels, ignore.case = TRUE))) {

  ancestor_mic <- mic_summary %>%
    filter(grepl("ancestor|ancestral", Lineage, ignore.case = TRUE)) %>%
    group_by(Experiment, Antibiotic, Sample, timepoint) %>%
    summarize(ancestor_mic = exp(mean(log(mic_geometric_mean), na.rm = TRUE)),
              .groups = "drop")

  relative_mic <- mic_summary %>%
    filter(grepl("evolved", Lineage, ignore.case = TRUE)) %>%
    left_join(ancestor_mic,
              by = c("Experiment", "Antibiotic", "Sample", "timepoint")) %>%
    mutate(relative_mic = mic_geometric_mean / ancestor_mic)
}

# =============================================================================
# PLOTS — same vocabulary as the kinetic script
# =============================================================================

sample_data <- all_data %>% filter(Type == "sample")

# 1. Dose-response at each MIC timepoint, ALL experiments overlaid.
#    Antibiotic + Sample down rows, timepoint across columns.
dose_resp <- sample_data %>%
  filter(Time_h != 0) %>%
  group_by(Experiment, Antibiotic, Sample, Lineage, Time_h, Concentration) %>%
  summarize(
    mean_OD = mean(OD_blanked, na.rm = TRUE),
    se_OD   = sd(OD_blanked,   na.rm = TRUE) / sqrt(sum(!is.na(OD_blanked))),
    .groups = "drop"
  ) %>%
  mutate(timepoint = factor(paste0(Time_h, "h"),
                            levels = paste0(sort(setdiff(mic_timepoints_h, 0)), "h")))

min_pos <- dose_resp %>%
  filter(Concentration > 0) %>%
  pull(Concentration) %>% { suppressWarnings(min(., na.rm = TRUE)) }
if (!is.finite(min_pos)) min_pos <- 0.01

p_dr <- ggplot(dose_resp,
               aes(x = pmax(Concentration, min_pos / 3),
                   y = mean_OD,
                   color = Lineage,
                   shape = Experiment,
                   group = interaction(Experiment, Lineage))) +
  geom_line(alpha = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_OD - se_OD, ymax = mean_OD + se_OD),
                width = 0.08) +
  geom_hline(yintercept = od_cutoff, linetype = "dashed", color = "gray40") +
  scale_x_log10() +
  facet_grid(Antibiotic + Sample ~ timepoint, scales = "free_y") +
  labs(title = "Dose-response at MIC timepoints (mean +/- SE)",
       subtitle = paste("Dashed line = od_cutoff =", od_cutoff,
                        "  |  zero-conc plotted at 1/3 of lowest non-zero"),
       x = "Concentration (ug/mL, log)", y = "OD600 (blanked)") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "dose_response_at_timepoints.pdf"),
       p_dr, width = 12, height = 10)

# 2. MIC summary grouped by sample cohort:
#    rows = sample_group (M-lines, T-lines, N-lines, PA strains),
#    cols = timepoint, x = Sample, color = Lineage, shape = Antibiotic.
#    Lets you read each cohort's behaviour on its own row.
mic_plot_data <- mic_summary %>%
  mutate(
    fully_censored = (n_censored == n_replicates),
    plot_y    = ifelse(fully_censored, max_conc_tested, mic_geometric_mean),
    plot_low  = ifelse(fully_censored, NA_real_,         mic_geo_se_lower),
    plot_high = ifelse(fully_censored, NA_real_,         mic_geo_se_upper),
    timepoint = factor(timepoint,
                       levels = paste0(sort(setdiff(mic_timepoints_h, 0)), "h"))
  )

lineage_colors <- c(ancestor = "darkblue", evolved = "firebrick")

# Drop combinations not present in the data so facet_wrap doesn't reserve
# blank panels (e.g. M-lines were never tested under TOB).
mic_plot_data <- mic_plot_data %>%
  mutate(panel = paste0(sample_group, " | ", Antibiotic)) %>%
  mutate(panel = factor(panel, levels = unique(panel[order(sample_group, Antibiotic)])))

p_mic <- ggplot(mic_plot_data,
                aes(x = Sample, y = plot_y,
                    color = Lineage, shape = timepoint,
                    group = interaction(Lineage, timepoint))) +
  geom_point(position = position_dodge(width = 0.6), size = 2.8) +
  geom_errorbar(aes(ymin = plot_low, ymax = plot_high),
                position = position_dodge(width = 0.6),
                width = 0.25, na.rm = TRUE) +
  scale_y_log10() +
  scale_color_manual(values = lineage_colors) +
  scale_shape_manual(values = c(`16h` = 16, `24h` = 17)) +
  facet_wrap(~ panel, scales = "free_x", drop = TRUE) +
  labs(title = "MIC by sample cohort (geometric mean +/- SE)",
       subtitle = paste("Cohorts: M-lines = MRP-evolved, T-lines = TOB-evolved,",
                        "N-lines = no-antibiotic, PA strains = PA01 reference.",
                        "\nPoints at max concentration = censored (no replicate reached cutoff)."),
       x = "Sample", y = "MIC (ug/mL)",
       shape = "Timepoint") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x     = element_text(angle = 30, hjust = 1))

ggsave(file.path(out_dir, "MIC_across_experiments.pdf"),
       p_mic, width = 11, height = 7)

# 3. Per-plate "as-laid-out" diagnostic: each well in its row x col position,
#    showing OD across the three timepoints with the cutoff line.
plot_plate_layout <- function(plate_data, plate_label) {
  pd <- plate_data %>%
    mutate(
      row_letter = factor(gsub("[0-9]+", "", Well), levels = LETTERS[1:8]),
      col_num    = factor(as.integer(gsub("[A-H]", "", Well)), levels = 1:12)
    )

  annot <- pd %>%
    distinct(Well, row_letter, col_num,
             Sample, Concentration, Lineage, Type, Antibiotic) %>%
    mutate(
      label = case_when(
        Type == "blank"  ~ "blank",
        Type == "empty"  ~ "empty",
        Type == "sample" & Concentration == 0 ~ paste0(Lineage, "\nctrl"),
        Type == "sample" ~ paste0(Sample, "\n", Concentration),
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

  pd <- pd %>% left_join(annot %>% dplyr::select(Well, well_color), by = "Well")

  ggplot(pd, aes(x = Time_h, y = OD_blanked, group = Well)) +
    geom_hline(yintercept = od_cutoff, linetype = "dashed",
               color = "firebrick", linewidth = 0.3) +
    geom_line(aes(color = well_color), linewidth = 0.5) +
    geom_point(aes(color = well_color), size = 0.8) +
    geom_text(data = annot, aes(label = label),
              x = Inf, y = Inf, hjust = 1.02, vjust = 1.15,
              size = 1.7, color = "gray20", inherit.aes = FALSE) +
    scale_color_manual(values = c(
      ancestor   = "#1f77b4",
      ancestral  = "#1f77b4",
      evolved    = "#d62728",
      control    = "darkgreen",
      blank      = "gray70",
      empty      = "gray85",
      other      = "black"
    ), name = NULL) +
    facet_grid(row_letter ~ col_num, switch = "y") +
    labs(title  = paste("Plate layout (endpoint reads):", plate_label),
         subtitle = paste("Dashed red = od_cutoff =", od_cutoff,
                          "  |  x-axis = Time_h (",
                          paste(sort(mic_timepoints_h), collapse = ","), ")"),
         x = "Time (h)", y = "OD600 (blanked)") +
    theme_bw() +
    theme(strip.text   = element_text(size = 7),
          axis.text    = element_text(size = 5),
          panel.spacing = unit(0.15, "lines"),
          legend.position = "bottom")
}

plate_groups <- all_data %>%
  distinct(Experiment, Plate, Antibiotic)

for (i in seq_len(nrow(plate_groups))) {
  pg <- plate_groups[i, ]
  pd <- all_data %>%
    filter(Experiment == pg$Experiment,
           Plate      == pg$Plate,
           Antibiotic == pg$Antibiotic)
  plate_label <- paste(pg$Experiment, pg$Plate, pg$Antibiotic, sep = " | ")
  fname <- paste0("plate_layout_",
                  pg$Experiment, "_", pg$Plate, "_", pg$Antibiotic, ".pdf")
  ggsave(file.path(out_dir, fname),
         plot_plate_layout(pd, plate_label),
         width = 16, height = 10)
}

# 4. Heatmap of OD at 24h per well (one panel per plate x antibiotic)
heat_data <- all_data %>%
  filter(Time_h == 24, Type == "sample") %>%
  mutate(
    row_letter = factor(gsub("[0-9]+", "", Well), levels = rev(LETTERS[1:8])),
    col_num    = factor(as.integer(gsub("[A-H]", "", Well)), levels = 1:12)
  )

p_heat <- ggplot(heat_data,
                 aes(x = col_num, y = row_letter, fill = OD_blanked)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(OD_blanked, 2)), size = 1.8) +
  scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue",
                       midpoint = 0.2, name = "OD600") +
  facet_wrap(~ Experiment + Plate + Antibiotic) +
  labs(title = "Plate heatmaps at 24h (blanked OD)",
       x = "Column", y = "Row") +
  theme_minimal() +
  coord_equal()

ggsave(file.path(out_dir, "MIC_plate_heatmaps_24h.pdf"),
       p_heat, width = 14, height = 10)

# =============================================================================
# PRINT + EXPORT
# =============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("MIC SUMMARY (across all experiments)\n", strrep("=", 60), "\n\n", sep = "")
print(mic_summary, n = Inf)

if (!is.null(relative_mic)) {
  cat("\nRELATIVE MIC (evolved / ancestor)\n", strrep("-", 40), "\n", sep = "")
  print(relative_mic %>%
          dplyr::select(Experiment, Antibiotic, Sample, timepoint,
                        mic_geometric_mean, ancestor_mic, relative_mic),
        n = Inf)
}

write_csv(all_data,    file.path(out_dir, "endpoint_all_wells_long.csv"))
write_csv(mic_per_rep, file.path(out_dir, "MIC_per_replicate.csv"))
write_csv(mic_summary, file.path(out_dir, "MIC_summary.csv"))
if (!is.null(relative_mic)) {
  write_csv(relative_mic, file.path(out_dir, "MIC_relative.csv"))
}

# Discrepancy audit: every cleanup applied during load lives in .discrepancies.
# Print at the end so it's the last thing on screen and easy to act on.
if (length(.discrepancies) > 0) {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("DESIGN / FILE DISCREPANCIES (",
      length(.discrepancies), " total)\n", sep = "")
  cat(strrep("=", 60), "\n", sep = "")
  audit <- bind_rows(.discrepancies) %>% distinct()
  print(audit, n = Inf)
  write_csv(audit, file.path(out_dir, "discrepancies.csv"))
}

cat("\nDone. Outputs in:", out_dir, "\n")
