# =========================================================
# 1. SETUP
# =========================================================
library(tidyverse)
library(ggpubr)
library(emmeans)
library(ggrepel)
library(ggplot2)

# Load Data
df <- read.csv("FebInsectNutrition.csv", check.names = FALSE)

# =========================================================
# 2. DATA CALCULATION -> Accounting for multiple individuals per sample
# =========================================================
df_clean <- df %>%
  select(Order, Family, Month, Quantity, Size, Sample_ID, Insect_Weight, Dry_Weight, Lipid_Free_Weight) %>%
  filter(!is.na(Order), !is.na(Family), !is.na(Month), 
         !is.na(Insect_Weight), !is.na(Dry_Weight), !is.na(Lipid_Free_Weight)) %>%
  mutate(
    # Format Month order
    Month = factor(Month, levels = c("May", "June", "July", "August", "September")),
    
    # Create Family Label (Preserving your Size Groups for plotting/analysis)
    Size_Label = ifelse(Size == "" | is.na(Size), Family, paste0(Family, " (", Size, ")")),
    Taxon_Group = Size_Label,
    
    # METRICS
    Mass_Ind = Insect_Weight / Quantity,
    Lipid_Total = (Dry_Weight - Lipid_Free_Weight) / Quantity,
    Lipid_Rel_Pct = ((Dry_Weight - Lipid_Free_Weight) / Insect_Weight)
  ) %>%
  # Filter valid data (>0 is required for Log transformation)
  filter(Mass_Ind > 0, Lipid_Total > 0, Lipid_Rel_Pct > 0) %>%
  mutate(
    # --- LOG TRANSFORMATIONS (FOR ALL 3) ---
    Log_Mass = log(Mass_Ind),
    Log_Lipid_Total = log(Lipid_Total),
    Logit_Lipid_Rel = log(Lipid_Rel_Pct)
  )

# =========================================================
# ADD LIFE HISTORY COLUMN
# =========================================================

df_clean <- df_clean %>%
  mutate(
    Life_History = case_when(
      # Slow Life History (Long-lived, high investment)
      Order %in% c("Lepidoptera", "Coleoptera") ~ "Slow",
      
      # Fast Life History (Short-lived, ephemeral)
      Order %in% c("Trichoptera", "Diptera", "Ephemeroptera") ~ "Fast",
      
      # Catch-all for any other orders not listed
      TRUE ~ "Other"
    )
  )

# Check to ensure it worked
print(table(df_clean$Order, df_clean$Life_History))

# =========================================================
# 3. NORMALITY TESTING (Raw vs. Log)
# =========================================================

print("=== SHAPIRO-WILK NORMALITY TESTS ===")

# A. Mass (Size)
print("--- 1. Individual Mass ---")
print(paste("Raw Data p:", format(shapiro.test(df_clean$Mass_Ind)$p.value, digits=3)))
print(paste("Log Data p:", format(shapiro.test(df_clean$Log_Mass)$p.value, digits=3)))

# B. Total Lipid Content
print("--- 2. Total Lipid Content ---")
print(paste("Raw Data p:", format(shapiro.test(df_clean$Lipid_Total)$p.value, digits=3)))
print(paste("Log Data p:", format(shapiro.test(df_clean$Log_Lipid_Total)$p.value, digits=3)))

# C. Relative Lipid Content
print("--- 3. Relative Lipid Content ---")
print(paste("Raw Data p:", format(shapiro.test(df_clean$Lipid_Rel_Pct)$p.value, digits=3)))
print(paste("Log Data p:", format(shapiro.test(df_clean$Logit_Lipid_Rel)$p.value, digits=3)))

##log didn't help, inspecting with qq
# Set up a 3x2 grid to view all plots at once
par(mfrow=c(3,2))

# --- 1. INDIVIDUAL MASS ---
# Raw
qqnorm(df_clean$Mass_Ind, main="Q-Q Plot: Raw Mass")
qqline(df_clean$Mass_Ind, col="red", lwd=2)

# Log
qqnorm(df_clean$Log_Mass, main="Q-Q Plot: Log Mass")
qqline(df_clean$Log_Mass, col="red", lwd=2)

# --- 2. TOTAL LIPID CONTENT ---
# Raw
qqnorm(df_clean$Lipid_Total, main="Q-Q Plot: Raw Total Lipid")
qqline(df_clean$Lipid_Total, col="red", lwd=2)

# Log
qqnorm(df_clean$Log_Lipid_Total, main="Q-Q Plot: Log Total Lipid")
qqline(df_clean$Log_Lipid_Total, col="red", lwd=2)

# --- 3. RELATIVE LIPID CONTENT ---
# Raw
qqnorm(df_clean$Lipid_Rel_Pct, main="Q-Q Plot: Raw Relative Lipid")
qqline(df_clean$Lipid_Rel_Pct, col="red", lwd=2)

# Log
qqnorm(df_clean$Logit_Lipid_Rel, main="Q-Q Plot: Log Relative Lipid")
qqline(df_clean$Logit_Lipid_Rel, col="red", lwd=2)

# Reset plot layout to default
par(mfrow=c(1,1))

###Relative lipid might beneift from logit, or could leave?

# Option B: Logit Transformation
# We add a tiny constant (0.01) to avoid infinite values if you have exact 0% or 100%
df_clean$Logit_Lipid_Rel <- log((df_clean$Lipid_Rel_Pct / 100 + 0.001) / (1 - (df_clean$Lipid_Rel_Pct / 100) + 0.001))

# 2. Run Shapiro-Wilk Tests
print("--- Shapiro-Wilk Results ---")
print(paste("Logit p:  ", format(shapiro.test(df_clean$Logit_Lipid_Rel)$p.value, digits=3)))

# Plot B: Logit
qqnorm(df_clean$Logit_Lipid_Rel, main="Q-Q Plot: Logit Transformed")
qqline(df_clean$Logit_Lipid_Rel, col="blue", lwd=2)


# =========================================================
# ADD Mouth column
# =========================================================
# 1. CREATE THE MOUTH LOOKUP TABLE
# Transcribed exactly from your updated image (image_2dfeff.png)
mouth_data <- tribble(
  ~Family,          ~Mouth_Type,
  "Geometridae",    "yes",
  "Crambidae",      "yes",
  "Noctuidae",      "yes",
  "Lasiocampidae",  "no",
  "Nematocera",     "yes",
  "Carabidae",      "yes",
  "Brachycera",     "yes",
  "Leptoceridae",   "yes",
  "Notodontidae",   "no",
  "Elateridae",     "yes",
  "Erebidae",       "mixed",
  "Pyralidae",      "yes",
  "Microlep",       "yes",
  "Gelechiidae",    "yes",
  "Caenidae",       "no",
  "Hydroptilidae",  "yes",
  "Limnephilidae",  "yes",
  "Heptageniidae",  "no",
  "Ephemeridae",    "no",
  "Cosmopterigidae","yes",
  "Heteroceridae",  "yes",
  "Staphylinidae",  "yes",
  "Tortricidae",    "yes",
  "Silphidae",      "yes",
  "Phryganeidae",   "no",
  "Helicopsychidae", "yes",
)

# =========================================================
# 4. TWO-WAY ANOVA (Order + Month)
# =========================================================
# --- Order Level Models ---
mod_ord_mass  <- lm(Log_Mass ~ Order + Mouth_Type + Month, data = df_clean)
mod_ord_total <- lm(Log_Lipid_Total ~ Order + Mouth_Type + Month, data = df_clean)
mod_ord_rel   <- lm(Logit_Lipid_Rel ~ Order + Mouth_Type + Month, data = df_clean)

# Print Results
print("=== ANOVA RESULTS (Two-Way: Taxonomy + Month) ===")
print("--- 1. Size (Log Mass) ---")
print(anova(mod_ord_mass))

print("--- 2. Total Lipid Content (Log) ---")
print(anova(mod_ord_total))

print("--- 3. Relative Lipid Content (Log) ---")
print(anova(mod_ord_rel))

# --- B. Create Filtered Data (Remove Singletons so we can run anova on families) ---
counts <- df_clean %>%
  group_by(Taxon_Group) %>%
  tally()

valid_groups <- counts %>%
  filter(n >= 2) %>%
  pull(Taxon_Group)

# Create the specific dataset for Family Analysis
df_fam_valid <- df_clean %>%
  filter(Taxon_Group %in% valid_groups)

#Check to see how many families were retained
print(paste("Original Families:", nrow(counts)))
print(paste("Valid Families (N>=2):", length(valid_groups)))

# --- C. Family Level Models (Run on FILTERED Data) ---
# Note: We name these '_valid' so they match your post-hoc code later
mod_fam_mass_valid  <- lm(Log_Mass ~  Mouth_Type + Taxon_Group + Month, data = df_fam_valid)
mod_fam_total_valid <- lm(Log_Lipid_Total ~ Mouth_Type + Taxon_Group + Month, data = df_fam_valid)
mod_fam_rel_valid   <- lm(Logit_Lipid_Rel ~ Mouth_Type + Taxon_Group + Month, data = df_fam_valid)

print("=== ANOVA RESULTS (Family Level: Valid N > 1) ===")
print("--- 2. Family Level ---")
print(anova(mod_fam_mass_valid))
print(anova(mod_fam_total_valid))
print(anova(mod_fam_rel_valid))
summary(mod_fam_mass_valid)

# =========================================================
# 5. SCIENTIFIC PLOTTING (Boxplots)
# =========================================================

# We use the RAW data for the axes, but apply 'scale_y_log10'

create_boxplot <- function(data, x_var, y_var, fill_var, title_text, y_label) {
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[fill_var]])) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    coord_flip() +
    scale_y_log10() + 
    theme_bw() +
    theme(legend.position = "none",
          axis.text.y = element_text(face = "bold")) +
    labs(title = title_text, y = y_label, x = "")
}

# 1. Mass (Size)
p_ord_mass <- create_boxplot_clean(
  data = df_clean, 
  x_var = "Order", 
  y_var = "Mass_Ind", 
  fill_var = "Order",
  title_text = "Order: Individual Size", 
  y_label = "Mass (g) - Log Scale",
  model = mod_ord_mass
)

# 2. Total Lipid
p_ord_total <- create_boxplot_clean(
  data = df_clean, 
  x_var = "Order", 
  y_var = "Lipid_Total", 
  fill_var = "Order",
  title_text = "Order: Total Lipid Content", 
  y_label = "Total Lipid (g) - Log Scale",
  model = mod_ord_total
)

# 3. Relative Lipid
p_ord_rel <- create_boxplot_clean(
  data = df_clean, 
  x_var = "Order", 
  y_var = "Lipid_Rel_Pct", 
  fill_var = "Order",
  title_text = "Order: Relative Lipid Content", 
  y_label = "Relative Lipid (Proportion) - Log Scale",
  model = mod_ord_rel
)

# Print
print(p_ord_mass)
print(p_ord_total)
print(p_ord_rel)

# --- Family Level Plots Unordered, better for comparisions between groups---
p_fam_mass <- create_boxplot(df_clean, "Sample_ID", "Mass_Ind", "Order", 
                             "Family: Individual Size", "Mass (g) - Log Scale") + 
  facet_grid(Order ~ ., scales = "free_y", space = "free_y")

p_fam_total <- create_boxplot(df_clean, "Sample_ID", "Lipid_Total", "Order", 
                              "Family: Total Lipid Content", "Total Lipid (g) - Log Scale") + 
  facet_grid(Order ~ ., scales = "free_y", space = "free_y")

p_fam_rel <- create_boxplot(df_clean, "Sample_ID", "Lipid_Rel_Pct", "Order", 
                            "Family: Relative Lipid Content", "Relative Lipid (Proportion) - Log Scale") + 
  facet_grid(Order ~ ., scales = "free_y", space = "free_y")

print(p_fam_mass)
print(p_fam_total)
print(p_fam_rel)

# =========================================================
# 6. POST-HOC COMPARISONS (Pairwise Tukey Tests)
# =========================================================

# --- Helper Function: Show ALL Pairs ---
print_all_contrasts <- function(model, factor_name, metric_label) {
  cat(paste0("\n>>> PAIRWISE COMPARISONS (ALL): ", metric_label, " <<<\n"))
  
  emm <- emmeans(model, specs = factor_name)
  
  # We use dplyr::select to prevent conflicts with other packages
  res <- as.data.frame(pairs(emm, adjust = "tukey")) %>%
    mutate(
      estimate = round(estimate, 3),
      p.value  = format.pval(p.value, digits = 3)
    ) %>%
    dplyr::select(contrast, estimate, p.value)
  
  print(res)
  
  cat("------------------------------------------------------\n")
}

# --- Helper Function: Show SIGNIFICANT Pairs Only ---
print_significant_contrasts <- function(model, factor_name, metric_label) {
  cat(paste0("\n>>> PAIRWISE COMPARISONS (SIGNIFICANT): ", metric_label, " <<<\n"))
  
  emm <- emmeans(model, specs = factor_name)
  
  res <- as.data.frame(pairs(emm, adjust = "tukey")) %>%
    filter(p.value < 0.05) %>%
    mutate(
      estimate = round(estimate, 3),
      p.value  = format.pval(p.value, digits = 3)
    ) %>%
    dplyr::select(contrast, estimate, p.value)
  
  if(nrow(res) > 0) {
    print(res)
  } else {
    cat("No significant pairwise differences found.\n")
  }
  cat("------------------------------------------------------\n")
}

# --- A. ORDER LEVEL (Show All Pairs) ---
print("=== 1. ORDER COMPARISONS ===")
print_all_contrasts(mod_ord_mass,  "Order", "Size")
print_all_contrasts(mod_ord_total, "Order", "Total Lipid")
print_all_contrasts(mod_ord_rel,   "Order", "Relative Lipid")

# --- B. FAMILY LEVEL (Show Significant Only) ---
print("=== 2. FAMILY COMPARISONS (Valid N>1) ===")
# Ensure you run this on the VALID family models from Step 5
print_significant_contrasts(mod_fam_mass_valid,  "Taxon_Group", "Size")
print_significant_contrasts(mod_fam_total_valid, "Taxon_Group", "Total Lipid")
print_significant_contrasts(mod_fam_rel_valid,   "Taxon_Group", "Relative Lipid")

# =========================================================
# 1. HELPER: EXTRACT SIGNIFICANT PAIRS
# =========================================================
get_tukey_data <- function(model, x_col) {
  # Get Emmeans
  emm <- emmeans(model, specs = x_col)
  
  # Get significant pairs only
  stats <- as.data.frame(pairs(emm, adjust = "tukey")) %>%
    filter(p.value < 0.05) %>%
    mutate(
      # Simplified Label: Just a "*" for any significant difference
      label = "*",
      # Ensure contrast is character to split it
      contrast = as.character(contrast)
    ) %>%
    separate(contrast, into = c("group1", "group2"), sep = " - ")
  
  return(stats)
}

# =========================================================
# 2. UPDATED PLOTTING FUNCTION (Auto-Stacking)
# =========================================================
create_boxplot_clean <- function(data, x_var, y_var, fill_var, title_text, y_label, model = NULL) {
  
  # 1. Calculate the starting height for brackets
  # We start drawing 10% above the highest data point
  max_val <- max(data[[y_var]], na.rm = TRUE)
  start_y <- max_val * 1.1 
  
  # 2. Base Plot
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[fill_var]])) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    
    # Log Scale
    scale_y_log10() + 
    
    # Theme
    theme_bw() +
    theme(legend.position = "none",
          axis.text.y = element_text(face = "bold", size = 10)) +
    labs(title = title_text, y = y_label, x = "") +
    coord_flip() 
  
  # 3. Add Brackets (If model provided)
  if (!is.null(model)) {
    stat.test <- get_tukey_data(model, x_var)
    
    if (nrow(stat.test) > 0) {
      # We assign the SAME y.position to all brackets initially
      stat.test$y.position <- start_y
      
      # We use 'step.increase' to automatically stack them 
      # 0.1 adds 10% spacing between each bracket
      p <- p + stat_pvalue_manual(
        stat.test, 
        label = "label", 
        y.position = "y.position",
        step.increase = 0.12,    # <--- This fixes the overlapping!
        coord.flip = TRUE,       # Essential for flipped plots
        size = 6,                # Size of the star
        bracket.size = 0.8,
        tip.length = 0.01
      )
      
      # Expand the plot limits so the highest bracket doesn't get cut off
      # We estimate how much space we need based on number of comparisons
      expansion_factor <- 1.1 + (nrow(stat.test) * 0.15)
      p <- p + expand_limits(y = max_val * expansion_factor)
    }
  }
  
  return(p)
}


##Paired family comparisons 
# --- C. Family Level Models (Run on FILTERED Data) ---
# Removing the life_history so we can visual straight forward comparisons 
mod_fam_mass_valid1  <- lm(Log_Mass ~ Taxon_Group + Month, data = df_fam_valid)
mod_fam_total_valid1 <- lm(Log_Lipid_Total ~  Taxon_Group + Month, data = df_fam_valid)
mod_fam_rel_valid1   <- lm(Logit_Lipid_Rel ~ Taxon_Group + Month, data = df_fam_valid)

# Function to create a clean Significance Heatmap
create_significance_heatmap <- function(model, factor_name, title_text) {
  
  # 1. Get Pairwise Comparisons
  emm <- emmeans(model, specs = factor_name)
  pairs_data <- as.data.frame(pairs(emm, adjust = "tukey")) %>%
    separate(contrast, into = c("Group1", "Group2"), sep = " - ") %>%
    mutate(
      p_bin = ifelse(p.value < 0.05, "Significant", "NS"),
      # Create a clean label for the tooltip/plot
      star = ifelse(p.value < 0.001, "***", 
                    ifelse(p.value < 0.01, "**", 
                           ifelse(p.value < 0.05, "*", "")))
    )
  
  # 2. Plot Matrix
  ggplot(pairs_data, aes(x = Group1, y = Group2, fill = p_bin)) +
    geom_tile(color = "white") +
    geom_text(aes(label = star), vjust = 0.7) + # Adds stars to significant squares
    scale_fill_manual(values = c("Significant" = "#FF6B6B", "NS" = "grey90"), 
                      name = "Tukey P-Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank()
    ) +
    labs(title = title_text, x = "", y = "") +
    coord_fixed()
}

# Generate Heatmaps for your Valid Family Models
p_matrix_mass <- create_significance_heatmap(mod_fam_mass_valid1, "Taxon_Group", "Mass Significance Matrix")
p_matrix_total <- create_significance_heatmap(mod_fam_total_valid1, "Taxon_Group", "Total Lipid Significance Matrix")
p_matrix_rel <- create_significance_heatmap(mod_fam_rel_valid1, "Taxon_Group", "Relative Lipid Content: Significance Matrix")

print(p_matrix_mass)
print(p_matrix_total)
print(p_matrix_rel)


# =========================================================
# Biplot
# =========================================================

# 1. PREPARE THE DATA 
df_tradeoff <- df_clean %>%
  # Add Life_History to the grouping so it stays in the summary
  group_by(Order, Taxon_Group, Life_History) %>%
  summarise(
    Mean_Mass = mean(Mass_Ind, na.rm = TRUE),
    Mean_Quality = mean(Lipid_Rel_Pct, na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  ) %>%
  filter(Count >= 2)

# 2. GENERATE THE PLOT
plot_tradeoff_circles <- ggplot(df_tradeoff, aes(x = Mean_Mass, y = Mean_Quality, color = Order)) +
  
  # A. The Grouping Circles (Ellipses)
  # We use 'stat_ellipse' to draw a circle around the Life History groups
  # alpha = 0.1 makes the fill very faint so it doesn't hide the points

  
  # B. The Points
  geom_point(size = 4, alpha = 0.9) +
  
  # C. The Labels
  geom_text_repel(aes(label = Taxon_Group), 
                  size = 3, 
                  box.padding = 0.5, 
                  max.overlaps = 30,
                  show.legend = FALSE) +
  
  # D. Scales & Colors
  scale_x_log10() + 
  
  # E. Reference Line
  geom_hline(yintercept = mean(df_clean$Lipid_Rel_Pct, na.rm=TRUE), 
             linetype = "dashed", color = "grey50") +
  
  # F. Formatting
  theme_bw() +
  labs(
    title = "The Trade-off: Quantity vs. Quality (Grouped by Life History)",
    x = "Quantity: Individual Mass (g) [Log Scale]",
    y = "Quality: Relative Lipid Content (Proportion)",
    fill = "Life History",
    linetype = "Life History"
  ) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# Print
print(plot_tradeoff_circles)

#####Size driver

# 1. RUN THE LINEAR MODEL
# We use the Log-Log version because biological scaling is almost always power-law.
mod_allometry <- lm(Log_Lipid_Total ~ Log_Mass, data = df_clean)

print("=== MODEL SUMMARY: DOES SIZE PREDICT FAT? ===")
summary(mod_allometry)

# 2. VISUALIZE THE RELATIONSHIP
ggplot(df_clean, aes(x = Log_Mass, y = Log_Lipid_Total)) +
  
  # A. The Points (Colored by Order to see if groups differ)
  geom_point(aes(color = Order), alpha = 0.6, size = 2) +
  
  # B. The Regression Line (The "Slope")
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  
  # C. Identity Line (Slope = 1) for comparison
  # If the black dashed line is steeper than this gray line, big bugs are "fatter."
  geom_abline(slope = 1, intercept = coef(mod_allometry)[1], color = "grey", linetype = "dotted") +
  
  theme_bw() +
  labs(
    title = "Allometry of Fat Storage",
    subtitle = "Black Line = Actual Trend. Dotted Grey = 1:1 Scaling (Isometric).",
    x = "Insect Size (Log Mass)",
    y = "Total Lipid Content (Log Mass)",
    color = "Order"
  )



##Mass predicting relative lipid


library(tidyverse)
library(scales) # For percentage labels

# 1. RUN THE LINEAR MODEL
# We use Log_Weight (Size) to predict Lipid_Rel_Pct (Quality).
# If the slope is positive, big bugs are better. If negative, big bugs are worse.
mod_tradeoff <- lm(Logit_Lipid_Rel ~ Log_Mass, data = df_clean)

print("=== MODEL SUMMARY: DOES MASS PREDICT QUALITY? ===")
summary(mod_tradeoff)

# 2. VISUALIZE THE RELATIONSHIP
ggplot(df_clean, aes(x = Log_Mass, y = Logit_Lipid_Rel)) +
  
  # A. The Points (Colored by Order to spot group differences)
  geom_point(aes(color = Order), alpha = 0.6, size = 2) +
  
  # B. The Regression Line (The "Slope")
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  
  # C. Scales & Formatting
  scale_y_continuous(labels = scales::percent) + # Show Y-axis as 5%, 10%
  
  theme_bw() +
  labs(
    title = "Quality-Quantity Trade-off",
    subtitle = "Does being larger (Mass) mean being fattier (Relative Lipid)?",
    x = "Insect Size (Log Mass)",
    y = "Nutritional Quality (Relative Lipid %)",
    color = "Order"
  )












#####Old code below

#add mouth part
library(tidyverse)

# 1. CREATE THE MOUTH LOOKUP TABLE
# Transcribed exactly from your updated image (image_2dfeff.png)
mouth_data <- tribble(
  ~Family,          ~Mouth_Type,
  "Geometridae",    "yes",
  "Crambidae",      "yes",
  "Noctuidae",      "yes",
  "Lasiocampidae",  "no",
  "Nematocera",     "yes",
  "Carabidae",      "yes",
  "Brachycera",     "yes",
  "Leptoceridae",   "yes",
  "Notodontidae",   "no",
  "Elateridae",     "yes",
  "Erebidae",       "mixed",
  "Pyralidae",      "yes",
  "Microlep",       "yes",
  "Gelechiidae",    "yes",
  "Caenidae",       "no",
  "Hydroptilidae",  "yes",
  "Limnephilidae",  "yes",
  "Heptageniidae",  "no",
  "Ephemeridae",    "no",
  "Cosmopterigidae","yes",
  "Heteroceridae",  "yes",
  "Staphylinidae",  "yes",
  "Tortricidae",    "yes",
  "Silphidae",      "yes",
  "Phryganeidae",   "no",
  "Helicopsychidae", "yes",
)

# 2. STANDARDIZE THE LABELS
# We map the raw "yes/no/mixed" to scientific terms.
mouth_data <- mouth_data %>%
  mutate(
    Feeding_Strategy = case_when(
      Mouth_Type == "no"    ~ "Capital (Non-Feeding)",
      Mouth_Type == "yes"   ~ "Income (Feeding)",
      Mouth_Type == "mixed" ~ "Mixed/Variable",
      TRUE ~ "Unknown"
    )
  )

# 3. MERGE WITH YOUR MAIN DATA
# We overwrite df_clean with the new version containing the 'Feeding_Strategy' column
df_clean <- df_clean %>%
  # Remove the old Feeding_Strategy column if it exists to avoid duplicates
  select(-any_of(c("Mouth_Type", "Feeding_Strategy"))) %>%
  left_join(mouth_data, by = "Family")

# Check for any families that didn't match (returns 0 rows if perfect)
print("=== FAMILIES MISSING MOUTH DATA ===")
print(df_clean %>% filter(is.na(Feeding_Strategy)) %>% distinct(Family))

# 4. RUN THE STATISTICAL TEST
# Does this biological trait predict fat content better than random chance?
mod_mouth <- lm(Logit_Lipid_Rel ~ Feeding_Strategy +Taxon_Group, data = df_clean)

print("=== ANOVA: FEEDING STRATEGY VS. FAT CONTENT ===")
print(anova(mod_mouth))
summary(mod_mouth)

# 5. VISUALIZE THE RESULTS
ggplot(df_clean %>% filter(!is.na(Feeding_Strategy)), 
       aes(x = Feeding_Strategy, y = Lipid_Rel_Pct, fill = Feeding_Strategy)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Capital (Non-Feeding)" = "#E69F00", 
                               "Income (Feeding)" = "#56B4E9",
                               "Mixed/Variable" = "grey70")) +
  theme_bw() +
  labs(
    title = "Does having a mouth predict fat content?",
    subtitle = "Comparison of Capital vs. Income vs. Mixed strategies",
    x = "Feeding Strategy",
    y = "Relative Lipid Content (%)"
  ) +
  theme(legend.position = "none")


####Mouth parts biplot
library(dplyr)

df_tradeoff <- df_clean %>%
  # Group by Family (so we get one dot per family)
  # But also include Order and Mouthparts so those columns stay in the result
  group_by(Family, Order, Mouth_Type) %>% 
  summarize(
    Mean_Mass = mean(Mass_Ind, na.rm = TRUE),
    Mean_Quality = mean(Lipid_Rel_Pct, na.rm = TRUE),
    Taxon_Group = first(Family), # Change label to Family name
    .groups = "drop"
  )

# Check the result: You should see many more rows now (one for each Family)
head(df_tradeoff)

plot_tradeoff_mouthparts <- ggplot(df_tradeoff, aes(x = Mean_Mass, y = Mean_Quality, color = Order)) +
  
  # A. The Grouping Circles (Ellipses)
  # We use 'Mouthparts' for the groups
  stat_ellipse(aes(group = Mouth_Type, fill = Mouth_Type), 
               geom = "polygon", alpha = 0.1, show.legend = TRUE) +
  
  # A2. Outlines for the circles
  stat_ellipse(aes(group = Mouth_Type, linetype = Mouth_Type), 
               geom = "path", color = "grey40", alpha = 0.6) +
  
  # B. The Points (Now representing Families)
  geom_point(size = 3, alpha = 0.8) +
  
  # C. The Labels (Family Names)
  geom_text_repel(aes(label = Taxon_Group), 
                  size = 3, 
                  box.padding = 0.4, 
                  max.overlaps = 20,
                  show.legend = FALSE) +
  
  # D. Scales & Colors
  scale_x_log10() + 
  
  # Colors for Mouthparts (Make sure these match your data: "Yes", "No", "Mixed")
  scale_fill_manual(values = c("yes" = "forestgreen", 
                               "no" = "firebrick", 
                               "mixed" = "orange")) +
  
  # E. Reference Line (Global Average)
  geom_hline(yintercept = mean(df_clean$Lipid_Rel_Pct, na.rm=TRUE), 
             linetype = "dashed", color = "grey50") +
  
  # F. Formatting
  theme_bw() +
  labs(
    title = "The Trade-off: Quantity vs. Quality (Family Averages)",
    subtitle = "Points = Families. Ellipses = Functional Mouthparts.",
    x = "Quantity: Individual Mass (g) [Log Scale]",
    y = "Quality: Relative Lipid Content (Proportion)",
    fill = "Mouthparts",
    linetype = "Mouthparts"
  ) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# Print
print(plot_tradeoff_mouthparts)

library(tidyverse)

# 1. FORCE CORRECT MONTH ORDER
# R defaults to alphabetical (April, August...), so we manually fix the levels.
df_clean <- df_clean %>%
  mutate(
    Month = factor(Month, levels = c("May", "June", "July", "August", "September"))
  ) %>%
  # Filter out any rows with missing months
  filter(!is.na(Month))

# 2. GENERATE THE PLOT
ggplot(df_clean, aes(x = Month, y = Lipid_Rel_Pct)) +
  
  # A. The Boxplot (Shows the spread/variance each month)
  geom_boxplot(aes(fill = Month), alpha = 0.6, outlier.shape = NA) +
  
  # B. The Raw Points (To see sample size)
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  
  # C. The Trend Line (Connecting the Means)
  # This helps answer "Does it increase?" by drawing a line through the average of each month.
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", size = 1, linetype = "dashed") +
  
  # D. Scales & Labels
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "YlOrRd") + # "Yellow-Orange-Red" palette looks like summer/autumn
  
  theme_bw() +
  labs(
    title = "Seasonal Trend: Does Fat Content Increase?",
    subtitle = "Diamonds/Dashed Line = Monthly Mean. Boxplots = Distribution.",
    x = "Month",
    y = "Relative Lipid Content (%)"
  ) +
  theme(legend.position = "none")

# 3. STATISTICAL CHECK
# We convert Month to a number (1, 2, 3...) to test for a linear trend (Slope).
mod_season <- lm(Logit_Lipid_Rel ~ as.numeric(Month), data = df_clean)

print("=== DOES MONTH PREDICT FAT? (Linear Trend) ===")
summary(mod_season)


