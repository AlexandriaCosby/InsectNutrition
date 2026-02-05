# =========================================================
# 1. SETUP
# =========================================================
library(tidyverse)
library(ggpubr)
library(emmeans)
library(ggrepel)

# Load Data
df <- read.csv("Insect_Nutrition_Prelim.csv", check.names = FALSE)

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
    Log_Lipid_Rel = log(Lipid_Rel_Pct)
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
print(paste("Log Data p:", format(shapiro.test(df_clean$Log_Lipid_Rel)$p.value, digits=3)))

##Confirms that log transform makes data normal 

# =========================================================
# 4. TWO-WAY ANOVA (Order + Month)
# =========================================================
# --- Order Level Models ---
mod_ord_mass  <- lm(Log_Mass ~ Order + Month, data = df_clean)
mod_ord_total <- lm(Log_Lipid_Total ~ Order + Month, data = df_clean)
mod_ord_rel   <- lm(Log_Lipid_Rel ~ Order + Month, data = df_clean)

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
mod_fam_mass_valid  <- lm(Log_Mass ~ Life_Hisotry + Taxon_Group + Month, data = df_fam_valid)
mod_fam_total_valid <- lm(Log_Lipid_Total ~ Life_History + Taxon_Group + Month, data = df_fam_valid)
mod_fam_rel_valid   <- lm(Log_Lipid_Rel ~ Life_History + Taxon_Group + Month, data = df_fam_valid)

print("=== ANOVA RESULTS (Family Level: Valid N > 1) ===")
print("--- 2. Family Level ---")
print(anova(mod_fam_mass_valid))
print(anova(mod_fam_total_valid))
print(anova(mod_fam_rel_valid))

# =========================================================
# 5. SCIENTIFIC PLOTTING (Boxplots)
# =========================================================
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
p_fam_mass <- create_boxplot(df_clean, "Taxon_Group", "Mass_Ind", "Order", 
                             "Family: Individual Size", "Mass (g) - Log Scale") + 
  facet_grid(Order ~ ., scales = "free_y", space = "free_y")

p_fam_total <- create_boxplot(df_clean, "Taxon_Group", "Lipid_Total", "Order", 
                              "Family: Total Lipid Content", "Total Lipid (g) - Log Scale") + 
  facet_grid(Order ~ ., scales = "free_y", space = "free_y")

p_fam_rel <- create_boxplot(df_clean, "Taxon_Group", "Lipid_Rel_Pct", "Order", 
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
mod_fam_rel_valid1   <- lm(Log_Lipid_Rel ~ Taxon_Group + Month, data = df_fam_valid)

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

# 1. PREPARE THE DATA (Now including Life_History!)
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
  stat_ellipse(aes(group = Life_History, fill = Life_History), 
               geom = "polygon", alpha = 0.1, show.legend = TRUE) +
  
  # A2. Optional: Dashed lines around the circles to make them pop
  stat_ellipse(aes(group = Life_History, linetype = Life_History), 
               geom = "path", color = "grey40", alpha = 0.6) +
  
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
  # We need to manually set fill colors for the ellipses (Fast vs Slow)
  scale_fill_manual(values = c("Fast" = "blue", "Slow" = "red", "Other" = "grey")) +
  
  # E. Reference Line
  geom_hline(yintercept = mean(df_clean$Lipid_Rel_Pct, na.rm=TRUE), 
             linetype = "dashed", color = "grey50") +
  
  # F. Formatting
  theme_bw() +
  labs(
    title = "The Trade-off: Quantity vs. Quality (Grouped by Life History)",
    subtitle = "Ellipses show the clustering of Fast vs. Slow life histories.",
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

# =========================================================
# ADD LIFESPAN DATA (Using 'Family' Column Only)
# =========================================================

df_clean <- df_clean %>%
  mutate(
    Avg_Lifespan_Days = case_when(
      
      # --- FLIES (Diptera) ---
      # Now looking directly for these names in the Family column
      Family == "Brachycera" ~ 21,
      Family == "Nematocera " ~ 14,
      Family %in% c("Muscidae", "Drosophilidae") ~ 21, 
      Family %in% c("Tipulidae", "Chironomidae", "Culicidae") ~ 14,
      
      # --- MOTHS (Lepidoptera) ---
      Family == "Pyralidae" ~ 30,
      Family == "Noctuidae" ~ 21,
      Family %in% c("Tortricidae", "Erebidae") ~ 15,
      Family %in% c("Crambidae", "Gelechiidae", "Cosmopterigidae") ~ 14,
      Family == "Microlep" ~ 14, 
      Family %in% c("Geometridae", "Lasiocampidae", "Notodontidae") ~ 7,
      
      # --- BEETLES (Coleoptera) ---
      Family %in% c("Carabidae", "Silphidae") ~ 365,
      Family %in% c("Staphylinidae", "Elateridae") ~ 40,
      Family == "Heteroceridae" ~ 30,
      
      # --- AQUATIC INSECTS ---
      Family == "Limnephilidae" ~ 60,
      Family %in% c("Leptoceridae", "Phryganeidae") ~ 30,
      Family == "Hydroptilidae" ~ 14,
      Family %in% c("Ephemeridae", "Heptageniidae") ~ 2,
      Family == "Caenidae" ~ 0.5,
      
      # --- DEFAULT ---
      TRUE ~ NA_real_
    )
  )

# Verify results
print(head(df_clean %>% select(Family, Avg_Lifespan_Days)))

# =========================================================
# 2. RUN THE COMBINED MODELS
# =========================================================

# Model A: Size (Mass)
mod_test_mass <- lm(Log_Mass ~ Avg_Lifespan_Days + Taxon_Group + Month, data = df_fam_valid)

# Model B: Nutritional Quality (% Lipid)
mod_test_rel <- lm(Log_Lipid_Rel ~ Avg_Lifespan_Days + Taxon_Group + Month, data = df_fam_valid)

# =========================================================
# 3. ANOVA RESULTS (Type I Sum of Squares)
# =========================================================
print("=== ANOVA: MASS (Size) ===")
print(anova(mod_test_mass))

print("=== ANOVA: QUALITY (Relative Lipid) ===")
print(anova(mod_test_rel))

# 1. PREPARE THE DATA & CREATE BINS
df_tradeoff_days <- df_clean %>%
  # Filter out rows with missing lifespan data
  filter(!is.na(Avg_Lifespan_Days)) %>%
  
  # Create Logical Bins for the Ellipses
  mutate(Lifespan_Bin = case_when(
    Avg_Lifespan_Days <= 7  ~ "Ephemeral (< 1 Wk)",
    Avg_Lifespan_Days <= 60 ~ "Seasonal (1-8 Wks)",
    Avg_Lifespan_Days > 60  ~ "Perennial (> 2 Mos)"
  )) %>%
  
  # Force the order of the bins for the legend
  mutate(Lifespan_Bin = factor(Lifespan_Bin, 
                               levels = c("Ephemeral (< 1 Wk)", "Seasonal (1-8 Wks)", "Perennial (> 2 Mos)"))) %>%
  
  # Group and Summarize
  group_by(Order, Family, Lifespan_Bin) %>%
  summarise(
    Mean_Mass = mean(Mass_Ind, na.rm = TRUE),
    Mean_Quality = mean(Lipid_Rel_Pct, na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  ) %>%
  filter(Count >= 2)

# 2. GENERATE THE PLOT (With Log-Transformed Y-Axis)
plot_tradeoff_days <- ggplot(df_tradeoff_days, aes(x = Mean_Mass, y = Mean_Quality)) +
  
  # A. The Grouping Circles (Ellipses)
  stat_ellipse(aes(fill = Lifespan_Bin, color = Lifespan_Bin), 
               geom = "polygon", alpha = 0.1, show.legend = TRUE) +
  
  # B. The Points 
  geom_point(aes(color = Lifespan_Bin), size = 4, alpha = 0.8) +
  
  # C. The Labels (Family Names)
  geom_text_repel(aes(label = Family), 
                  size = 3, 
                  box.padding = 0.5, 
                  max.overlaps = 30,
                  show.legend = FALSE) +
  
  # D. SCALES (Both Axes are now Log-10)
  scale_x_log10() + 
  scale_y_log10() +  # <--- THIS WAS THE MISSING PIECE
  
  # E. Colors (Blue -> Yellow -> Red)
  scale_fill_manual(values = c("Ephemeral (< 1 Wk)" = "#3498db", 
                               "Seasonal (1-8 Wks)" = "#f1c40f", 
                               "Perennial (> 2 Mos)" = "#e74c3c"),
                    name = "Adult Lifespan") +
  
  scale_color_manual(values = c("Ephemeral (< 1 Wk)" = "#3498db", 
                                "Seasonal (1-8 Wks)" = "#f39c12", 
                                "Perennial (> 2 Mos)" = "#c0392b"),
                     name = "Adult Lifespan") +
  
  # F. Reference Line (Global Mean)
  # Since we are using scale_y_log10, we can still plot the raw mean, 
  # and ggplot will position it correctly on the log axis.
  geom_hline(yintercept = mean(df_clean$Lipid_Rel_Pct, na.rm=TRUE), 
             linetype = "dashed", color = "grey50") +
  
  # G. Formatting
  theme_bw() +
  labs(
    title = "Trade-off: Quantity vs. Quality (Grouped by Adult Lifespan)",
    subtitle = "Both axes are Log-Transformed to show allometric scaling.",
    x = "Quantity: Individual Mass (g) [Log Scale]",
    y = "Quality: Relative Lipid Content [Log Scale]"
  ) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

# Print
print(plot_tradeoff_days)

# =========================================================
# LM Comparing lifespan and lipid content
# =========================================================

# 1. SUMMARIZE DATA BY FAMILY (Same as before)
df_family_means <- df_fam_valid %>%
  group_by(Order, Taxon_Group, Avg_Lifespan_Days) %>%
  summarise(
    Mean_Lipid = mean(Lipid_Rel_Pct, na.rm = TRUE),
    SD_Lipid = sd(Lipid_Rel_Pct, na.rm = TRUE),
    Count = n(),
    SE_Lipid = SD_Lipid / sqrt(Count), 
    .groups = "drop"
  ) %>%
  filter(!is.na(Avg_Lifespan_Days)) %>%
  filter(!is.na(Mean_Lipid))

# 2. GENERATE THE PLOT (No Shading)
plot_lifespan_means <- ggplot(df_family_means, aes(x = Avg_Lifespan_Days, y = Mean_Lipid)) +
  
  # A. The Regression Line (Trend ONLY, no shading)
  # se = FALSE removes the gray ribbon
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  
  # B. Error Bars (Mean +/- SE)
  geom_errorbar(aes(ymin = Mean_Lipid - SE_Lipid, 
                    ymax = Mean_Lipid + SE_Lipid, 
                    color = Order), 
                width = 0.1, size = 0.6, alpha = 0.6) +
  
  # C. The Points (Mean Values)
  geom_point(aes(color = Order), size = 4) +
  
  # D. Labels
  geom_text_repel(aes(label = Taxon_Group), 
                  size = 3.5, 
                  box.padding = 0.5,
                  max.overlaps = 20,
                  show.legend = FALSE) +
  
  # E. Scales
  scale_y_log10() + 
  scale_x_log10()+
  # F. Formatting
  theme_bw() +
  labs(
    title = "Does Adult Lifespan Predict Nutritional Quality?",
    subtitle = "Points = Family Means. Dashed Line = Linear Trend (Family Level).",
    x = "Average Adult Lifespan (Days) [Log Scale]",
    y = "Mean Nutritional Quality (Relative Lipid %)",
    color = "Order"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Print
print(plot_lifespan_means)

# 3. STAT CHECK (Run the regression on these MEANS)
# Does the trend hold up when we treat the Family as the unit of analysis?
mod_means <- lm(log(Mean_Lipid) ~ log(Avg_Lifespan_Days), data = df_family_means)
print("=== REGRESSION RESULTS (FAMILY MEANS) ===")
summary(mod_means)


install.packages("pwr")
library(pwr)

# 1. INPUT THE VALUES FROM YOUR MODEL SUMMARY
r_squared <-  0.049

# 2. CALCULATE EFFECT SIZE (Cohen's f2)
f2_value <- r_squared / (1 - r_squared)

print(paste("Observed Effect Size (f2):", round(f2_value, 5)))

# 3. RUN THE POWER TEST
power_result <- pwr.f2.test(u = 1, 
                            f2 = f2_value, 
                            sig.level = 0.05, 
                            power = 0.80)

print("=== SAMPLES NEEDED ===")
print(power_result)

# 4. CALCULATE TOTAL N
total_n <- round(power_result$v + 1 + 1)
print(paste("Total Samples Required:", total_n))

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
  "Nematocera",     "mixed",
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
  "Phryganeidae",   "no"
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
mod_mouth <- lm(Log_Lipid_Rel ~ Feeding_Strategy +Taxon_Group, data = df_clean)

print("=== ANOVA: FEEDING STRATEGY VS. FAT CONTENT ===")
print(anova(mod_mouth))

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
