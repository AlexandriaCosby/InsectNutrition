# =========================================================
# 1. SETUP
# =========================================================
library(tidyverse)
library(ggpubr)
library(emmeans)

# Load Data
df <- read.csv("Insect_Nutrition_Prelim.csv", check.names = FALSE)

# =========================================================
# 2. DATA CALCULATION
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
# We use case_when to map the Order to the specific strategy

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
# Using base R shapiro.test() to avoid library conflicts

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

# =========================================================
# 4. TWO-WAY ANOVA (Order + Month)
# =========================================================

# --- Order Level Models ---
mod_ord_mass  <- lm(Log_Mass ~ Life_History + Order + Month, data = df_clean)
mod_ord_total <- lm(Log_Lipid_Total ~ Life_History + Order + Month, data = df_clean)
mod_ord_rel   <- lm(Log_Lipid_Rel ~ Life_History + Order + Month, data = df_clean)

# Print Results
print("=== ANOVA RESULTS (Two-Way: Taxonomy + Month) ===")
print("--- 1. Size (Log Mass) ---")
print(anova(mod_ord_mass))

print("--- 2. Total Lipid Content (Log) ---")
print(anova(mod_ord_total))

print("--- 3. Relative Lipid Content (Log) ---")
print(anova(mod_ord_rel))

# --- B. Create Filtered Data (Remove Singletons) ---
# We must count samples and remove families with N < 2
counts <- df_clean %>%
  group_by(Taxon_Group) %>%
  tally()

valid_groups <- counts %>%
  filter(n >= 2) %>%
  pull(Taxon_Group)

# Create the specific dataset for Family Analysis
df_fam_valid <- df_clean %>%
  filter(Taxon_Group %in% valid_groups)

print(paste("Original Families:", nrow(counts)))
print(paste("Valid Families (N>=2):", length(valid_groups)))

# --- C. Family Level Models (Run on FILTERED Data) ---
# Note: We name these '_valid' so they match your post-hoc code later
mod_fam_mass_valid  <- lm(Log_Mass ~ Life_History + Taxon_Group + Month, data = df_fam_valid)
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

# --- Order Level Plots ---
p_ord_mass <- create_boxplot(df_clean, "Order", "Mass_Ind", "Order", 
                             "Order: Individual Size", "Mass (g) - Log Scale")

p_ord_total <- create_boxplot(df_clean, "Order", "Lipid_Total", "Order", 
                              "Order: Total Lipid Content", "Total Lipid (g) - Log Scale")

p_ord_rel <- create_boxplot(df_clean, "Order", "Lipid_Rel_Pct", "Order", 
                            "Order: Relative Lipid Content", "Relative Lipid (Proportion) - Log Scale")

print(p_ord_mass)
print(p_ord_total)
print(p_ord_rel)

# --- Family Level Plots ---
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

# =========================================================
# 3. GENERATE THE CLEAN PLOTS
# =========================================================
# This removes " Slow", " Fast", or any other extra text, keeping only the first word.
df_clean <- df_clean %>%
  mutate(
    # Extract just the first word (The actual Order name)
    Order = word(Order, 1)
  )

# Verify it worked (Should just see standard names now)
print(unique(df_clean$Order))



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


# --- Helper Function for Sorting Family ---
create_sorted_boxplot <- function(data, x_col, y_col, fill_col, title_text, y_label) {
  ggplot(data, aes(x = reorder(.data[[x_col]], .data[[y_col]], median), 
                   y = .data[[y_col]], 
                   color = .data[[fill_col]])) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    
    # Facet by Order (scales="free_y" allows each facet to fit its own families)
    facet_grid(Order ~ ., scales = "free_y", space = "free_y") +
    
    scale_y_log10() +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.y = element_text(face = "bold", size = 9),
          strip.text = element_text(face = "bold")) +
    labs(title = title_text, subtitle = "Sorted by Median", y = y_label, x = "")
}

# --- Family Level Plots (Sorted) ---

p_fam_mass <- create_sorted_boxplot(df_fam_valid, "Taxon_Group", "Mass_Ind", "Order", 
                                    "Family: Individual Size", "Mass (g) - Log Scale")

p_fam_total <- create_sorted_boxplot(df_fam_valid, "Taxon_Group", "Lipid_Total", "Order", 
                                     "Family: Total Lipid Content", "Total Lipid (g) - Log Scale")

p_fam_rel <- create_sorted_boxplot(df_fam_valid, "Taxon_Group", "Lipid_Rel_Pct", "Order", 
                                   "Family: Relative Lipid Content", "Relative Lipid (Prop) - Log Scale")

print(p_fam_mass)
print(p_fam_total)
print(p_fam_rel)

#Family comparision heat map

library(tidyverse)
library(emmeans)

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
p_matrix_mass <- create_significance_heatmap(mod_fam_mass_valid, "Taxon_Group", "Mass Significance Matrix")
p_matrix_total <- create_significance_heatmap(mod_fam_total_valid, "Taxon_Group", "Total Lipid Significance Matrix")
p_matrix_rel <- create_significance_heatmap(mod_fam_rel_valid, "Taxon_Group", "Relative Lipid Content: Significance Matrix")

print(p_matrix_mass)
print(p_matrix_total)
print(p_matrix_rel)


####biplot
library(tidyverse)
library(ggrepel)

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

##biplot linear scale
library(tidyverse)
library(ggrepel)

# 1. PREPARE THE DATA (Must include Life_History)
df_tradeoff_linear <- df_clean %>%
  # Group by ALL factors to keep them in the summary
  group_by(Order, Taxon_Group, Life_History) %>%
  summarise(
    Mean_Mass = mean(Mass_Ind, na.rm = TRUE),         # Quantity
    Mean_Quality = mean(Lipid_Rel_Pct, na.rm = TRUE), # Quality
    Count = n(),
    .groups = "drop"
  ) %>%
  filter(Count >= 2) # Remove singletons for a cleaner chart

# 2. PLOT (LINEAR SCALE)
plot_tradeoff_linear <- ggplot(df_tradeoff_linear, aes(x = Mean_Mass, y = Mean_Quality)) +
  
  # A. The Grouping Circles (Ellipses)
  # Fill = Life History (Red/Blue)
  # Color = "grey50" (Neutral border so it doesn't fight with the Order colors)
  stat_ellipse(aes(fill = Life_History), 
               geom = "polygon", 
               color = "grey50",  # Neutral border color
               alpha = 0.1, 
               show.legend = TRUE) +
  
  # B. The Points (COLORED BY ORDER IS BACK!)
  geom_point(aes(color = Order), size = 4, alpha = 0.8) +
  
  # C. The Labels
  geom_text_repel(aes(label = Taxon_Group), 
                  size = 3.5, 
                  box.padding = 0.4,
                  point.padding = 0.3,
                  max.overlaps = 20, 
                  show.legend = FALSE) +
  
  # D. Scales 
  # Manual colors for Life History FILLS (Ellipses)
  scale_fill_manual(values = c("Fast" = "blue", "Slow" = "red", "Other" = "grey"), 
                    name = "Life History") +
  
  # Standard colors for Order POINTS (Let R pick them, or customize if you want)
  # We add a guide to ensure the legend looks nice
  guides(color = guide_legend(override.aes = list(linetype = 0, size = 4))) +
  
  # E. Reference Line (Global Average Fatness)
  geom_hline(yintercept = mean(df_clean$Lipid_Rel_Pct, na.rm=TRUE), 
             linetype = "dashed", color = "grey50") +
  
  # F. Formatting
  theme_bw() +
  labs(
    title = "Trade-off: Prey Size vs. Nutritional Quality",
    subtitle = "Linear Scale: Points colored by Order. Ellipses show Life History.",
    x = "Quantity: Individual Dry Mass (g)",
    y = "Quality: Relative Lipid Content (Proportion)"
  ) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

# Print
print(plot_tradeoff_linear)
######trying with a lm 
library(tidyverse)
library(sjPlot) # Excellent for visualizing LM coefficients (Install if needed: install.packages("sjPlot"))
library(sjmisc)

# =========================================================
# 1. RUN THE LINEAR MODELS (Functional Approach)
# =========================================================
# We use Life_History instead of Order to test the "Fast vs Slow" hypothesis
# Reference Level: R chooses one automatically (likely "Fast"), but we can set it.

# Model 1: Mass (Size)
lm_lh_mass <- lm(Log_Mass ~ Life_History + Month + Order, data = df_clean)

# Model 2: Total Lipid (Caloric Reward)
lm_lh_total <- lm(Log_Lipid_Total ~ Life_History + Month + Order, data = df_clean)

# Model 3: Relative Lipid (Quality)
lm_lh_rel <- lm(Log_Lipid_Rel ~ Life_History + Month + Order, data = df_clean)

# =========================================================
# 2. VIEW THE RESULTS (Coefficients)
# =========================================================
# The "Estimate" column is the most important.
# It tells you how much bigger/fattier the "Slow" group is compared to the "Fast" baseline.

print("=== LINEAR MODEL: SIZE (MASS) ===")
summary(lm_lh_mass)

print("=== LINEAR MODEL: TOTAL CALORIC REWARD ===")
summary(lm_lh_total)

print("=== LINEAR MODEL: NUTRITIONAL QUALITY (%) ===")
summary(lm_lh_rel)

# =========================================================
# 3. VISUALIZE THE EFFECT SIZES (Forest Plot)
# =========================================================
# This is the standard way to report LMs.
# Points to the RIGHT of the line = Positive effect (Bigger/Fattier than baseline)
# Points to the LEFT of the line = Negative effect (Smaller/Leaner than baseline)
# If the error bar crosses Zero, it is not significant.

plot_models(
  lm_lh_mass, lm_lh_total, lm_lh_rel,
  std.est = NULL,
  show.values = TRUE,
  show.p = TRUE,
  title = "Effect of Life History Strategy on Nutritional Metrics",
  m.labels = c("Mass", "Total Lipid", "Relative Lipid"),
  legend.title = "Metric"
) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(y = "Effect Size (Contrast vs. 'Fast' Life History)")


####More specific life spans

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
# We put 'Avg_Lifespan_Days' FIRST.
# This asks: "Does Lifespan explain the variance?"
# Then 'Family' asks: "Is there still variation between families of the SAME lifespan?"

# Model A: Size (Mass)
mod_test_mass <- lm(Log_Mass ~ Avg_Lifespan_Days + Family + Month, data = df_fam_valid)

# Model B: Nutritional Quality (% Lipid)
mod_test_rel <- lm(Log_Lipid_Rel ~ Avg_Lifespan_Days + Family + Month, data = df_fam_valid)

# =========================================================
# 3. ANOVA RESULTS (Type I Sum of Squares)
# =========================================================
print("=== ANOVA: MASS (Size) ===")
print(anova(mod_test_mass))
summary(mod_test_mass)

print("=== ANOVA: QUALITY (Relative Lipid) ===")
print(anova(mod_test_rel))

#Plot
library(tidyverse)
library(ggrepel)

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

# 2. GENERATE THE PLOT
plot_tradeoff_days <- ggplot(df_tradeoff_days, aes(x = Mean_Mass, y = Mean_Quality)) +
  
  # A. The Grouping Circles (Ellipses)
  # Grouped by our new "Lifespan Bins"
  stat_ellipse(aes(fill = Lifespan_Bin, color = Lifespan_Bin), 
               geom = "polygon", alpha = 0.1, show.legend = TRUE) +
  
  # B. The Points 
  # Colored by the Bin as well to make the pattern super clear
  geom_point(aes(color = Lifespan_Bin), size = 4, alpha = 0.8) +
  
  # C. The Labels (Family Names)
  geom_text_repel(aes(label = Family), 
                  size = 3, 
                  box.padding = 0.5, 
                  max.overlaps = 30,
                  show.legend = FALSE) +
  
  # D. Scales & Colors
  scale_x_log10() + 
  
  # Logical Color Scheme: Blue (Fast/Short) -> Orange (Medium) -> Red (Long/Slow)
  scale_fill_manual(values = c("Ephemeral (< 1 Wk)" = "#3498db",  # Blue
                               "Seasonal (1-8 Wks)" = "#f1c40f",  # Yellow/Orange
                               "Perennial (> 2 Mos)" = "#e74c3c"), # Red
                    name = "Adult Lifespan") +
  
  scale_color_manual(values = c("Ephemeral (< 1 Wk)" = "#3498db", 
                                "Seasonal (1-8 Wks)" = "#f39c12", # Darker Orange for dots
                                "Perennial (> 2 Mos)" = "#c0392b"), # Darker Red for dots
                     name = "Adult Lifespan") +
  
  # E. Reference Line
  geom_hline(yintercept = mean(df_clean$Lipid_Rel_Pct, na.rm=TRUE), 
             linetype = "dashed", color = "grey50") +
  
  # F. Formatting
  theme_bw() +
  labs(
    title = "Trade-off: Quantity vs. Quality (Grouped by Adult Lifespan)",
    subtitle = "Do longer-lived insects provide better nutrition?",
    x = "Quantity: Individual Mass (g) [Log Scale]",
    y = "Quality: Relative Lipid Content (Proportion)"
  ) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

# Print
print(plot_tradeoff_days)

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
  scale_x_log10() + 
  
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
mod_means <- lm(Mean_Lipid ~ Avg_Lifespan_Days, data = df_family_means)
print("=== REGRESSION RESULTS (FAMILY MEANS) ===")
summary(mod_means)


#####forgot to include size in my anovas, redoing this here. 

# =========================================================
# 1. RUN THE ANOVA (Does Sample_ID predict Quality?)
# =========================================================
# We filter only to ensure we have data for the model
df_stats <- df_clean %>% 
  filter(!is.na(Log_Lipid_Rel), !is.na(Sample_ID))

mod_anova <- lm(Log_Lipid_Rel ~ Sample_ID + Avg_Lifespan_Days + Month, data = df_stats)

print("=== ANOVA: PREDICTING QUALITY BY SAMPLE_ID ===")
print(anova(mod_anova))

# Optional: See who is significantly different
# print(TukeyHSD(aov(Log_Lipid_Rel ~ Sample_ID, data = df_stats)))

# =========================================================
# 2. RUN THE LINEAR REGRESSION (Lifespan Trend)
# =========================================================
# Does lifespan predict lipid content across all samples?
mod_lifespan <- lm(Log_Lipid_Rel ~ Avg_Lifespan_Days, data = df_clean)

print("=== REGRESSION: LIFESPAN PREDICTING QUALITY ===")
summary(mod_lifespan)

# =========================================================
# 3. GENERATE THE FIGURE
# =========================================================
# A. Summarize Data by Sample_ID (for the dots)
df_means <- df_clean %>%
  # Only exclude rows that are missing the visual coordinates
  filter(!is.na(Avg_Lifespan_Days), !is.na(Log_Lipid_Rel)) %>%
  group_by(Order, Sample_ID, Avg_Lifespan_Days) %>%
  summarise(
    Mean_Lipid = mean(Lipid_Rel_Pct, na.rm = TRUE),
    SD_Lipid = sd(Lipid_Rel_Pct, na.rm = TRUE),
    Count = n(),
    SE_Lipid = SD_Lipid / sqrt(Count),
    .groups = "drop"
  )

# B. Plot
ggplot(df_means, aes(x = Avg_Lifespan_Days, y = Mean_Lipid)) +
  
  # 1. Trend Line (Dashed)
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  
  # 2. Error Bars (Mean +/- SE)
  geom_errorbar(aes(ymin = Mean_Lipid - SE_Lipid, 
                    ymax = Mean_Lipid + SE_Lipid, 
                    color = Order), 
                width = 0.1, size = 0.6, alpha = 0.6) +
  
  # 3. Points (One dot per Sample_ID)
  geom_point(aes(color = Order), size = 4) +
  
  # 4. Labels (Using Sample_ID)
  geom_text_repel(aes(label = Sample_ID), 
                  size = 3.5, 
                  max.overlaps = 50,
                  show.legend = FALSE) +
  
  # 5. Scales & Theme
  scale_x_log10() + 
  theme_bw() +
  labs(
    title = "Does Adult Lifespan Predict Nutritional Quality?",
    subtitle = "Analysis grouped by Sample_ID.",
    x = "Adult Lifespan (Days) [Log Scale]",
    y = "Mean Nutritional Quality (Relative Lipid %)",
    color = "Order"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )


# This asks: "Do longer-lived Species Groups tend to be fattier?"
# This prevents families with huge sample sizes (like Beetles) from biasing the result.

# First, calculate the means (if you haven't already)
df_means_stats <- df_clean %>%
  filter(!is.na(Avg_Lifespan_Days), !is.na(Log_Lipid_Rel)) %>%
  group_by(Sample_ID) %>%
  summarise(
    Mean_Lipid_Log = mean(Log_Lipid_Rel, na.rm = TRUE),
    Lifespan = first(Avg_Lifespan_Days) # All samples in ID have same lifespan
  )

mod_lifespan_means <- lm(Mean_Lipid_Log ~ Lifespan, data = df_means_stats)

print("=== RESULTS: GROUP MEANS (Species Level) ===")
summary(mod_lifespan_means)

install.packages("pwr")
library(pwr)

# 1. INPUT THE VALUES FROM YOUR MODEL SUMMARY
r_squared <-  0.006023

# 2. CALCULATE EFFECT SIZE (Cohen's f2)
# Formula: R2 / (1 - R2)
f2_value <- r_squared / (1 - r_squared)

print(paste("Observed Effect Size (f2):", round(f2_value, 5)))
# Result will be ~0.0129 (Tiny)

# 3. RUN THE POWER TEST
# u = 1 (Because you have 1 predictor: Lifespan)
# power = 0.80 (Standard 80% chance to find the effect)
# sig.level = 0.05 (Standard P-value cutoff)
power_result <- pwr.f2.test(u = 1, 
                            f2 = f2_value, 
                            sig.level = 0.05, 
                            power = 0.80)

print("=== SAMPLES NEEDED ===")
print(power_result)

# 4. CALCULATE TOTAL N
# N = v (denominator df) + u (numerator df) + 1
total_n <- round(power_result$v + 1 + 1)

print(paste("Total Samples Required:", total_n))

###Do bigger bugs have more fat?

library(tidyverse)
library(ggrepel)

# =========================================================
# 1. RUN THE MODEL
# =========================================================
# Predictor: Log_Mass (Size)
# Response: Log_Lipid_Rel (Quality/Concentration)
mod_allometry <- lm(Log_Lipid_Total ~ Log_Mass + Month, data = df_clean)

print("=== REGRESSION RESULTS: DOES SIZE PREDICT FAT CONTENT? ===")
summary(mod_allometry)

# =========================================================
# 2. VISUALIZE THE RELATIONSHIP
# =========================================================
ggplot(df_clean, aes(x = Log_Mass, y = Log_Lipid_Total)) +
  
  # A. The Points (Colored by Order to see phylogenetic clusters)
  geom_point(aes(color = Order), size = 3, alpha = 0.7) +
  
  # B. The Regression Line (The biological trend)
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
  
  # C. Labels (Optional: Label only the outliers/extremes)
  # using sample_ID if available, or Family
  geom_text_repel(aes(label = Family), size = 3, max.overlaps = 10) +
  
  # D. Formatting
  theme_bw() +
  labs(
    title = "Allometry of Nutritional Quality",
    subtitle = "Does getting bigger (Mass) mean getting fattier (% Lipid)?",
    x = "Insect Size (Log Mass)",
    y = "Nutritional Quality (Log Relative Lipid)",
    color = "Order"
  )




library(scales) # Needed for the 'percent' or 'exp' formatting

ggplot(df_clean, aes(x = Log_Mass, y = Log_Lipid_Rel)) +
  
  # A. The Points
  geom_point(aes(color = Order), size = 3, alpha = 0.8) +
  
  # B. The Trend Line
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
  
  # C. Labels
  geom_text_repel(aes(label = Sample_ID), 
                  size = 3, 
                  max.overlaps = 20,
                  show.legend = FALSE) +
  
  # D. Formatting the Y-Axis to look like Real Numbers
  # This keeps the data logarithmic but labels it as the original decimal/percent
  scale_y_continuous(labels = function(x) paste0(round(exp(x) * 100, 1), "%")) +
  
  theme_bw() +
  labs(
    title = "Does Size Predict Nutritional Quality?",
    subtitle = "Higher on Y-axis = Fattier Insect.",
    x = "Insect Size (Log Mass)",
    # Now the label matches what the reader actually sees
    y = "Nutritional Quality (Estimated % Lipid)", 
    color = "Order"
  )

######Try a captial vs income approach instead of days alive 


library(tidyverse)

# 1. ASSIGN FEEDING GUILD
# I am using the EXACT list from your screenshot 'image_e32ff9.png'
df_guilds <- df_clean %>%
  mutate(
    Feeding_Guild = case_when(
      # === CAPITAL BREEDERS (Adults do not feed) ===
      # Moths
      Family == "Lasiocampidae" ~ "Capital",
      Family == "Notodontidae" ~ "Capital",
      
      # Caddisflies (Trichoptera)
      Family == "Leptoceridae" ~ "Capital",
      Family == "Limnephilidae" ~ "Capital",
      Family == "Phryganeidae" ~ "Capital",
      Family == "Hydroptilidae" ~ "Capital",
      
      # Mayflies (Ephemeroptera)
      Family == "Ephemeridae" ~ "Capital",
      Family == "Heptageniidae" ~ "Capital",
      Family == "Caenidae" ~ "Capital",
      
      # === INCOME BREEDERS (Adults feed) ===
      # Moths
      Family == "Noctuidae" ~ "Income",
      Family == "Erebidae" ~ "Income",
      Family == "Geometridae" ~ "Income",
      Family == "Pyralidae" ~ "Income",
      Family == "Crambidae" ~ "Income",
      Family == "Tortricidae" ~ "Income",
      Family == "Gelechiidae" ~ "Income",
      Family == "Cosmopterigidae" ~ "Income",
      
      # Beetles
      Family == "Carabidae" ~ "Income",
      Family == "Silphidae" ~ "Income",
      Family == "Staphylinidae" ~ "Income",
      Family == "Elateridae" ~ "Income",
      Family == "Heteroceridae" ~ "Income",
      
      # Flies
      Family == "Brachycera" ~ "Income",
      
      # === EXCLUDED / MIXED ===
      # Nematocera (Mosquitoes eat, Midges don't - too risky to group)
      # Microlep (Too distinct to group)
      TRUE ~ "Excluded"
    )
  ) %>%
  # Filter out anything that didn't fit the strict definitions
  filter(Feeding_Guild != "Excluded")

# 2. RUN THE ANOVA
mod_guild <- lm(Log_Lipid_Rel ~ Feeding_Guild + Month, data = df_guilds)

print("=== ANOVA RESULTS ===")
print(anova(mod_guild))

print("=== ESTIMATES (Direction of effect) ===")
summary(mod_guild)

# 3. VISUALIZE
ggplot(df_guilds, aes(x = Feeding_Guild, y = Lipid_Rel_Pct)) +
  geom_boxplot(aes(fill = Feeding_Guild), alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(
    title = "Does Feeding Strategy Predict Fat?",
    subtitle = "Capital (Non-feeding) vs. Income (Feeding)",
    y = "Lipid %",
    x = "Strategy"
  ) +
  theme(legend.position = "none")


library(tidyverse)
library(scales) # For percentage labels

# 1. ASSIGN FEEDING GUILD (Using your exact family list)
df_biplot <- df_clean %>%
  mutate(
    Feeding_Guild = case_when(
      # === CAPITAL BREEDERS (Non-Feeding) ===
      Family %in% c("Lasiocampidae", "Notodontidae", "Saturniidae", 
                    "Leptoceridae", "Limnephilidae", "Phryganeidae", "Hydroptilidae",
                    "Ephemeridae", "Heptageniidae", "Caenidae") ~ "Capital (Non-Feeding)",
      
      # === INCOME BREEDERS (Feeding) ===
      Family %in% c("Noctuidae", "Erebidae", "Geometridae", 
                    "Pyralidae", "Crambidae", "Tortricidae", "Gelechiidae", "Cosmopterigidae",
                    "Carabidae", "Silphidae", "Staphylinidae", "Elateridae", "Heteroceridae",
                    "Brachycera", "Muscidae", "Drosophilidae") ~ "Income (Feeding)",
      
      # === EXCLUDED ===
      TRUE ~ "Excluded"
    )
  ) %>%
  # Filter to only show the two groups we are comparing
  filter(Feeding_Guild != "Excluded")

# 2. GENERATE THE BIPLOT
ggplot(df_biplot, aes(x = Log_Mass, y = Log_Lipid_Rel, color = Feeding_Guild, fill = Feeding_Guild)) +
  
  # A. The Points
  geom_point(size = 3, alpha = 0.7) +
  
  # B. The "Circles" (Confidence Ellipses)
  # level = 0.95 draws a circle around 95% of the data for that group
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, show.legend = FALSE) +
  
  # C. Trend Lines (Optional - shows the slope difference)
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 0.8) +
  
  # D. Formatting Axes to be Human-Readable
  # We use the log data for the math, but label it as percentages
  scale_y_continuous(labels = function(x) paste0(round(exp(x) * 100, 1), "%")) +
  
  # E. Custom Colors (Orange for Capital, Blue for Income)
  scale_color_manual(values = c("Capital (Non-Feeding)" = "#E69F00", 
                                "Income (Feeding)" = "#56B4E9")) +
  scale_fill_manual(values = c("Capital (Non-Feeding)" = "#E69F00", 
                               "Income (Feeding)" = "#56B4E9")) +
  
  # F. Theme & Labels
  theme_bw() +
  labs(
    title = "Nutritional Niche: Capital vs. Income Breeders",
    subtitle = "Ellipses = 95% Confidence Intervals.",
    x = "Insect Size (Log Mass)",
    y = "Nutritional Quality (Estimated % Lipid)",
    color = "Strategy",
    fill = "Strategy"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )


# 1. RUN THE LINEAR MODEL
# We test if Family + Strategy + Month predicts Fat
mod_combined <- lm(Log_Lipid_Rel ~ Feeding_Guild + Family + Month, data = df_biplot)

print("=== MODEL SUMMARY ===")
summary(mod_combined)

# 2. RUN THE ANOVA (To see which variable carries the weight)
# The order matters here: We ask "Does Strategy matter?" first, then "Does Family matter after that?"
print("=== ANOVA RESULTS ===")
print(anova(mod_combined))
#garbage

#####Just trying to visualize shit

library(tidyverse)
library(scales)

# 1. PREPARE THE DATA
# Ensure we have the necessary columns and filter out missing IDs
df_vis <- df_clean %>%
  filter(!is.na(Sample_ID), !is.na(Lipid_Rel_Pct))

# 2. GENERATE THE RANKED PLOT
ggplot(df_vis, aes(x = reorder(Sample_ID, Lipid_Rel_Pct, FUN = median), y = Lipid_Rel_Pct)) +
  
  # A. The Boxplot (Shows the spread of data for each group)
  geom_boxplot(aes(fill = Order), alpha = 0.7, outlier.shape = NA) +
  
  # B. The Raw Points (Jittered so you can see sample size)
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  
  # C. Formatting
  scale_y_continuous(labels = scales::percent) + # Shows 5%, 10% instead of 0.05, 0.1
  coord_flip() + # FLIPS the plot so names are readable on the left
  theme_bw() +
  
  # D. Labels
  labs(
    title = "Lipid Content by Sample ID (Ranked)",
    subtitle = "Ordered from Leanest (Bottom) to Fattiest (Top)",
    x = "Sample ID",
    y = "Relative Lipid Content (%)",
    fill = "Order" # Legend for colors
  ) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"), # Make names easy to read
    legend.position = "right"
  )

