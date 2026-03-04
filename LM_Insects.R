# =========================================================
# 1. SETUP
# =========================================================
library(tidyverse)
library(ggpubr)
library(emmeans)
library(ggrepel)
library(ggplot2)
library(ggsignif)
library(patchwork)
library(lme4)
library(glmmTMB)  # For Generalized Linear Mixed Models (GLMMs), specifically Beta regression
library(DHARMa)
library(car)
library(performance)
library(ggsignif)

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
    Log_Lipid_Rel = log(Lipid_Rel_Pct)
  )


# =========================================================
# 3. NORMALITY TESTING (Raw vs. Log)
# =========================================================
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
qqnorm(df_clean$Log_Lipid_Rel, main="Q-Q Plot: Log Relative Lipid")
qqline(df_clean$Log_Lipid_Rel, col="red", lwd=2)

# Reset plot layout to default
par(mfrow=c(1,1))

##Moving forward with logged data

# =========================================================
# Add Mouth column
# =========================================================
# 1. CREATE THE MOUTH LOOKUP TABLE
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
  "Microlep",       "mixed",#to be safe
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
  "Phryganeidae",   "yes", 
  "Helicopsychidae", "yes",
)

# 2. STANDARDIZE THE LABELS
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
df_clean <- df_clean %>%
  select(-any_of(c("Mouth_Type", "Feeding_Strategy"))) %>%
  left_join(mouth_data, by = "Family")

# =========================================================
# Add habitat column
# =========================================================
habitat_data <- tribble(
  ~Family,           ~Habitat_Type,
  "Geometridae",     "terrestrial",
  "Crambidae",       "mixed",        
  "Noctuidae",       "terrestrial",
  "Lasiocampidae",   "terrestrial",
  "Nematocera",      "semi-aquatic", 
  "Carabidae",       "terrestrial",  
  "Brachycera",      "mixed",        
  "Leptoceridae",    "semi-aquatic", 
  "Notodontidae",    "terrestrial",
  "Elateridae",      "terrestrial",
  "Erebidae",        "terrestrial",
  "Pyralidae",       "terrestrial",
  "Microlep",        "terrestrial", 
  "Gelechiidae",     "terrestrial",
  "Caenidae",        "semi-aquatic", 
  "Hydroptilidae",   "semi-aquatic", 
  "Limnephilidae",   "semi-aquatic", 
  "Heptageniidae",   "semi-aquatic", 
  "Ephemeridae",     "semi-aquatic", 
  "Cosmopterigidae", "terrestrial",
  "Heteroceridae",   "semi-aquatic", 
  "Staphylinidae",   "terrestrial",  
  "Tortricidae",     "terrestrial",
  "Silphidae",       "terrestrial",
  "Phryganeidae",    "semi-aquatic", 
  "Helicopsychidae", "semi-aquatic"  
)


df_clean <- df_clean %>%
  select(-any_of(c("Habitat_Type"))) %>%
  left_join(habitat_data, by = "Family")


# =========================================================
# Seeing if it's important to retain sample month
# =========================================================
##### Mass
# Model without Month
mod_no_month_mass <- lmer(Log_Mass ~ Order + Mouth_Type + (1 | Family), data = df_clean, REML = FALSE)
# Model with Month
mod_with_month_mass <- lmer(Log_Mass ~ Order + Mouth_Type + (1 | Family) + (1 | Month), data = df_clean, REML = FALSE)
# Compare them
anova(mod_no_month_mass, mod_with_month_mass)

##### Lipid (unsure if I'm going to keep lipid)
# Model without Month
mod_no_month_lipid <- lmer(Log_Lipid_Total ~ Order + Mouth_Type + (1 | Family), data = df_clean, REML = FALSE)
# Model with Month
mod_with_month_lipid <- lmer(Log_Lipid_Total ~ Order + Mouth_Type + (1 | Family) + (1 | Month), data = df_clean, REML = FALSE)
# Compare them
anova(mod_no_month_lipid, mod_with_month_lipid)

##### Rel Lipid
# Model without Month
mod_no_month <- lmer(Log_Lipid_Rel ~ Order + Mouth_Type + (1 | Family), data = df_clean, REML = FALSE)
# Model with Month
mod_with_month <- lmer(Log_Lipid_Rel ~ Order + Mouth_Type + (1 | Family) + (1 | Month), data = df_clean, REML = FALSE)
# Compare them
anova(mod_no_month, mod_with_month)

## --> Don't bother including month

# =========================================================
# Set up boxplot criteria for figures
# =========================================================
create_phd_boxplot <- function(data, x_var, y_var, fill_var, title_text, y_label, model = NULL) {
  
  # 1. Force the x-variable to be a factor and capture its order
  data[[x_var]] <- as.factor(data[[x_var]])
  levs <- levels(data[[x_var]])
  
  # 2. Base Plot
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[fill_var]])) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.3, color = "black") +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, color = "darkgray") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
          axis.title.y = element_text(face = "bold")) +
    labs(title = title_text, y = y_label, x = "")
  
  # 3. Bracket Logic
  if (!is.null(model)) {
    emm <- emmeans(model, specs = x_var)
    pairs_data <- as.data.frame(pairs(emm, adjust = "tukey"))
    sig_pairs <- pairs_data[pairs_data$p.value < 0.05, ]
    
    if (nrow(sig_pairs) > 0) {
      xmin_vec <- numeric(nrow(sig_pairs))
      xmax_vec <- numeric(nrow(sig_pairs))
      
      for(i in 1:nrow(sig_pairs)) {
        # Split the contrast (e.g., "mixed - terrestrial")
        pair <- trimws(strsplit(as.character(sig_pairs$contrast[i]), " - ")[[1]])
        
        # FUZZY MATCH: find which axis level is IN the model's pair strings
        # This handles hidden spaces or "Habitat_Type" prefixes automatically
        xmin_match <- which(sapply(levs, function(l) grepl(tolower(trimws(l)), tolower(pair[1]))))
        xmax_match <- which(sapply(levs, function(l) grepl(tolower(trimws(l)), tolower(pair[2]))))
        
        # Only assign if a match was found to prevent the "length zero" error
        if(length(xmin_match) > 0) xmin_vec[i] <- xmin_match[1]
        if(length(xmax_match) > 0) xmax_vec[i] <- xmax_match[1]
      }
      
      # Filter out any pairs that failed to match
      valid <- xmin_vec > 0 & xmax_vec > 0
      
      if(any(valid)) {
        max_y <- max(data[[y_var]], na.rm = TRUE)
        y_range <- max_y - min(data[[y_var]], na.rm = TRUE)
        ladder <- seq(from = max_y + (y_range * 0.1), 
                      by = y_range * 0.15, 
                      length.out = sum(valid))
        
        p <- p + geom_signif(
          xmin = xmin_vec[valid], 
          xmax = xmax_vec[valid],
          y_position = ladder,
          annotations = rep("*", sum(valid)),
          tip_length = 0.02, 
          color = "black", 
          textsize = 7, 
          vjust = 0.5
        ) +
          # Ensures brackets aren't clipped off the top of the plot
          scale_y_continuous(expand = expansion(mult = c(0.05, 0.4)))
      }
    }
  }
  return(p)
}


# =========================================================
# Looking at factors that best predict Mass
# =========================================================
#1. Taxonomy
# Nested model: Family is nested within Order (Order / Family)
mod_tax_mass <- lm(Log_Mass ~ Order / Sample_ID, data = df_clean)
Anova(mod_tax_mass, type = "II")

# 2. Functional Model
mod_func_mass <- lm(Log_Mass ~ Order + Mouth_Type, data = df_clean)
Anova(mod_func_mass, type = "II")

# 3. Environmental Model
mod_habitat_mass <- lm(Log_Mass ~ Order + Habitat_Type, data = df_clean)
Anova(mod_habitat_mass, type = "II")


#####Visualize mass

#pulling out order differences in mass
order_means_mass <- emmeans(mod_tax_mass, ~ Order)
order_pairs_mass <- contrast(order_means_mass, method = "pairwise")
order_pairs_df_mass <- as.data.frame(order_pairs_mass)
significant_orders_mass <- subset(order_pairs_df_mass, p.value < 0.05)
print(significant_orders_mass[, c("contrast", "estimate", "p.value")])

#pulling out functional differences in mass
func_means_mass <- emmeans(mod_func_mass, ~ Mouth_Type)
func_pairs_mass <- contrast(func_means_mass, method = "pairwise")
func_pairs_df_mass <- as.data.frame(func_pairs_mass)
significant_func_mass <- subset(func_pairs_df_mass, p.value < 0.05)
print(significant_func_mass[, c("contrast", "estimate", "p.value")])

#pulling out habitat differences in mass
hab_means_mass <- emmeans(mod_habitat_mass, ~ Habitat_Type )
hab_pairs_mass <- contrast(hab_means_mass, method = "pairwise")
hab_pairs_df_mass <- as.data.frame(hab_pairs_mass)
significant_hab_mass <- subset(hab_pairs_df_mass, p.value < 0.05)
print(significant_hab_mass[, c("contrast", "estimate", "p.value")])

##Visualize order differences in mass
#Order (having to do this because the nesting doesn't yield brackets)
# 1. Define the 5 significant pairs 
mass_comparisons <- list(
  c("Coleoptera", "Diptera"),
  c("Coleoptera", "Ephemeroptera"),
  c("Coleoptera", "Trichoptera"),
  c("Diptera", "Ephemeroptera"),
  c("Diptera", "Lepidoptera"),
  c("Diptera", "Trichoptera"),
  c("Ephemeroptera", "Lepidoptera"),
  c("Lepidoptera", "Trichoptera")
)

# 2. Create the final plot
plot_mass_nested <- ggplot(df_clean, aes(x = Order, y = Log_Mass, fill = Order)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.3, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, color = "darkgray") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
        axis.title.y = element_text(face = "bold")) +
  labs(title = "Insect Biomass by Taxonomic Order", 
       y = "Log Total Mass (mg)", 
       x = "") +
  
  # 3. Add the significance brackets manually
  geom_signif(
    comparisons = mass_comparisons,
    annotations = rep("*", length(mass_comparisons)),
    step_increase = 0.12, 
    color = "black",
    textsize = 7,
    vjust = 0.5
  ) +
  
  # 4. Expand the Y-axis so the brackets don't get cut off
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.5)))

# View the plot
plot_mass_nested

#Functional
plot_func_mass <- create_phd_boxplot(
  data = df_clean, 
  x_var = "Mouth_Type", 
  y_var = "Log_Mass", 
  fill_var = "Mouth_Type", 
  title_text = "Insect Biomass by Mouth Type", 
  y_label = "Log Total Mass (mg)", 
  model = mod_func_mass # This triggers the brackets!
)
plot_func_mass

#Habitat
plot_hab_mass <- create_phd_boxplot(
  data = df_clean, 
  x_var = "Habitat_Type", 
  y_var = "Log_Mass", 
  fill_var = "Habitat_Type", 
  title_text = "Insect Biomass by Habitat Type", 
  y_label = "Log Total Mass (mg)", 
  model = mod_habitat_mass
)
plot_hab_mass

# =========================================================
# Which best describes size drivers of insect orders?
# =========================================================
compare_models <- function(tax_mod, func_mod, hab_mod, label) {
  cat("\n--- Comparison for:", label, "---\n")
  print(compare_performance(tax_mod, func_mod, rank = TRUE))
}

# Comparison for Mass
compare_performance(mod_tax_mass, mod_func_mass, mod_habitat_mass, rank = TRUE)

# =========================================================
# Looking at factors that best predict Relative Lipid content
# =========================================================
#1. Taxonomy
# Nested model: Family is nested within Order (Order / Family)
mod_tax_rel_lipid <- lm(Log_Lipid_Rel ~ Order / Sample_ID, data = df_clean)
Anova(mod_tax_rel_lipid, type = "II")

# 2. Functional Model
mod_func_rel_lipid <- lm(Log_Lipid_Rel ~ Order + Mouth_Type, data = df_clean)
Anova(mod_func_rel_lipid, type = "II")

# 3. Environmental Model
mod_hab_rel_lipid <- lm(Log_Lipid_Rel ~ Order + Habitat_Type, data = df_clean)
Anova(mod_hab_rel_lipid, type = "II")


#####Visualize, pulling our order differences
#Order
order_means_rel_lipid <- emmeans(mod_tax_rel_lipid, ~ Order)
order_pairs_rel_lipid <- contrast(order_means_rel_lipid, method = "pairwise")
order_pairs_df_rel_lipid <- as.data.frame(order_pairs_rel_lipid)
significant_orders_rel_lipid <- subset(order_pairs_df_rel_lipid, p.value < 0.05)
print(significant_orders_rel_lipid[, c("contrast", "estimate", "p.value")])

#Functional
func_means_rel_lipid <- emmeans(mod_func_rel_lipid, ~ Mouth_Type)
func_pairs_rel_lipid <- contrast(func_means_rel_lipid, method = "pairwise")
func_pairs_df_rel_lipid <- as.data.frame(func_pairs_rel_lipid)
significant_func_rel_lipid <- subset(func_pairs_df_rel_lipid, p.value < 0.05)
print(significant_func_rel_lipid[, c("contrast", "estimate", "p.value")])

#Habitat
hab_means_rel_lipid <- emmeans(mod_hab_rel_lipid, ~ Habitat_Type)
hab_pairs_rel_lipid <- contrast(hab_means_rel_lipid, method = "pairwise")
hab_pairs_df_rel_lipid <- as.data.frame(hab_pairs_rel_lipid)
significant_hab_rel_lipid <- subset(hab_pairs_df_rel_lipid, p.value < 0.05)
print(significant_hab_rel_lipid[, c("contrast", "estimate", "p.value")])

##Visualize order differences in relative lipid
#Order (have to manually add brackets because it doesn't like the nesting)
manual_comparisons <- list(
  c("Coleoptera", "Diptera"),
  c("Coleoptera", "Ephemeroptera"),
  c("Diptera", "Lepidoptera"),
  c("Diptera", "Trichoptera"),
  c("Ephemeroptera", "Lepidoptera"),
  c("Ephemeroptera", "Trichoptera")
)

# 2. Create the plot using your nested model data
plot_manual_brackets <- ggplot(df_clean, aes(x = Order, y = Log_Lipid_Rel, fill = Order)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.3, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, color = "darkgray") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
        axis.title.y = element_text(face = "bold")) +
  labs(title = "Nutritional Quality by Taxonomic Order", 
       y = "Log Relative Lipid Content", 
       x = "") +
  
  # 3. Manually add the brackets
  geom_signif(
    comparisons = manual_comparisons,
    annotations = rep("*", length(manual_comparisons)), # Adds one star per bracket
    step_increase = 0.12, # Spreads the 6 brackets out vertically so they don't overlap
    color = "black",
    textsize = 7,
    vjust = 0.5
  ) +
  
  # 4. Expand the Y-axis so the 6-rung ladder fits on your poster
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.6)))

plot_manual_brackets

#Functional
plot_func_mass <- create_phd_boxplot(
  data = df_clean, 
  x_var = "Mouth_Type", 
  y_var = "Log_Mass", 
  fill_var = "Mouth_Type", 
  title_text = "Insect Biomass by Mouth Type", 
  y_label = "Log Total Mass (mg)", 
  model = mod_func_mass # This triggers the brackets!
)
plot_func_mass

#Habitat
plot_hab_mass <- create_phd_boxplot(
  data = df_clean, 
  x_var = "Habitat_Type", 
  y_var = "Log_Mass", 
  fill_var = "Habitat_Type", 
  title_text = "Insect Biomass by Habitat Type", 
  y_label = "Log Total Mass (mg)", 
  model = mod_habitat_mass
)
plot_hab_mass

# =========================================================
# Which best describes size drivers of insect orders?
# =========================================================
compare_models_rel_lipid <- function(mod_tax_rel_lipid, mod_func_rel_lipid, mod_hab_rel_lipid, label) {
  cat("\n--- Comparison for:", label, "---\n")
  print(compare_performance(mod_tax_rel_lipid, mod_func_rel_lipid, rank = TRUE))
}

# Comparison for Mass
compare_performance(mod_tax_rel_lipid, mod_func_rel_lipid, mod_hab_rel_lipid, rank = TRUE)


# =========================================================
# Creating plots to comapre groups
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

# 1. Mass (Size)
p_fam_mass <- create_boxplot_clean(
  data = df_clean, 
  x_var = "Sample_ID", 
  y_var = "Mass_Ind", 
  fill_var = "Order",
  title_text = "A) Individual Size", 
  y_label = "Mass (g) - Log Scale"
) + 
  facet_grid(Order ~ ., scales = "free_y", space = "free_y")

# 2. Relative Lipid
p_fam_rel <- create_boxplot_clean(
  data = df_clean, 
  x_var = "Sample_ID", 
  y_var = "Lipid_Rel_Pct", 
  fill_var = "Order",
  title_text = "B) Relative Lipid Content", 
  y_label = "Relative Lipid (Proportion) - Log Scale"
) + 
  facet_grid(Order ~ ., scales = "free_y", space = "free_y")

# Print to verify
print(p_fam_mass)
print(p_fam_rel)

multi_panel_fam= p_fam_mass | p_fam_rel
print(multi_panel_fam)


# =========================================================
# Heatmap comparing families
# =========================================================
# Same fit as mod_tax_mass, but "flat" for the heatmap
mod_heatmap_mass <- lm(Log_Mass ~ Sample_ID, data = df_clean)

# Same fit as mod_tax_rel_lipid, but "flat" for the heatmap
mod_heatmap_rel <- lm(Log_Lipid_Rel ~ Sample_ID, data = df_clean)

create_significance_heatmap <- function(model, factor_name, title_text) {
  
  # 1. Get Pairwise Comparisons
  # No nesting argument needed here because Sample_ID is now the main effect
  emm <- emmeans(model, specs = factor_name)
  pairs_data <- as.data.frame(pairs(emm, adjust = "tukey")) %>%
    separate(contrast, into = c("Group1", "Group2"), sep = " - ") %>%
    mutate(
      Group1 = trimws(Group1),
      Group2 = trimws(Group2),
      # Filter out any lingering NAs
      p_bin = ifelse(!is.na(p.value) & p.value < 0.05, "Significant", "NS"),
      star = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ ""
      )
    ) %>%
    filter(!is.na(p_bin)) # Remove the "NA" noise
  
  # 2. Plot Matrix
  ggplot(pairs_data, aes(x = Group1, y = Group2, fill = p_bin)) +
    geom_tile(color = "white") +
    geom_text(aes(label = star), vjust = 0.7, size = 3) + 
    scale_fill_manual(values = c("Significant" = "#FF6B6B", "NS" = "grey90"), 
                      name = "Significance") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 14)
    ) +
    labs(title = title_text, 
         subtitle = "Pairwise family comparisons across all orders",
         x = "", y = "") +
    coord_fixed()
}

# Generate the fixed heatmaps
p_matrix_mass_fixed <- create_significance_heatmap(mod_heatmap_mass, "Sample_ID", "Mass Significance: Family Comparisons")
p_matrix_rel_fixed <- create_significance_heatmap(mod_heatmap_rel, "Sample_ID", "Lipid Significance: Family Comparisons")

print(p_matrix_mass_fixed)
print(p_matrix_rel_fixed)



# =========================================================
# Does size predict lipid content?
# =========================================================
mod_allometry <- lm(Log_Lipid_Total ~ Log_Mass, data = df_clean)

print("=== MODEL SUMMARY: DOES SIZE PREDICT FAT? ===")
summary(mod_allometry)

# 2. VISUALIZE THE RELATIONSHIP
ggplot(df_clean, aes(x = Log_Mass, y = Log_Lipid_Total)) +
  
  # A. The Points (Colored by Order to see if groups differ)
  geom_point(aes(color = Order), alpha = 0.6, size = 2) +
  
  # B. The Regression Line (The "Slope")
  geom_smooth(method = "lm", color = "black", linetype = "solid") +
  
  # C. Identity Line (Slope = 1) for comparison
  # If the black dashed line is steeper than this gray line, big bugs are "fatter."
  geom_abline(slope = 1, intercept = coef(mod_allometry)[1], color = "black", linetype = "dashed") +
  
  theme_bw() +
  labs(
    title = "Allometry of Fat Storage",
    subtitle = "Black Line = Actual Trend. Dotted = 1:1 Scaling (Isometric).",
    x = "Insect Size (Log Mass)",
    y = "Total Lipid Content (Log Mass)",
    color = "Order"
  )

# =========================================================
# Biplot of size vs quality
# =========================================================

# 1. PREPARE THE DATA 
df_tradeoff <- df_clean %>%
  # Add Life_History to the grouping so it stays in the summary
  group_by(Order, Sample_ID) %>%
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
  geom_text_repel(aes(label = Sample_ID), 
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
    y = "Quality: Relative Lipid Content (Proportion)"
  ) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# Print
print(plot_tradeoff_circles)


# =========================================================
# ####Adding predator brackets to the lipid vs weight
# =========================================================
# --- 1. SETUP ---
# Helper to make axes readable (Log -> Grams)
inv_log_mass <- function(x) { format(exp(x), digits = 2, scientific = FALSE) }

# Bracket Data (Using the "Safe Pair" Staggering)
diet_brackets <- data.frame(
  Species = c("E. Small-footed", "Tri-colored", "Little Brown", 
              "N. Long-eared", "Eastern Red", "Big Brown", 
              "Hoary Bat", "Nighthawk", "Whip-poor-will"),
  # X-axis (Natural Log of prey weights)
  xmin = log(c(0.002, 0.003, 0.005, 0.005, 0.01, 0.02, 0.05, 0.05, 0.1)), 
  xmax = log(c(0.015, 0.02, 0.04, 0.06, 0.10, 0.12, 0.25, 0.5, 1.0)),
  
  # Staggering Levels (1 = Lowest, 5 = Highest)
  y_level = c(1, 2, 3, 4, 5, 1, 2, 3, 4) 
)

# Position Calculations (Adapted for Log_Lipid_Total)
# We find the max Y (Total Lipid) to float the brackets above the data
max_y_lipid <- max(df_clean$Log_Lipid_Total, na.rm = TRUE)

# Step size needs to be in Log Units now. 
# 0.5 log units is roughly a 1.6x increase in mass, which is a good gap.
y_start <- max_y_lipid + 0.2
step <- 0.4
tick <- 0.1

diet_brackets$y_bar <- y_start + (diet_brackets$y_level * step)

# --- 2. THE PLOT ---
ggplot(df_clean, aes(x = Log_Mass, y = Log_Lipid_Total)) +
  
  # A. Identity Line (Isometric Scaling)
  # Note: Ensure 'mod_allometry' is defined in your environment!
  geom_abline(slope = 1, intercept = coef(mod_allometry)[1], 
              color = "grey50", linetype = "dotted", size = 0.8) +
  
  # B. The Data Points
  geom_point(aes(fill = Order), shape = 21, color = "black", size = 2.5, stroke = 0.6, alpha = 0.8) +
  
  # C. The Regression Line (Actual Trend)
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  
  # D. The Brackets
  geom_segment(data = diet_brackets, aes(x = xmin, xend = xmax, y = y_bar, yend = y_bar, color = Species), linewidth = 0.7, inherit.aes = FALSE) +
  geom_segment(data = diet_brackets, aes(x = xmin, xend = xmin, y = y_bar - tick, yend = y_bar + tick, color = Species), linewidth = 0.7, inherit.aes = FALSE) +
  geom_segment(data = diet_brackets, aes(x = xmax, xend = xmax, y = y_bar - tick, yend = y_bar + tick, color = Species), linewidth = 0.7, inherit.aes = FALSE) +
  geom_text(data = diet_brackets, aes(x = (xmin + xmax)/2, y = y_bar + 0.15, label = Species, color = Species), fontface = "bold", size = 3, show.legend = FALSE, inherit.aes = FALSE) +
  
  # --- SCALES ---
  # X-Axis: Insect Mass (Grams)
  scale_x_continuous(labels = inv_log_mass, breaks = log(c(0.002, 0.01, 0.05, 0.2, 1.0))) + 
  
  # Y-Axis: Total Lipid Mass (Grams)
  # We use inv_log_mass here too so the axis shows "0.001", "0.01" etc. instead of -6, -4
  scale_y_continuous(
    labels = inv_log_mass,
    breaks = log(c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)), # Customize these based on your data range
    limits = c(NA, max(diet_brackets$y_bar) + 0.2)
  ) +
  
  # Colors
  scale_color_brewer(palette = "Paired", name = "Predator") + 
  scale_fill_viridis_d(option = "turbo", name = "Insect Order") + 
  
  theme_bw() +
  theme(
    legend.position = "right",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  labs(
    title = "Allometry of Fat Storage",
    subtitle = "Are bigger bugs fatter? (Slope > 1 = Hyperallometry)",
    x = "Insect Mass (grams)",
    y = "Total Lipid Content (grams)"
  )

# =========================================================
# ####Relative lipid by weight
# =========================================================
# --- 1. SETUP ---
# Helper for X-axis labels (converts log mass back to grams)
inv_log_mass <- function(x) { format(exp(x), digits = 3, scientific = FALSE) }

# --- 2. THE PLOT ---
ggplot(df_clean, aes(x = Log_Mass, y = Lipid_Rel_Pct)) + 
  
  # A. The Data Points
  # Points are colored by Order to show taxonomic nutritional clustering
  geom_point(aes(fill = Order), shape = 21, color = "black", size = 3, stroke = 0.6, alpha = 0.8) +
  
  # B. The Trend Line
  # Shows the overall relationship between size and fat content
  geom_smooth(method = "lm", color = "black", linetype = "solid", linewidth = 1, se = TRUE) +
  
  # --- SCALES ---
  # Y-Axis: Percentage of lipid content
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1), 
    breaks = seq(0, 1, by = 0.1),
    expand = expansion(mult = c(0.05, 0.1))
  ) + 
  
  # X-Axis: Readable grams instead of log values
  scale_x_continuous(
    labels = inv_log_mass, 
    breaks = log(c(0.002, 0.01, 0.05, 0.25, 1.0))
  ) + 
  
  # Colors: Using the turbo palette for distinct Order identification
  scale_fill_viridis_d(option = "turbo", name = "Insect Order") + 
  
  # --- THEME & LABELS ---
  theme_classic() + 
  theme(
    legend.position = "right",
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10, color = "black"),
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_line(color = "grey95")
  ) +
  labs(
    title = "Insect Nutritional Quality vs. Body Mass",
    subtitle = "Analysis of prey availability on Manitoulin Island",
    x = "Insect Mass (grams)",
    y = "Relative Lipid Content (%)"
  )


# =========================================================
# ####Relative lipid by weight with prey backets
# =========================================================
# --- 1. SETUP ---
# Helper for X-axis
inv_log_mass <- function(x) { format(exp(x), digits = 2, scientific = FALSE) }

# Bracket Data (Optimized to 5 'Safe' Levels)
diet_brackets <- data.frame(
  Species = c("E. Small-footed", "Tri-colored", "Little Brown", 
              "N. Long-eared", "Eastern Red", "Big Brown", 
              "Hoary Bat", "Nighthawk", "Whip-poor-will"),
  # X-axis (Natural Log of prey weights)
  xmin = log(c(0.002, 0.003, 0.005, 0.005, 0.01, 0.02, 0.05, 0.05, 0.1)), 
  xmax = log(c(0.015, 0.02, 0.04, 0.06, 0.10, 0.12, 0.25, 0.5, 1.0)),
  
  # SMART STAGGERING (1 = Lowest, 5 = Highest)
  # We pair small species (left) with large species (right) to share rows
  y_level = c(1,  # Small-footed (Ends -4.2) -> Fits with Big Brown
              2,  # Tri-colored (Ends -3.9) -> Fits with Hoary
              3,  # Little Brown (Ends -3.2) -> Fits with Nighthawk
              4,  # N. Long-eared (Ends -2.8) -> Fits with Whip-poor-will
              5,  # Eastern Red (The "Middle Child" - needs its own row)
              1,  # Big Brown (Starts -3.9)
              2,  # Hoary Bat (Starts -3.0)
              3,  # Nighthawk (Starts -3.0)
              4)  # Whip-poor-will (Starts -2.3)
)

# Position Calculations
# Find the top of your actual data
max_lipid <- max(df_clean$Lipid_Rel_Pct, na.rm = TRUE)

# Start brackets just 2% above the highest dot
y_start <- max_lipid + 0.02 
step <- 0.05   # 5% gap between rows for clear readability
tick <- 0.015  # 1.5% tall ticks

diet_brackets$y_bar <- y_start + (diet_brackets$y_level * step)

# --- 2. THE PLOT ---
predator_plot = ggplot(df_clean, aes(x = Log_Mass, y = Lipid_Rel_Pct)) + 
  
  # A. The Data Points
  geom_point(aes(fill = Order), shape = 21, color = "black", size = 2.5, stroke = 0.6, alpha = 0.8) +
  geom_smooth(method = "lm", color = "grey40", linetype = "dashed", se = FALSE) +
  
  # B. The Brackets
  # Horizontal Bar
  geom_segment(data = diet_brackets, 
               aes(x = xmin, xend = xmax, y = y_bar, yend = y_bar, color = Species),
               linewidth = 0.7, inherit.aes = FALSE) +
  # Left Tick
  geom_segment(data = diet_brackets, 
               aes(x = xmin, xend = xmin, y = y_bar - tick, yend = y_bar + tick, color = Species),
               linewidth = 0.7, inherit.aes = FALSE) +
  # Right Tick
  geom_segment(data = diet_brackets, 
               aes(x = xmax, xend = xmax, y = y_bar - tick, yend = y_bar + tick, color = Species),
               linewidth = 0.7, inherit.aes = FALSE) +
  # Labels
  geom_text(data = diet_brackets,
            aes(x = (xmin + xmax)/2, y = y_bar + 0.02, label = Species, color = Species),
            fontface = "bold", size = 3, show.legend = FALSE, inherit.aes = FALSE) +
  
  # --- SCALES ---
  # Y-Axis: Shows 0%, 10%, 20%... 
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1), 
    breaks = seq(0, 1, by = 0.1), 
    limits = c(0, max(diet_brackets$y_bar) + 0.03) # Auto-adjusts to fit top bracket
  ) + 
  
  scale_x_continuous(labels = inv_log_mass, breaks = log(c(0.002, 0.01, 0.05, 0.2, 1.0))) + 
  
  # Colors
  scale_color_brewer(palette = "Paired", name = "Predator") + 
  scale_fill_viridis_d(option = "turbo", name = "Insect Order") + 
  
  theme_classic() + 
  theme(
    legend.position = "right",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  labs(
    title = "Insect Nutritional Quality vs. Predator Dietary Scope",
    subtitle = "Staggered brackets show prey size range for each species",
    x = "Insect Mass (grams)",
    y = "Relative Lipid Content (%)"
  )

# Capture the Relative Lipid Scatter Plot (Clean Version)
plot_clean_nutrition <- ggplot(df_clean, aes(x = Log_Mass, y = Lipid_Rel_Pct)) + 
  geom_point(aes(fill = Order), shape = 21, color = "black", size = 3, stroke = 0.6, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", linetype = "solid", linewidth = 1, se = TRUE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 1, by = 0.1)) + 
  scale_x_continuous(labels = inv_log_mass, breaks = log(c(0.002, 0.01, 0.05, 0.25, 1.0))) + 
  scale_fill_viridis_d(option = "turbo", name = "Insect Order") + 
  theme_classic() + 
  labs(title = "Insect Nutritional Quality vs. Body Mass", x = "Insect Mass (grams)", y = "Relative Lipid Content (%)")

# Capture the Allometry Scatter Plot (Clean Version)
plot_clean_allometry <- ggplot(df_clean, aes(x = Log_Mass, y = Log_Lipid_Total)) +
  geom_abline(slope = 1, intercept = coef(mod_allometry)[1], color = "grey60", linetype = "dotted", linewidth = 1) +
  geom_point(aes(fill = Order), shape = 21, color = "black", size = 3, stroke = 0.6, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", linetype = "solid", linewidth = 1.2, se = TRUE) +
  scale_x_continuous(labels = inv_log_mass, breaks = log(c(0.002, 0.01, 0.05, 0.25, 1.0))) + 
  scale_y_continuous(labels = inv_log_mass, breaks = log(c(0.0001, 0.001, 0.01, 0.1))) +
  scale_fill_viridis_d(option = "turbo", name = "Insect Order") + 
  theme_classic() +
  labs(title = "Allometry of Fat Storage", x = "Insect Mass (grams)", y = "Total Lipid Content (grams)")





#save plots 
# Define your specific folder name
output_folder <- "InsectNutrition"

# Create the folder if it doesn't already exist in your working directory
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Updated list including predator_plot
export_list <- list(
  "Taxon_Mass_Nested"     = plot_mass_nested,
  "Taxon_Lipid_Nested"    = plot_manual_brackets,
  "Family_Heatmap_Mass"   = p_matrix_mass_fixed,
  "Family_Heatmap_Lipid"  = p_matrix_rel_fixed,
  "Nutritional_Scatter"   = plot_clean_nutrition,
  "Fat_Allometry"         = plot_clean_allometry,
  "Predator_Diet_Scope"   = predator_plot,     # <--- Added this back in!
  "Tradeoff_Biplot"       = plot_tradeoff_circles
)

# Loop and Save to InsectNutrition
for (name in names(export_list)) {
  file_path <- file.path("InsectNutrition", paste0(name, ".png"))
  
  ggsave(
    filename = file_path,
    plot = export_list[[name]],
    width = 9, 
    height = 6.5,
    dpi = 300, # Professional conference quality
    bg = "white"
  )
  
  message("Successfully exported: ", file_path)
}
