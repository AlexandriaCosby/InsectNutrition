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
qqnorm(df_clean$Log_Lipid_Rel, main="Q-Q Plot: Log Relative Lipid")
qqline(df_clean$Log_Lipid_Rel, col="red", lwd=2)

# Reset plot layout to default
par(mfrow=c(1,1))

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

##Aquatic addtion
habitat_data <- tribble(
  ~Family,           ~Habitat_Type,
  "Geometridae",     "terrestrial",
  "Crambidae",       "mixed",        # Many are terrestrial, but Acentropinae are aquatic
  "Noctuidae",       "terrestrial",
  "Lasiocampidae",   "terrestrial",
  "Nematocera",      "semi-aquatic", # Broad suborder, but overwhelmingly aquatic/semi-aquatic in larval stages (e.g., midges, mosquitoes)
  "Carabidae",       "terrestrial",  # Mostly terrestrial, though some are riparian
  "Brachycera",      "mixed",        # Broad suborder; includes both terrestrial flies and semi-aquatic (e.g., Tabanidae, Stratiomyidae)
  "Leptoceridae",    "semi-aquatic", # Caddisfly: aquatic larvae, terrestrial adults
  "Notodontidae",    "terrestrial",
  "Elateridae",      "terrestrial",
  "Erebidae",        "terrestrial",
  "Pyralidae",       "terrestrial",
  "Microlep",        "terrestrial",  # Broad grouping, overwhelmingly terrestrial
  "Gelechiidae",     "terrestrial",
  "Caenidae",        "semi-aquatic", # Mayfly: aquatic larvae, terrestrial adults
  "Hydroptilidae",   "semi-aquatic", # Caddisfly: aquatic larvae, terrestrial adults
  "Limnephilidae",   "semi-aquatic", # Caddisfly: aquatic larvae, terrestrial adults
  "Heptageniidae",   "semi-aquatic", # Mayfly: aquatic larvae, terrestrial adults
  "Ephemeridae",     "semi-aquatic", # Mayfly: aquatic larvae, terrestrial adults
  "Cosmopterigidae", "terrestrial",
  "Heteroceridae",   "semi-aquatic", # Mud-loving beetles, live in riparian zones
  "Staphylinidae",   "terrestrial",  # Mostly terrestrial, some live near water
  "Tortricidae",     "terrestrial",
  "Silphidae",       "terrestrial",
  "Phryganeidae",    "semi-aquatic", # Caddisfly: aquatic larvae, terrestrial adults
  "Helicopsychidae", "semi-aquatic"  # Caddisfly: aquatic larvae, terrestrial adults
)


df_clean <- df_clean %>%
  # 1. Remove old columns if they exist to prevent name clutter (.x, .y)
  select(-any_of(c("Habitat_Type"))) %>%
  
  # 2. Left join the new habitat_data table
  # This keeps all rows in df_clean and adds the habitat info where it matches
  left_join(habitat_data, by = "Family")


########################
#Seeing if I should include month in my models
#######################
#Mass
# Model without Month
mod_no_month_mass <- lmer(Log_Mass ~ Order + Mouth_Type + (1 | Family), data = df_clean, REML = FALSE)
# Model with Month
mod_with_month_mass <- lmer(Log_Mass ~ Order + Mouth_Type + (1 | Family) + (1 | Month), data = df_clean, REML = FALSE)
# Compare them
anova(mod_no_month_mass, mod_with_month_mass)

#Lipid
# Model without Month
mod_no_month_lipid <- lmer(Log_Lipid_Total ~ Order + Mouth_Type + (1 | Family), data = df_clean, REML = FALSE)
# Model with Month
mod_with_month_lipid <- lmer(Log_Lipid_Total ~ Order + Mouth_Type + (1 | Family) + (1 | Month), data = df_clean, REML = FALSE)
# Compare them
anova(mod_no_month_lipid, mod_with_month_lipid)

#Rel Lipid
# Model without Month
mod_no_month <- lmer(Log_Lipid_Rel ~ Order + Mouth_Type + (1 | Family), data = df_clean, REML = FALSE)
# Model with Month
mod_with_month <- lmer(Log_Lipid_Rel ~ Order + Mouth_Type + (1 | Family) + (1 | Month), data = df_clean, REML = FALSE)
# Compare them
anova(mod_no_month, mod_with_month)

## Dont bother including month

# =========================================================
# Justify size separation SAYS Seperate by size 
# =========================================================
library(dplyr)
library(stringr)

# 1. Prepare the data without aggressive filtering
df_direct <- df_clean %>%
  mutate(
    Base_Family = str_trim(str_remove(Family, "\\(>1\\)|\\(<1\\)")),
    Size_Tag = case_when(
      str_detect(Family, ">1") ~ ">1",
      str_detect(Family, "<1") ~ "<1",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Size_Tag))

# 2. Get the list of all families that have at least two entries
all_families <- unique(df_direct$Base_Family)

# 3. Force the comparison for every family
final_comparison <- data.frame()

for (fam in all_families) {
  sub <- df_direct %>% filter(Base_Family == fam)
  
  # Only attempt if there are at least 2 distinct size groups present
  if (length(unique(sub$Size_Tag)) == 2) {
    
    # Using a simple t-test; if it fails (e.g., n=1 in a group), it skips silently
    test <- try(t.test(Log_Lipid_Rel ~ Size_Tag, data = sub), silent = TRUE)
    
    if (!inherits(test, "try-error")) {
      final_comparison <- rbind(final_comparison, data.frame(
        Family = fam,
        P_Value = round(test$p.value, 4),
        Mean_Small = round(as.numeric(test$estimate[1]), 3),
        Mean_Large = round(as.numeric(test$estimate[2]), 3)
      ))
    }
  }
}

# 4. View results and make the decision for your Ontario Bat Network abstract
if (nrow(final_comparison) > 0) {
  final_comparison$Significant <- ifelse(final_comparison$P_Value < 0.05, "KEEP SEPARATE", "COMBINE")
  print(final_comparison)
} else {
  print("Check your 'Family' column strings; no matches for '<1' or '>1' were found.")
}


#############THIS FEELS VERY LEGIT 
#1. Taxonomy
# Nested model: Family is nested within Order (Order / Family)
mod_tax <- lm(Log_Lipid_Rel ~ Order / Sample_ID, data = df_clean)
Anova(mod_tax, type = "II")

# 2. Functional Model
mod_func <- lm(Log_Lipid_Rel ~ Order + Mouth_Type, data = df_clean)
Anova(mod_func, type = "II")

# 3. Environmental Model
mod_habitat <- lm(Log_Lipid_Rel ~ Order + Habitat_Type, data = df_clean)
Anova(mod_habitat, type = "II")


#####Visualize, pulling our order differences
library(emmeans)

order_means <- emmeans(mod_tax, ~ Order)

order_pairs <- contrast(order_means, method = "pairwise")

order_pairs_df <- as.data.frame(order_pairs)

significant_orders <- subset(order_pairs_df, p.value < 0.05)

print(significant_orders[, c("contrast", "estimate", "p.value")])

###set up box plot style 
library(ggplot2)
library(ggsignif)
library(emmeans)

create_boxplot_clean <- function(data, x_var, y_var, fill_var, title_text, y_label, model = NULL) {
  
  # 1. Base Plot (Vertical)
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[fill_var]])) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.3, color = "black") +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, color = "darkgray") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
          axis.title.y = element_text(face = "bold")) +
    labs(title = title_text, y = y_label, x = "")
  
  # 2. Manual Nuclear Option
  if (!is.null(model)) {
    # Get the order of categories on the plot
    levs <- levels(as.factor(data[[x_var]]))
    
    # Manually defining your 6 significant pairs as numeric coordinates
    # 1=Coleoptera, 2=Diptera, 3=Ephemeroptera, 4=Lepidoptera, 5=Trichoptera
    # These match your emmeans printout exactly
    xmin_vec <- c(1, 1, 2, 2, 3, 3) 
    xmax_vec <- c(2, 3, 4, 5, 4, 5)
    
    # Set fixed heights for the 6-rung ladder
    max_y <- max(data[[y_var]], na.rm = TRUE)
    ladder <- seq(from = max_y + 0.5, by = 0.5, length.out = 6)
    
    p <- p + geom_signif(
      xmin = xmin_vec,
      xmax = xmax_vec,
      y_position = ladder,
      annotations = rep("*", 6),
      tip_length = 0.02,
      color = "black",
      textsize = 7, # Larger stars for the poster
      vjust = 0.5
    ) +
      # Force the Y-axis tall enough so no brackets get clipped
      expand_limits(y = max(ladder) + 0.5)
  }
  
  return(p)
}
# This creates the missing object
mod_viz_order <- lm(Log_Lipid_Rel ~ Order, data = df_clean)

plot_order <- create_boxplot_clean(data = df_clean, 
                                   x_var = "Order", 
                                   y_var = "Log_Lipid_Rel", 
                                   fill_var = "Order", 
                                   title_text = "Relative Lipid Content by Taxonomic Order", 
                                   y_label = "Log Relative Lipid Content", 
                                   model = mod_viz_order)

plot_order


















library(ggplot2)

# Calculate means and standard errors per Sample_ID
df_summary <- df_clean %>%
  group_by(Order, Sample_ID) %>%
  summarize(
    mean_lipid = mean(Log_Lipid_Rel, na.rm = TRUE),
    se_lipid = sd(Log_Lipid_Rel, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(df_summary, aes(x = reorder(Sample_ID, mean_lipid), y = mean_lipid, color = Order)) +
  geom_pointrange(aes(ymin = mean_lipid - se_lipid, ymax = mean_lipid + se_lipid), size = 0.8) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, face = "bold"),
        legend.position = "bottom") +
  labs(title = "Ranked Lipid Content by Family (Sample_ID)",
       x = "Family / Size Category",
       y = "Mean Log Relative Lipid (+/- SE)")



























#############################
# Taxonomic Models
mod_mass_tax <- lm(Log_Mass ~ Order / Family, data = df_clean)
mod_lipid_tax <- lm(Log_Lipid_Total ~ Order / Family, data = df_clean)
mod_rel_tax <- lm(Log_Lipid_Rel ~ Order / Family, data = df_clean)

# Functional Models
mod_mass_func <- lm(Log_Mass ~ Order + Mouth_Type, data = df_clean)
mod_lipid_func <- lm(Log_Lipid_Total ~ Order + Mouth_Type, data = df_clean)
mod_rel_func <- lm(Log_Lipid_Rel ~ Order + Mouth_Type, data = df_clean)

#Habitat model
mod_mass_habitat <- lm(Log_Mass ~ Order + Habitat_Type, data = df_clean)
mod_lipid_habitat <- lm(Log_Lipid_Total ~ Order + Habitat_Type, data = df_clean)
mod_rel_habitat <- lm(Log_Lipid_Rel ~ Order + Habitat_Type, data = df_clean)



library(performance) # Useful for comparing models

# Function to compare the two models for a single metric
compare_models <- function(tax_mod, func_mod, hab_mod, label) {
  cat("\n--- Comparison for:", label, "---\n")
  print(compare_performance(tax_mod, func_mod, rank = TRUE))
}

# Run comparisons
# Comparison for Mass
compare_performance(mod_mass_tax, mod_mass_func, mod_mass_habitat, rank = TRUE)

# Comparison for Total Lipid
compare_performance(mod_lipid_tax, mod_lipid_func, mod_lipid_habitat, rank = TRUE)

# Comparison for Relative Lipid
compare_performance(mod_rel_tax, mod_rel_func, mod_rel_habitat, rank = TRUE)



# Boxplot comparing Mouth vs No Mouth for Relative Lipid
ggplot(df_clean, aes(x = Mouth_Type, y = Log_Lipid_Rel, fill = Mouth_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  scale_fill_manual(values = c("no" = "#D55E00", "yes" = "#0072B2")) +
  theme_classic() +
  labs(title = "Energy Density by Functional Mouth Type",
       subtitle = "Capital breeders (no mouth) show higher energy density",
       x = "Mouthparts Present",
       y = "Log Relative Lipid Content") +
  theme(legend.position = "none")

