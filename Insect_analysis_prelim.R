# Load necessary libraries
library(tidyverse)
library(emmeans)

# 1. Load Data
df <- read.csv("Insect_Nutrition_Prelim.csv", check.names = FALSE)

# 2. Data Cleaning & Metric Calculation
df_clean <- df %>%
  select(Order, Quantity, Insect_Weight, Dry_Weight, Lipid_Free_Weight) %>%
  filter(!is.na(Order), !is.na(Insect_Weight), 
         !is.na(Dry_Weight), !is.na(Lipid_Free_Weight)) %>%
  mutate(
    # Fix typos
    Order = ifelse(Order == "Lepedoptera", "Lepidoptera", Order),
    
    # Metric 1: Individual Biomass
    Biomass_Ind = Insect_Weight / Quantity,
    
    # Metric 2: Lipid Content (Wet Basis)
    Lipid_Mass = Dry_Weight - Lipid_Free_Weight,
    Lipid_Wet_Pct = Lipid_Mass / Insect_Weight
  ) %>%
  # Filter out invalid calculations
  filter(Biomass_Ind > 0, Lipid_Wet_Pct >= 0, Lipid_Wet_Pct <= 1)

# 3. Run Statistical Models
# We weight by Quantity to prioritize samples with more individuals
model_bio <- lm(Biomass_Ind ~ Order, data = df_clean, weights = Quantity)
model_lipid <- lm(Lipid_Wet_Pct ~ Order, data = df_clean, weights = Quantity)

# ---------------------------------------------------------
# 4. PRINT ANOVA RESULTS (Statistical Significance)
# ---------------------------------------------------------
print("--- ANOVA: Individual Biomass (Size) ---")
print(anova(model_bio))

print("--- ANOVA: Lipid Content (Wet %) ---")
print(anova(model_lipid))

# ---------------------------------------------------------
# 5. VISUALIZATION
# ---------------------------------------------------------

# Extract Estimated Means for Plotting
# Biomass
emm_bio <- emmeans(model_bio, ~ Order)
df_bio <- as.data.frame(emm_bio) %>% 
  mutate(Metric = "Individual Mass (g)") %>%
  rename(Value = emmean)

# Lipid
emm_lipid <- emmeans(model_lipid, ~ Order)
df_lipid <- as.data.frame(emm_lipid) %>% 
  mutate(Metric = "Lipid Content (Wet %)") %>%
  rename(Value = emmean)

# Combine datasets
all_order_data <- bind_rows(df_bio, df_lipid)

# Create the Plot
ggplot(all_order_data, aes(x = Order, y = Value, color = Order)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, linewidth = 1) +
  
  # LAYOUT:
  # . ~ Metric : Puts Metrics (Mass, Lipid) in Columns (Labels on Top)
  # scales = "free_y" : Allows the Value axis (which starts as Y) to vary between metrics
  facet_grid(. ~ Metric, scales = "free_y") + 
  
  # FLIP:
  # Puts Order on the Left (Y-axis) and Value on the Bottom (X-axis)
  coord_flip() +
  
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.position = "none") +
  
  labs(title = "Comparison of Insect Orders",
       subtitle = "Estimated Means with 95% Confidence Intervals",
       y = "Estimated Mean Value",
       x = "")

###Family level consideration
# Load necessary libraries
library(tidyverse)
library(emmeans)

# 1. Load Data
df <- read.csv("Insect_Nutrition_Prelim.csv", check.names = FALSE)

# 2. Data Cleaning & Metric Calculation
df_clean <- df %>%
  select(Order, Family, Month, Quantity, Insect_Weight, Dry_Weight, Lipid_Free_Weight) %>%
  filter(!is.na(Family), !is.na(Month), !is.na(Insect_Weight), 
         !is.na(Dry_Weight), !is.na(Lipid_Free_Weight)) %>%
  mutate(
    # Fix typos
    Order = ifelse(Order == "Lepedoptera", "Lepidoptera", Order),
    
    # METRIC 1: Individual Biomass (Wet Weight / Count)
    Biomass_Ind = Insect_Weight / Quantity,
    
    # METRIC 2: Lipid Content (Wet Basis)
    # (Dry - LipidFree) / Wet
    Lipid_Mass = Dry_Weight - Lipid_Free_Weight,
    Lipid_Wet_Pct = Lipid_Mass / Insect_Weight
  ) %>%
  # Filter out invalid calculations
  filter(Biomass_Ind > 0, Lipid_Wet_Pct >= 0, Lipid_Wet_Pct <= 1)

# 3. ANALYSIS 1: Biomass per Individual by Family
# We control for Month to account for seasonal growth.
model_bio_fam <- lm(Biomass_Ind ~ Family + Month, 
                    data = df_clean, 
                    weights = Quantity)

cat("\n--- ANOVA: Biomass per Individual ---\n")
anova(model_bio_fam)

# Get Family Means (averaged over Month)
emm_bio_fam <- emmeans(model_bio_fam, ~ Family)

# 4. ANALYSIS 2: Lipid Content (Wet) by Family
model_lipid_fam <- lm(Lipid_Wet_Pct ~ Family + Month, 
                      data = df_clean, 
                      weights = Quantity)

cat("\n--- ANOVA: Lipid Content (Wet) ---\n")
anova(model_lipid_fam)

# Get Family Means
emm_lipid_fam <- emmeans(model_lipid_fam, ~ Family)

# 5. Visualization
# We'll plot the estimated means. Faceting by Order helps organize the many families.

# Convert emmeans to data frames for plotting
bio_plot_data <- as.data.frame(emm_bio_fam)
bio_plot_data$Metric <- "Biomass (g)"

lipid_plot_data <- as.data.frame(emm_lipid_fam)
lipid_plot_data$Metric <- "Lipid % (Wet)"

# Combine
all_plot_data <- bind_rows(bio_plot_data, lipid_plot_data) %>%
  # Join with Order info so we can facet
  left_join(unique(df_clean[, c("Family", "Order")]), by = "Family")

# Plot
ggplot(all_plot_data, aes(x = Family, y = emmean, color = Order)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  facet_grid(Metric ~ Order, scales = "free", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  labs(title = "Comparison of Families (Controlled for Month)",
       y = "Estimated Mean Value",
       x = "Family")

#####Adding size for relavent families

# Load necessary libraries
library(tidyverse)
library(emmeans)

# 1. Load Data
df <- read.csv("Insect_Nutrition_Prelim.csv", check.names = FALSE)

# 2. Data Cleaning & Group Creation
df_clean <- df %>%
  select(Order, Family, Month, Quantity, Size, Insect_Weight, Dry_Weight, Lipid_Free_Weight) %>%
  filter(!is.na(Family), !is.na(Month), !is.na(Insect_Weight), 
         !is.na(Dry_Weight), !is.na(Lipid_Free_Weight)) %>%
  mutate(
    # Fix typos
    Order = ifelse(Order == "Lepedoptera", "Lepidoptera", Order),
    
    # Create the label: "Family" or "Family (Size)"
    Size_Label = ifelse(Size == "" | is.na(Size), Family, paste0(Family, " (", Size, ")")),
    Taxon_Group = Size_Label,
    
    # Metric 1: Individual Biomass
    Biomass_Ind = Insect_Weight / Quantity,
    
    # Metric 2: Lipid Content (Wet Basis)
    Lipid_Mass = Dry_Weight - Lipid_Free_Weight,
    Lipid_Wet_Pct = Lipid_Mass / Insect_Weight
  ) %>%
  filter(Biomass_Ind > 0, Lipid_Wet_Pct >= 0, Lipid_Wet_Pct <= 1)

# 3. Run Statistical Models
# We include Month to control for season
model_bio <- lm(Biomass_Ind ~ Taxon_Group + Month, data = df_clean, weights = Quantity)
model_lipid <- lm(Lipid_Wet_Pct ~ Taxon_Group + Month, data = df_clean, weights = Quantity)

# ---------------------------------------------------------
# 4. PRINT ANOVA RESULTS
# ---------------------------------------------------------
print("--- ANOVA: Individual Biomass (Size) ---")
print(anova(model_bio))

print("--- ANOVA: Lipid Content (Wet %) ---")
print(anova(model_lipid))

# ---------------------------------------------------------
# 5. VISUALIZATION
# ---------------------------------------------------------

# Extract Estimated Means
emm_bio <- emmeans(model_bio, ~ Taxon_Group)
df_bio <- as.data.frame(emm_bio) %>% 
  mutate(Metric = "Individual Mass (g)") %>%
  rename(Value = emmean)

emm_lipid <- emmeans(model_lipid, ~ Taxon_Group)
df_lipid <- as.data.frame(emm_lipid) %>% 
  mutate(Metric = "Lipid Content (Wet %)") %>%
  rename(Value = emmean)

# Combine and Add Order Info for Grouping
all_data <- bind_rows(df_bio, df_lipid) %>%
  left_join(unique(df_clean[, c("Taxon_Group", "Order")]), by = "Taxon_Group")

# Create the Plot
ggplot(all_data, aes(x = Taxon_Group, y = Value, color = Order)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3) +
  
  # LAYOUT:
  # Order ~ Metric : Rows = Order (to group families), Columns = Metric
  # scales = "free" : Allows Mass and Lipid axes to be independent
  # space = "free_y" : Makes the panels bigger/smaller based on how many families are in that Order
  facet_grid(Order ~ Metric, scales = "free", space = "free_y") +
  
  coord_flip() +
  
  theme_bw() +
  labs(title = "Nutritional Profile by Family and Size",
       subtitle = "Controlled for Month and Quantity",
       x = "", y = "Estimated Mean Value") +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "none")

