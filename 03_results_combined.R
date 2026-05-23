require(pmsampsize)

# 1. Define the vectors for your scenarios
target_auc <- c(0.7, 0.75, 0.8, 0.85, 0.9)
target_prev <- c(0.05, 0.1, 0.2)
n_params <- 10  # Keeping parameters constant as per your example

# 2. Create a grid of all combinations
pm_table <- expand.grid(Target_AUC = target_auc, 
                        Target_Prev = target_prev)

# 3. Initialize a column to store the results
pm_table$Sample_Size <- NA

# 4. Loop through each row to calculate sample size
for(i in 1:nrow(pm_table)) {
  
  # Calculate pmsampsize for the current row's parameters
  result <- pmsampsize(type = "b", 
                       cstatistic = pm_table$Target_AUC[i], 
                       parameters = n_params, 
                       prevalence = pm_table$Target_Prev[i])
  
  # Extract and store the sample size
  pm_table$Sample_Size[i] <- result$sample_size
  
  
}

# 5. View the final table sorted for easier reading
pm_table <- pm_table[order(pm_table$Target_Prev, pm_table$Target_AUC),]
print(pm_table)
save(pm_table, file = "pm_table.RData")


##############
load("table_2.RData")
load("pm_table.RData")

library(dplyr)


final_data <- table_2 %>%
  # 1. Rename column in the main table
  rename(n_new = required_n) %>%
  # 2. Join with the second table (renaming its column on the fly)
  left_join(
    pm_table %>% rename(n_riley = Sample_Size),
    by = c("Target_AUC", "Target_Prev")
  ) %>%
  # 3. Move n_riley to be immediately after n_new
  relocate(n_riley, .after = n_new)


table_3 <- final_data
save(table_3, file = "table_3.RData")



