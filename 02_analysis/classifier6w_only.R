set.seed(123)
dir.create("SixWeek", showWarnings = FALSE)
dir.create("SixWeekResults", showWarnings = FALSE)

# Data Preparation — verbatim from SixWeekModelDetermination.R
# (readRDS replaced with pre-loaded combinedcongo)
meco_dataset <- phyloseq2meco(combinedcongo)
meco_dataset$tidy_dataset()

meco_dataset$tax_table <- base::subset(meco_dataset$tax_table, Kingdom %in% c("k__Archaea", "k__Bacteria"))
meco_dataset$filter_pollution(taxa = c("mitochondria", "chloroplast"))
st <- subset(meco_dataset$sample_table, Type == "sample")
st$BeforeMalaria <- as.character(st$BeforeMalaria)
st$BeforeMalaria[is.na(st$BeforeMalaria)] <- "Never"
st <- subset(st, BeforeMalaria %in% c("Yes", "Never"))
st <- subset(st, Time == "6W")
st$BeforeMalaria <- factor(st$BeforeMalaria, levels = c("Yes", "Never"))
meco_dataset$sample_table <- st
meco_dataset$tidy_dataset()
## --------------------------------------------------------------------
## 2. Determine Available Taxonomic Levels
## --------------------------------------------------------------------
get_available_taxonomic_levels <- function(meco_obj) {
  # Get all taxonomic columns from tax_table
  tax_cols <- colnames(meco_obj$tax_table)
  
  # Remove 'Taxa' column if it exists (this is typically the row identifier)
  tax_cols <- tax_cols[tax_cols != "Taxa"]
  
  # Common taxonomic level mappings
  level_mappings <- c(
    "Kingdom" = "Kingdom",
    "Phylum" = "Phylum", 
    "Class" = "Class",
    "Order" = "Order",
    "Family" = "Family",
    "Genus" = "Genus",
    "Species" = "Species"
  )
  
  # Find which levels are available
  available_levels <- c()
  for(level in names(level_mappings)) {
    if(level %in% tax_cols) {
      available_levels <- c(available_levels, level)
    }
  }
  
  # If no standard names found, use all columns except Taxa
  if(length(available_levels) == 0) {
    available_levels <- tax_cols
    cat("Warning: Using all tax_table columns as taxonomic levels\n")
  }
  
  return(available_levels)
}

# Get available taxonomic levels
taxonomic_levels <- get_available_taxonomic_levels(meco_dataset)
cat("Available taxonomic levels:\n")
cat(paste(taxonomic_levels, collapse = ", "), "\n\n")

## --------------------------------------------------------------------
## 3. Robust Performance Extraction Function
## --------------------------------------------------------------------
extract_microeco_performance <- function(classifier_obj) {
  results <- list(
    accuracy = NA,
    sensitivity = NA,
    specificity = NA,
    precision = NA,
    f1_score = NA,
    balanced_accuracy = NA,
    confusion_matrix = NULL,
    success = FALSE
  )
  
  tryCatch({
    # Check multiple possible locations for confusion matrix
    cm <- NULL
    
    # Method 1: res_predict$confusionMatrix
    if(!is.null(classifier_obj$res_predict) && is.list(classifier_obj$res_predict)) {
      if("confusionMatrix" %in% names(classifier_obj$res_predict)) {
        cm <- classifier_obj$res_predict$confusionMatrix
      }
    }
    
    # Method 2: Direct access to confusionMatrix
    if(is.null(cm) && !is.null(classifier_obj$res_confusion_fit)) {
      cm <- classifier_obj$res_confusion_fit
    }
    
    # Method 3: Check res_confusion_stats
    if(is.null(cm) && !is.null(classifier_obj$res_confusion_stats)) {
      cm <- classifier_obj$res_confusion_stats
    }
    
    if(!is.null(cm)) {
      # Extract overall metrics
      if(!is.null(cm$overall) && is.vector(cm$overall)) {
        if("Accuracy" %in% names(cm$overall)) {
          results$accuracy <- as.numeric(cm$overall["Accuracy"])
        }
      }
      
      # Extract by-class metrics
      if(!is.null(cm$byClass) && is.vector(cm$byClass)) {
        metric_names <- names(cm$byClass)
        
        if("Sensitivity" %in% metric_names) {
          results$sensitivity <- as.numeric(cm$byClass["Sensitivity"])
        }
        if("Specificity" %in% metric_names) {
          results$specificity <- as.numeric(cm$byClass["Specificity"])
        }
        if("Pos Pred Value" %in% metric_names) {
          results$precision <- as.numeric(cm$byClass["Pos Pred Value"])
        }
        if("F1" %in% metric_names) {
          results$f1_score <- as.numeric(cm$byClass["F1"])
        }
        if("Balanced Accuracy" %in% metric_names) {
          results$balanced_accuracy <- as.numeric(cm$byClass["Balanced Accuracy"])
        }
        
        # Calculate balanced accuracy if not available
        if(is.na(results$balanced_accuracy) && !is.na(results$sensitivity) && !is.na(results$specificity)) {
          results$balanced_accuracy <- (results$sensitivity + results$specificity) / 2
        }
      }
      
      results$confusion_matrix <- cm
      results$success <- TRUE
      
      cat(sprintf("  Performance extracted: Accuracy=%.3f, Sensitivity=%.3f, Specificity=%.3f\n",
                  ifelse(is.na(results$accuracy), 0, results$accuracy),
                  ifelse(is.na(results$sensitivity), 0, results$sensitivity),
                  ifelse(is.na(results$specificity), 0, results$specificity)))
    } else {
      cat("  Warning: No confusion matrix found in classifier object\n")
    }
    
  }, error = function(e) {
    cat(sprintf("  Error extracting performance: %s\n", e$message))
  })
  
  return(results)
}

## --------------------------------------------------------------------
## 4. Enhanced Single Classifier Training Function
## --------------------------------------------------------------------
train_microeco_classifier <- function(dataset, method_name, method_display_name, 
                                      taxonomic_level, with_feature_selection = TRUE, verbose = TRUE) {
  
  if(verbose) cat(sprintf("\nTraining %s (%s) at %s level%s...\n", 
                          method_display_name, method_name, taxonomic_level,
                          ifelse(with_feature_selection, " with Feature Selection", "")))
  
  result <- list(
    method_code = method_name,
    method_name = method_display_name,
    taxonomic_level = taxonomic_level,
    with_fs = with_feature_selection,
    model = NULL,
    performance = NULL,
    success = FALSE,
    error_message = NULL
  )
  
  tryCatch({
    # Create classifier
    classifier <- trans_classifier$new(
      dataset = dataset, 
      y.response = "BeforeMalaria", 
      x.predictors = taxonomic_level  # Now using the specified taxonomic level
    )
    
    # Check if we have enough features to proceed
    if(!is.null(classifier$data_feature)) {
      initial_feature_count <- ncol(classifier$data_feature)
      if(verbose) cat(sprintf("  Available features at %s level: %d\n", taxonomic_level, initial_feature_count))
      
      # Skip if too few features for meaningful analysis
      if(initial_feature_count < 2) {
        stop(sprintf("Insufficient features (%d) at %s level for classification", initial_feature_count, taxonomic_level))
      }
    }
    
    # Split data
    classifier$cal_split(prop.train = 3/4)
    
    # Feature selection if requested
    if(with_feature_selection) {
      # Check number of features available
      initial_features <- ncol(classifier$data_feature)
      if(verbose) cat(sprintf("  Initial features available: %d\n", initial_features))
      
      # Skip feature selection if too few features
      if(initial_features <= 3) {
        if(verbose) cat(sprintf("  Skipping feature selection (only %d features available)\n", initial_features))
      } else {
        classifier$cal_feature_sel(boruta.maxRuns = 200, boruta.pValue = 0.01)
        
        # Check if any features were selected
        remaining_features <- ncol(classifier$data_train)
        if(verbose) cat(sprintf("  Features after selection: %d\n", remaining_features))
        
        # If no features selected, skip this model
        if(remaining_features == 0) {
          stop("No features selected by Boruta - model cannot be trained")
        }
      }
    }
    
    # Set training control
    classifier$set_trainControl(method = "cv", number = 5, savePredictions = TRUE, classProbs = TRUE)
    
    # Train model with method-specific parameters
    if(method_name == "rf") {
      classifier$cal_train(method = "rf")
    } else if(method_name == "svmRadial") {
      classifier$cal_train(method = "svmRadial", tuneLength = 5)
    } else if(method_name == "svmLinear") {
      classifier$cal_train(method = "svmLinear", tuneLength = 3)
    } else if(method_name == "glm") {
      classifier$cal_train(method = "glm")
    } else if(method_name == "knn") {
      classifier$cal_train(method = "knn", tuneGrid = expand.grid(k = c(3, 5, 7)))
    } else if(method_name == "nb") {
      classifier$cal_train(method = "nb")
    } else {
      classifier$cal_train(method = method_name)
    }
    
    # Make predictions
    classifier$cal_predict()
    
    # Calculate ROC (optional, with error handling)
    tryCatch({
      classifier$cal_ROC()
    }, error = function(e) {
      if(verbose) cat(sprintf("  Note: ROC calculation failed: %s\n", e$message))
    })
    
    # Store model
    result$model <- classifier
    
    # Extract performance
    performance <- extract_microeco_performance(classifier)
    result$performance <- performance
    result$success <- performance$success
    
    if(result$success) {
      if(verbose) cat(sprintf("[OK] %s (%s) training completed successfully\n", method_display_name, taxonomic_level))
    } else {
      if(verbose) cat(sprintf("[WARNING] %s (%s) training completed but performance extraction failed\n", method_display_name, taxonomic_level))
    }
    
  }, error = function(e) {
    result$error_message <- e$message
    if(verbose) cat(sprintf("[FAIL] Error training %s (%s): %s\n", method_display_name, taxonomic_level, e$message))
  })
  
  return(result)
}

## --------------------------------------------------------------------
## 5. Main Analysis Pipeline
## --------------------------------------------------------------------
cat("\n=== COMPREHENSIVE MULTI-TAXONOMIC MICROECO CLASSIFIER ANALYSIS ===\n")
cat("Dataset info:\n")
cat(sprintf("- Total samples: %d\n", nrow(meco_dataset$sample_table)))
cat(sprintf("- Response variable levels: %s\n", paste(levels(meco_dataset$sample_table$BeforeMalaria), collapse = ", ")))
cat(sprintf("- Total features available: %d\n", nrow(meco_dataset$otu_table)))
cat(sprintf("- Taxonomic levels to analyze: %s\n", paste(taxonomic_levels, collapse = ", ")))

# Define methods to test
methods_to_test <- list(
  "rf" = "Random Forest",
  "svmRadial" = "SVM (Radial Kernel)",
  "svmLinear" = "SVM (Linear Kernel)", 
  "glm" = "Logistic Regression",
  "knn" = "k-Nearest Neighbors",
  "nb" = "Naive Bayes"
)

# Storage for all results
all_results <- list()

## --------------------------------------------------------------------
## 6. Run All Classifiers for All Taxonomic Levels
## --------------------------------------------------------------------
cat("\n=== TRAINING ALL CLASSIFIERS ACROSS ALL TAXONOMIC LEVELS ===\n")

total_combinations <- length(taxonomic_levels) * length(methods_to_test) * 2  # *2 for with/without FS
current_combination <- 0

for(tax_level in taxonomic_levels) {
  cat(sprintf("\n--- ANALYZING TAXONOMIC LEVEL: %s ---\n", tax_level))
  
  # Quick check of feature availability at this level
  tryCatch({
    temp_classifier <- trans_classifier$new(
      dataset = meco_dataset, 
      y.response = "BeforeMalaria", 
      x.predictors = tax_level
    )
    feature_count <- ncol(temp_classifier$data_feature)
    cat(sprintf("Features available at %s level: %d\n", tax_level, feature_count))
    
    if(feature_count < 2) {
      cat(sprintf("[WARNING] Skipping %s level - insufficient features for classification\n", tax_level))
      # Still increment counters to maintain progress tracking
      current_combination <- current_combination + (length(methods_to_test) * 2)
      next
    }
    
  }, error = function(e) {
    cat(sprintf("[WARNING] Error checking %s level: %s\n", tax_level, e$message))
    current_combination <- current_combination + (length(methods_to_test) * 2)
    next
  })
  
  for(method_code in names(methods_to_test)) {
    method_name <- methods_to_test[[method_code]]
    
    # With feature selection
    current_combination <- current_combination + 1
    cat(sprintf("\n[%d/%d] ", current_combination, total_combinations))
    
    result_with_fs <- train_microeco_classifier(
      dataset = meco_dataset,
      method_name = method_code,
      method_display_name = method_name,
      taxonomic_level = tax_level,
      with_feature_selection = TRUE
    )
    all_results[[paste0(tax_level, "_", method_code, "_with_fs")]] <- result_with_fs
    
    # Without feature selection
    current_combination <- current_combination + 1
    cat(sprintf("\n[%d/%d] ", current_combination, total_combinations))
    
    result_without_fs <- train_microeco_classifier(
      dataset = meco_dataset,
      method_name = method_code,
      method_display_name = method_name,
      taxonomic_level = tax_level,
      with_feature_selection = FALSE
    )
    all_results[[paste0(tax_level, "_", method_code, "_without_fs")]] <- result_without_fs
    
    # Brief pause
    Sys.sleep(1)
  }
}

## --------------------------------------------------------------------
## 7. Create Comprehensive Results Summary
## --------------------------------------------------------------------
cat("\n=== CREATING RESULTS SUMMARY ===\n")

create_results_dataframe <- function(results_list) {
  df <- data.frame(
    Taxonomic_Level = character(),
    Method = character(),
    Feature_Selection = character(),
    Success = logical(),
    Accuracy = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    Precision = numeric(),
    F1_Score = numeric(),
    Balanced_Accuracy = numeric(),
    Error_Message = character(),
    stringsAsFactors = FALSE
  )
  
  for(result_name in names(results_list)) {
    result <- results_list[[result_name]]
    
    # Extract performance metrics
    acc <- ifelse(!is.null(result$performance) && !is.na(result$performance$accuracy), 
                  result$performance$accuracy, 0)
    sens <- ifelse(!is.null(result$performance) && !is.na(result$performance$sensitivity), 
                   result$performance$sensitivity, 0)
    spec <- ifelse(!is.null(result$performance) && !is.na(result$performance$specificity), 
                   result$performance$specificity, 0)
    prec <- ifelse(!is.null(result$performance) && !is.na(result$performance$precision), 
                   result$performance$precision, 0)
    f1 <- ifelse(!is.null(result$performance) && !is.na(result$performance$f1_score), 
                 result$performance$f1_score, 0)
    bal_acc <- ifelse(!is.null(result$performance) && !is.na(result$performance$balanced_accuracy), 
                      result$performance$balanced_accuracy, 0)
    
    df <- rbind(df, data.frame(
      Taxonomic_Level = result$taxonomic_level,
      Method = result$method_name,
      Feature_Selection = ifelse(result$with_fs, "With FS", "Without FS"),
      Success = result$success,
      Accuracy = acc,
      Sensitivity = sens,
      Specificity = spec,
      Precision = prec,
      F1_Score = f1,
      Balanced_Accuracy = bal_acc,
      Error_Message = ifelse(is.null(result$error_message), "", result$error_message),
      stringsAsFactors = FALSE
    ))
  }
  
  return(df)
}

# Create summary dataframe
results_summary <- create_results_dataframe(all_results)

# Print summary by taxonomic level
cat("\n=== PERFORMANCE SUMMARY BY TAXONOMIC LEVEL ===\n")
for(tax_level in taxonomic_levels) {
  cat(sprintf("\n--- %s Level ---\n", tax_level))
  level_results <- results_summary[results_summary$Taxonomic_Level == tax_level, ]
  print(level_results[, c("Method", "Feature_Selection", "Success", "Accuracy", "Sensitivity", "Specificity", "Balanced_Accuracy")])
}

## --------------------------------------------------------------------
## 8. Identify Best Models
## --------------------------------------------------------------------
cat("\n=== BEST MODEL IDENTIFICATION ===\n")

successful_results <- results_summary[results_summary$Success == TRUE, ]

if(nrow(successful_results) > 0) {
  # Sort by balanced accuracy
  successful_results <- successful_results[order(successful_results$Balanced_Accuracy, decreasing = TRUE), ]
  
  cat("Top 10 performing models across all taxonomic levels:\n")
  top_n <- min(10, nrow(successful_results))
  for(i in 1:top_n) {
    cat(sprintf("%d. %s - %s (%s): Balanced Accuracy = %.3f, Accuracy = %.3f\n", 
                i, 
                successful_results$Taxonomic_Level[i],
                successful_results$Method[i],
                successful_results$Feature_Selection[i],
                successful_results$Balanced_Accuracy[i],
                successful_results$Accuracy[i]))
  }
  
  # Best overall model
  best_model <- successful_results[1, ]
  cat(sprintf("\n>>> BEST OVERALL MODEL: %s - %s (%s)\n", 
              best_model$Taxonomic_Level, best_model$Method, best_model$Feature_Selection))
  cat(sprintf("   Performance: Accuracy=%.3f, Sensitivity=%.3f, Specificity=%.3f, Balanced Accuracy=%.3f\n",
              best_model$Accuracy, best_model$Sensitivity, best_model$Specificity, best_model$Balanced_Accuracy))
  
  # Best model for each taxonomic level
  cat("\n=== BEST MODEL PER TAXONOMIC LEVEL ===\n")
  for(tax_level in taxonomic_levels) {
    level_results <- successful_results[successful_results$Taxonomic_Level == tax_level, ]
    if(nrow(level_results) > 0) {
      best_for_level <- level_results[1, ]
      cat(sprintf("%s: %s (%s) - Balanced Accuracy = %.3f\n",
                  tax_level, best_for_level$Method, best_for_level$Feature_Selection, 
                  best_for_level$Balanced_Accuracy))
    } else {
      cat(sprintf("%s: No successful models\n", tax_level))
    }
  }
  
} else {
  cat("[WARNING] No successful models found. Check for errors in the analysis.\n")
}

## --------------------------------------------------------------------
## 9. Generate Enhanced Visualizations
## --------------------------------------------------------------------
cat("\n=== GENERATING VISUALIZATIONS ===\n")

if(nrow(successful_results) > 0) {
  
  # Performance comparison plot by taxonomic level
  library(reshape2)
  plot_data <- successful_results[, c("Taxonomic_Level", "Method", "Feature_Selection", "Accuracy", "Sensitivity", "Specificity", "Balanced_Accuracy")]
  melted_data <- melt(plot_data, id.vars = c("Taxonomic_Level", "Method", "Feature_Selection"))
  
  performance_plot <- ggplot(melted_data, aes(x = Method, y = value, fill = Feature_Selection)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    facet_grid(variable ~ Taxonomic_Level, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text = element_text(size = 9)) +
    labs(title = "Classifier Performance Across All Taxonomic Levels",
         subtitle = sprintf("Analysis of %d samples with %s outcome classification", 
                            nrow(meco_dataset$sample_table),
                            paste(levels(meco_dataset$sample_table$BeforeMalaria), collapse = " vs ")),
         x = "Classification Method",
         y = "Performance Metric",
         fill = "Feature Selection") +
    scale_fill_manual(values = c("With FS" = "#2E86C1", "Without FS" = "#E74C3C")) +
    ylim(0, 1)
  
  ggsave("SixWeek/SixWeekmicrobiome_multi_taxonomic_performance.pdf", performance_plot, width = 20, height = 12, dpi = 300)
  cat("[OK] Multi-taxonomic performance comparison plot saved to SixWeek folder\n")
  
  # Summary heatmap showing best performance per taxonomic level
  summary_for_heatmap <- successful_results %>%
    group_by(Taxonomic_Level, Method) %>%
    summarise(Best_Balanced_Accuracy = max(Balanced_Accuracy), .groups = 'drop')
  
  heatmap_plot <- ggplot(summary_for_heatmap, aes(x = Method, y = Taxonomic_Level, fill = Best_Balanced_Accuracy)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.3f", Best_Balanced_Accuracy)), color = "white", fontface = "bold") +
    scale_fill_gradient2(low = "#E74C3C", mid = "#F39C12", high = "#27AE60", 
                         midpoint = 0.6, name = "Balanced\nAccuracy") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Best Performance Heatmap by Taxonomic Level and Method",
         subtitle = "Values shown are best balanced accuracy (with or without feature selection)",
         x = "Classification Method",
         y = "Taxonomic Level")
  
  ggsave("SixWeek/SixWeekmicrobiome_performance_heatmap.pdf", heatmap_plot, width = 12, height = 8, dpi = 300)
  cat("[OK] Performance heatmap saved to SixWeek folder\n")
  
}

## --------------------------------------------------------------------
## 10. Save All Results
## --------------------------------------------------------------------
cat("\n=== SAVING RESULTS ===\n")

# Save summary table to SixWeek folder
write.csv(results_summary, "SixWeek/SixWeekmicrobiome_multi_taxonomic_results.csv", row.names = FALSE)
cat("[OK] Results summary saved to SixWeek folder\n")

# Save R workspace to root folder
save(all_results, results_summary, meco_dataset, successful_results, taxonomic_levels,
     file = "SixWeekmicrobiome_multi_taxonomic_analysis.RData")
cat("[OK] Complete analysis saved as RData file in root folder\n")

## --------------------------------------------------------------------
## 11. Final Summary Report
## --------------------------------------------------------------------
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("MULTI-TAXONOMIC MICROBIOME CLASSIFICATION ANALYSIS - FINAL REPORT\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

cat(sprintf("Dataset: Congo Malaria Microbiome Study\n"))
cat(sprintf("Target Variable: Before Malaria (Yes/Never)\n"))
cat(sprintf("Total Samples: %d\n", nrow(meco_dataset$sample_table)))
cat(sprintf("Taxonomic Levels Analyzed: %s\n", paste(taxonomic_levels, collapse = ", ")))
cat(sprintf("Feature Selection: Boruta Algorithm\n"))
cat(sprintf("Cross-Validation: 5-fold CV\n\n"))

cat(sprintf("Total Models Trained: %d\n", length(all_results)))
cat(sprintf("Models per Taxonomic Level: %d (%d with FS, %d without FS)\n", 
            length(methods_to_test) * 2,
            length(methods_to_test), 
            length(methods_to_test)))

cat(sprintf("Successful Models: %d\n", nrow(successful_results)))
cat(sprintf("Failed Models: %d\n", length(all_results) - nrow(successful_results)))

if(nrow(successful_results) > 0) {
  cat(sprintf("\nBest Overall Model: %s - %s (%s)\n", 
              best_model$Taxonomic_Level, best_model$Method, best_model$Feature_Selection))
  cat(sprintf("Best Model Accuracy: %.1f%%\n", best_model$Accuracy * 100))
  cat(sprintf("Best Model Balanced Accuracy: %.1f%%\n", best_model$Balanced_Accuracy * 100))
  
  # Performance interpretation
  if(best_model$Balanced_Accuracy > 0.75) {
    interpretation <- "Excellent predictive performance - strong microbiome signal for malaria history"
  } else if(best_model$Balanced_Accuracy > 0.65) {
    interpretation <- "Good predictive performance - moderate microbiome signal detected"
  } else if(best_model$Balanced_Accuracy > 0.55) {
    interpretation <- "Moderate predictive performance - weak but detectable microbiome signal"
  } else {
    interpretation <- "Limited predictive performance - microbiome signal for malaria history is weak"
  }
  
  cat(sprintf("\nInterpretation: %s\n", interpretation))
  
  # Taxonomic level analysis
  cat(sprintf("\nTaxonomic Level Performance Summary:\n"))
  for(tax_level in taxonomic_levels) {
    level_results <- successful_results[successful_results$Taxonomic_Level == tax_level, ]
    if(nrow(level_results) > 0) {
      avg_performance <- mean(level_results$Balanced_Accuracy)
      best_performance <- max(level_results$Balanced_Accuracy)
      cat(sprintf("- %s: Avg=%.3f, Best=%.3f (%d successful models)\n", 
                  tax_level, avg_performance, best_performance, nrow(level_results)))
    } else {
      cat(sprintf("- %s: No successful models\n", tax_level))
    }
  }
}

cat("\nGenerated Files:\n")
cat("=== SixWeek Folder ===\n")
cat("- SixWeekmicrobiome_multi_taxonomic_results.csv: Detailed performance metrics\n")
cat("- SixWeekmicrobiome_multi_taxonomic_performance.pdf: Performance comparison plot\n")
cat("- SixWeekmicrobiome_performance_heatmap.pdf: Summary heatmap\n")
cat("\n=== Root Folder ===\n")
cat("- SixWeekmicrobiome_multi_taxonomic_analysis.RData: Complete R workspace\n")

