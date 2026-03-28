# =========================================================================
# Script 1: Bulk RNA-seq Machine Learning Feature Selection
# Cohort: GSE12021 (Discovery)
# =========================================================================

library(glmnet)
library(randomForest)
library(Boruta)
library(VennDiagram)

# 1. LASSO Regression (10-fold Cross-Validation)
set.seed(123)
cv_fit <- cv.glmnet(x = as.matrix(expr_train), y = group_train, 
                    family = "binomial", alpha = 1, nfolds = 10)
lasso_coef <- coef(cv_fit, s = "lambda.min")
lasso_genes <- rownames(lasso_coef)[which(lasso_coef != 0)][-1] # 排除截距

# 2. Random Forest (ntree = 500)
set.seed(123)
rf_model <- randomForest(x = expr_train, y = as.factor(group_train), 
                         ntree = 500, importance = TRUE)
rf_importance <- importance(rf_model)
rf_genes <- rownames(rf_importance)[order(rf_importance[, "MeanDecreaseGini"], decreasing = TRUE)[1:13]]

# 3. Boruta Algorithm (MaxRuns = 500, p < 0.01)
set.seed(123)
boruta_model <- Boruta(x = expr_train, y = as.factor(group_train), 
                       pValue = 0.01, maxRuns = 500)
boruta_genes <- getSelectedAttributes(boruta_model, withTentative = FALSE)

# 4. Final Core Panel Intersection
core_targets <- intersect(intersect(lasso_genes, rf_genes), boruta_genes)
print(paste("Core diagnostic biomarkers identified:", paste(core_targets, collapse = ", ")))

# Result: CD83, TNFAIP3, NFKBIA, MAFF (MAFF subsequently filtered out in validation)