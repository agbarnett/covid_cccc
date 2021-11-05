# 1_obese_impute.R
# impute missing bmi in patient data
# called by 1_obese.Rmd
# May 2021

# variable lists
imputerVars <- c("country","ethnicity","height","weight","age","sex",'apache_ii','sofa') # variables used for imputation of BMI (some also have missing)
vars_to_keep = c('pin') # other variables not used for imputation, but that need to be kept
imputedOnlyVars <- c("bmi") # could add apache, but too much missing?

# flag what variables to use in imputation (exclude exit/to from imputation, but keep them in the data); see https://rpubs.com/kaz_yos/mice-exclude
to_impute = select(patients, all_of(vars_to_keep), all_of(imputerVars), all_of(imputedOnlyVars))
allVars <- names(to_impute)
predictorMatrix <- matrix(0, ncol = length(allVars), nrow = length(allVars))
rownames(predictorMatrix) <- allVars
colnames(predictorMatrix) <- allVars
## Keep variables that actually exist in dataset
imputerVars <- intersect(unique(imputerVars), allVars)
imputerMatrix <- predictorMatrix
imputerMatrix[,imputerVars] <- 1
#
missVars <- names(to_impute)[colSums(is.na(to_impute)) > 0]
## Imputers that have missingness must be imputed.
imputedVars <- intersect(unique(c(imputedOnlyVars, imputerVars)), missVars)
imputedMatrix <- predictorMatrix
imputedMatrix[imputedVars,] <- 1
#
predictorMatrix <- imputerMatrix * imputedMatrix
## Diagonals must be zeros (a variable cannot impute itself)
diag(predictorMatrix) <- 0
#predictorMatrix

# imputation using MICE
N.impute = 5
imputed_data = mice(to_impute, m=N.impute, seed=4003, printFlag=FALSE, predictorMatrix = predictorMatrix)
