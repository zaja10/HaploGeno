devtools::load_all()

# Create dummy object
obj <- HaploObject$new(tempfile())

# Mock data to make summary interesting
obj$geno <- matrix(1:100, 10, 10)
obj$blocks <- data.table::data.table(BlockID = 1:5, Start = c(1, 3, 5, 7, 9), End = c(2, 4, 6, 8, 10))
obj$significance <- data.table::data.table(
    BlockID = 1:5,
    Variance = runif(5),
    PVE = c(0.1, 0.05, 0.02, 0.2, 0.01),
    PVE_Adj = c(0.11, 0.055, 0.022, 0.22, 0.011) # Adjusted is higher
)
obj$model_info <- list(lambda = 5.5)

# Run summary
cat("Running summary()...\n")
summary(obj)
