tryCatch(
    {
        devtools::load_all("/Users/zaldiss/Library/Mobile Documents/com~apple~CloudDocs/Antigravity/Elrond/HaploGeno", quiet = TRUE)

        obj <- HaploObject$new("test")

        # Check new methods exist
        methods <- c("plot_manhattan", "plot_factor_heatmap", "plot_block_sizes", "plot_hrm", "plot_volcano")

        for (m in methods) {
            if (!is.function(obj[[m]])) {
                stop(paste("Method", m, "not found!"))
            } else {
                message(paste("Method", m, "exists"))
            }
        }

        message("Verification Successful!")
    },
    error = function(e) {
        message("Verification Failed: ", e$message)
        quit(status = 1)
    }
)
