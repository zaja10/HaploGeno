devtools::load_all(".")

message("Testing load_demo_data()...")
haplo <- load_demo_data()

message("Checking object...")
print(haplo)

message("Checking visualization...")
tryCatch(
    {
        haplo$plot_manhattan()
        message("Manhattan plot: OK")
    },
    error = function(e) stop(e)
)

message("Success!")
