# Define the path to the .Rprofile file in your project folder
rprofile_path <- file.path(getwd(), ".Rprofile")

# Define the code that needs to run every time R starts
startup_code <- c(
  "local({",
  "  # 1. Only run this on macOS/Linux systems where SDKs are an issue",
  "  if (.Platform$OS.type == 'unix') {",
  "    ",
  "    # 2. Find the current location of the Xcode SDK",
  "    sdk_path <- tryCatch(",
  "      system('xcrun --show-sdk-path', intern = TRUE, ignore.stderr = TRUE),",
  "      error = function(e) NULL",
  "    )",
  "    ",
  "    # 3. If an SDK is found, create a local Makevars file",
  "    if (!is.null(sdk_path) && length(sdk_path) > 0) {",
  "      # Find clang library path for linker",
  "      clang_lib_path <- tryCatch({",
  "          files <- list.files(\"/Library/Developer/CommandLineTools/usr/lib/clang\", ",
  "                             pattern = \"libclang_rt.osx.a\", ",
  "                             recursive = TRUE, ",
  "                             full.names = TRUE)",
  "          if (length(files) > 0) dirname(files[1]) else NULL",
  "      }, error = function(e) NULL)",
  "      ",
  "      ldflags <- paste0('-isysroot ', sdk_path)",
  "      if (!is.null(clang_lib_path)) {",
  "          ldflags <- paste0(ldflags, ' -L', clang_lib_path)",
  "      }",
  "      makevars_content <- paste0(",
  "        'CPPFLAGS += -isysroot ', sdk_path, ' -I', sdk_path, '/usr/include/c++/v1\\n',",
  "        'CXXFLAGS += -isysroot ', sdk_path, ' -I', sdk_path, '/usr/include/c++/v1\\n',",
  "        'LDFLAGS  += ', ldflags",
  "      )",
  "      ",
  "      # Save it to 'Makevars.local' in the project folder",
  "      local_makevars <- file.path(getwd(), 'Makevars.local')",
  "      writeLines(makevars_content, local_makevars)",
  "      ",
  "      # 4. Force R to use this file instead of the locked system file",
  "      Sys.setenv(R_MAKEVARS_USER = local_makevars)",
  "      message('✅ Build environment configured: Using Makevars.local')",
  "    }",
  "  }",
  "})"
)

# Write this code into the .Rprofile file
writeLines(startup_code, rprofile_path)

message("Success! .Rprofile created. Please restart R now.")
