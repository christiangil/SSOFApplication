## Setup
using Pkg
Pkg.activate("NEID")
Pkg.instantiate()

## Importing
using CSV, DataFrames
import SSOFApplication as SSOFA

# parsing inputs
base_folder = ARGS[1]  # something like "C:\\user\\test"

list_of_filenames = SSOFA.collect_filenames(base_folder)
df = DataFrame(Filename=list_of_filenames)
CSV.write(base_folder*"/merged_manifest.csv", df)

