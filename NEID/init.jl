## Setup
#using Pkg
#Pkg.activate("NEID")
#Pkg.instantiate()

## Importing data with Eric's code
import StellarSpectraObservationFitting as SSOF
using EchelleInstruments, EchelleInstruments.NEID
using CSV, DataFrames, Query
using FITSIO
import SSOFApplication as SSOFA

# parsing inputs
if length(ARGS) != 0; ENV["GKSwstype"] = "100" end
csv_with_filenames = ARGS[1]
path_to_save_to = ARGS[2]
max_spectra_to_use = parse(Int, ARGS[3])

# creating manifest of files
df_filenames = CSV.read(csv_with_filenames, DataFrame)
# df_filenames = make_manifest(df_filenames)
if !(size(df_filenames,1)>=1)
   println("# Exiting early because no manifest files found.")
   exit(0)
end
df_files = SSOFA.make_manifest_neid(df_filenames[:, 1])
df_files_use = SSOFA.filter_manifest_neid(df_files, max_spectra_to_use=max_spectra_to_use)
if !(size(df_files_use,1)>=4)
   println("# Exiting early because only " * string(size(df_files_use,1)) * " files found for use by SSOF.")
   exit(0)
end

# defines NEIDLSF.NEID_lsf()
include("lsf.jl") 

mkpath(path_to_save_to)

# reformating spectra into SSOF Data objects
SSOFA.reformat_spectra(
	df_files_use,
	path_to_save_to,
	NEID,
	min_order(NEID2D()):118;
	lsf_f = NEIDLSF.neid_lsf,
	interactive=length(ARGS)==0,
	min_snr=5,
	λ_thres=4000)

# saving some of the NEID metadata
SSOFA.neid_extras(df_files_use, path_to_save_to)
