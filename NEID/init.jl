## Setup
using Pkg
Pkg.activate("NEID")
Pkg.instantiate()

## Importing data with Eric's code
import StellarSpectraObservationFitting as SSOF
using EchelleInstruments, EchelleInstruments.NEID
using CSV, DataFrames, Query
using FITSIO
import SSOFApplication as SSOFA

# parsing inputs
csv_with_filenames = ARGS[1]
path_to_save_to = ARGS[2]

# creating manifest of files
df_filenames = CSV.read(csv_with_filenames, DataFrame)
# df_filenames = make_manifest(df_filenames)
df_files = SSOFA.make_manifest_neid(df_filenames[:, 1])
df_files_use = SSOFA.filter_manifest_neid(df_files, max_spectra_to_use=200)

# defines NEIDLSF.NEID_lsf()
include("lsf.jl") 

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