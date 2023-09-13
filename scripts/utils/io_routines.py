from pathlib import Path
import pandas as pd
import datetime

def read_table(file, removeNA=False, sep=","):
	file_path = Path(file)
	#try:
	if not file_path.exists():
		raise IOError("Error while reading csv file. File does not exist: \n" + str(file_path))
	table = pd.read_csv(file_path, sep=sep)
	if(removeNA):
		table = table.dropna()
	return(table)

def write_table(table, file, sep=","):
	file_path = Path(file)
	if table.empty:
		raise IOError("Error while writing table " + str(file_path) + ".\nTable does not exist.")
	file_path.parent.mkdir(parents=True, exist_ok=True)  
	table.to_csv(file, index=False, sep=sep, header=True)

def write_phi_per_bin_table(phi_per_bin_table, path, suffix=""):
  file = Path().joinpath(path +  "/phi_estimates_per_bin" + suffix + ".csv")
  print("Write phi estimates per bin " + str(file))
  write_table(phi_per_bin_table, file = file)

def write_sequence_info_per_day_table(seq_info_perDay_table, path, suffix=""):
  file = Path().joinpath(path +  "/sequence_stats_per_day" + suffix + ".csv")
  print("Write sequence info table " + str(file))
  write_table(seq_info_perDay_table, file = file)

def write_smoothed_phi_table(smoothed_phi_table, path, suffix=""):
  file = Path().joinpath(path +  "/smoothed_phi_estimates" + suffix + ".csv")
  print("Write smoothed phi table " + str(file))
  write_table(smoothed_phi_table, file = file)

def write_min_incidence_table(mi_table, path, suffix=""):
  file = Path().joinpath(path +  "/minimal_incidence" + suffix + ".csv")
  print("Write minimal incidence table " + str(file))
  write_table(mi_table, file = file)

	
def read_and_extract_snv_file(snv_file):
	snv_file_path = Path(snv_file)
	if not snv_file_path.is_file():
		raise IOError("Error while reading snv file " + str(snv_file_path) + ".\nFile does not exist.")
	
	print("Read table\n")

	### read table with SNVs
	seq_info_table = pd.read_csv(snv_file_path)
	
	# turn NaNs into empty string
	seq_info_table["dna_profile"] = seq_info_table["dna_profile"].fillna("")

	# Assert the needed columns are there (as in covSonar output)
	required_colnames = {"date", "dna_profile"}
	if not required_colnames.issubset(seq_info_table.columns):
		raise ValueError("Error in snv table: \nThe comma-separated snv table requires the columns date and dna_profile!")

	# Assert that that the dates are in the correct format
	try:
		seq_info_table["date"].apply(datetime.date.fromisoformat) 
	except ValueError:
		raise ValueError("Error in date format: Please provide dates in the format yyyy-mm-dd!")

	
	### extract important columns 
	minDate = datetime.date.fromisoformat(min(seq_info_table['date']))
	print("Extract columns\n")

	seq_info_short_table = seq_info_table.rename(columns = {"dna_profile": "snvs"})[['date', 'snvs']]
	seq_info_short_table['t'] = pd.to_timedelta(seq_info_table["date"].apply(datetime.date.fromisoformat) - minDate).dt.days
	seq_info_short_table = seq_info_short_table.sort_values("t")
	seq_info_short_table = seq_info_short_table.reset_index()

	return(seq_info_short_table)
