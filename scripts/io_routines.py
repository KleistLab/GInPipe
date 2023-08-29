from pathlib import Path
import pandas as pd
import datetime


# TODO: io routines to handle more speaking error?
# TODO: Add help message
def read_table(file, removeNA=False):
	file_path = Path(file)
	try:
		if not file_path.exists():
			raise IOError("Error while reading csv file. File does not exist: \n" + str(file_path))
		table = pd.read_csv(file_path)
		if(removeNA):
			table = table.dropna()
	except Exception as e:
		print(e)
	else:
		return(table)

def write_table(table, file):
	file_path = Path(file)
	try:
		if table.empty():
			raise IOError("Error while writing table " + str(file_path) + ".\nTable does not exist.")
		file_path.parent.mkdir(parents=True, exist_ok=True)  
		table.to_csv(file, index=False, sep=",", header=True)
	except Exception as e:
		print(e)
	
def read_and_extract_snv_file(snv_file):
	snv_file_path = Path(snv_file)
	try:
		if not snv_file_path.is_file():
			raise IOError("Error while reading snv file " + str(snv_file_path) + ".\nFile does not exist.")
		
		print("Read table\n")
		#tic("Read table")
		### read table with SNVs
		seq_info_table = pd.read_csv(snv_file_path)
		
		# Assert the needed columns are there (as in covSonar output)
		required_colnames = {"date", "dna_profile"}
		if not required_colnames.issubset(seq_info_table.columns):
			raise ValueError("Error in snv table: \nThe comma-separated snv table requires the columns date and dna_profile!")
		#if(sum(!required_colnames %in% colnames(seq_info.table)) > 0)
		#stop("The comma-separated snv table requires the columns date and dna_profile!")
		
		# Assert that that the dates are in the correct format
		try:
			seq_info_table["date"].apply(datetime.date.fromisoformat) 
		except ValueError:
			raise ValueError("Error in date format: Please provide dates in the format yyyy-mm-dd!")
		#if(class(try_parsing) == "try-error")
		#stop(try_parsing, "\nPlease provide dates in the format yyyy-mm-dd!")

		
		### extract important columns 
		minDate = datetime.date.fromisoformat(min(seq_info_table['date']))
		
		#tic("make small table")
		print("Extract columns\n")
		## 24 Sekunden aufm Mac
		seq_info_short_table = seq_info_table.rename(columns = {"dna_profile": "snvs"})[['date', 'snvs']]
		seq_info_short_table['t'] = (seq_info_table["date"].apply(datetime.date.fromisoformat) - minDate).dt.days
		seq_info_short_table = seq_info_short_table.sort_values("t")

	except Exception as e:
		print(e)
	else:
	  return(seq_info_short_table)
