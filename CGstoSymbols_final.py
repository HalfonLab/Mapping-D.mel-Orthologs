#!/usr/bin/env python

def create_FbgnDict(fbgn_annotation_file):
	FbgnDict = {}
	with open(fbgn_annotation_file, 'r') as f:
		for line in f:
			if not line.startswith('#') and line.strip():
				columns = line.strip().split('\t')
				symbol = columns[0]  # Assuming symbols is the 3rd column
				annotation_ID = columns[4]  # Assuming annotation_ID is the 5th column
				FbgnDict[annotation_ID] = symbol
	return FbgnDict

def process_second_file(input_file, output_file, FbgnDict):
	with open(input_file, 'r') as f:
		lines = f.readlines()

	with open(output_file, 'w') as f:
		for line in lines:
			columns = line.strip().split('\t')
			if len(columns) >= 18:
				annotation_ids_7th_col = columns[6].split(',')  # Split multiple IDs by comma
				annotation_ids_12th_col = columns[11].split(',')  # Split multiple IDs by comma
				
				new_ids_7th_col = [FbgnDict.get(annot_id, annot_id) for annot_id in annotation_ids_7th_col]
				new_ids_12th_col = [FbgnDict.get(annot_id, annot_id) for annot_id in annotation_ids_12th_col]
				
				columns[6] = ','.join(new_ids_7th_col)
				columns[11] = ','.join(new_ids_12th_col)
				
				new_line = '\t'.join(columns)
				f.write(new_line + '\n')

def main():
	fbgn_annotation_file = input("Enter the path to the first file (fbgn_annotation_ID_fb_2023_04.tsv): ")
	input_file = input("Enter the path to the SCRM file containing CGs(second_file.scrm_3.bed): ")
	output_file = input("Enter the desired output file name with symbols (output_aceph_symb.scrm.bed): ")
	
	FbgnDict = create_FbgnDict(fbgn_annotation_file)
	process_second_file(input_file, output_file, FbgnDict)
	
	print("Processing completed.")

if __name__ == "__main__":
	main()
