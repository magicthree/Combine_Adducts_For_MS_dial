import os
import pandas as pd

from method import extract_file, rm_dup, blank_filter, merge_file, peak_finder, alignment

file_path = r""
file_name = ''
with_head = False
file_data = pd.read_csv(os.path.join(file_path, file_name), sep='\t')

qc_columns = file_data.iloc[0] == 'QC'
desired_columns = file_data.columns[qc_columns]
QC = file_data.loc[3, desired_columns].tolist()
bk_columns = file_data.iloc[0] == 'Blank'
desired_columns = file_data.columns[bk_columns]
BK = file_data.loc[3, desired_columns].tolist()

extract_file(file_name, file_path, with_head)
rm_dup(file_name[:-4] + '-key.csv', file_path, 0.2, 0.02)
alignment(file_name[:-4] + '-key-filter.csv', file_path, RTolerance=0.2)
peak_finder(file_name=file_name[:-4] + '-key-filter-aligned.csv',
            file_unknown_name=file_name[:-4] + '-key-filter_with_unknown.csv', file_path=file_path)
blank_filter(file_name[:-4] + '-key-filter-aligned-find.csv', file_path, 3000, 10, QC, BK)
merge_file(file_name[:-4] + '-key-filter-aligned-find-bk.csv', file_path, QC, BK)
