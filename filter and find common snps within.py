import os
import pandas as pd

input = r"C:\Users\sayeraselvan\Desktop\output1\bio"
output = r"C:\Users\sayeraselvan\Desktop\output1\bio\filter"

def apply_threshold(file_path, threshold_value, chrom):
    df = pd.read_csv(file_path, sep='\t')
    column_name = 'p_wald'
    chrom_col = 'chr'
    mask = df[column_name] < threshold_value
    filtered_df = df[mask]

    output_filename = f"pos_1e7_{os.path.basename(file_path)}"
    output_file = os.path.join(output, output_filename)
    filtered_df.to_csv(output_file, sep='\t', index=False)
'''    
    chromo = filtered_df[chrom_col] == chrom
    chromosome = filtered_df[chromo]
    chrom_output_filename = f"chr_pthreshold_{os.path.basename(file_path)}"
    chrom_output_file = os.path.join(output, chrom_output_filename)
    chromosome.to_csv(chrom_output_file, sep='\t', index=False)

'''


threshold_value = 1e-7
chrom = 1

# Create the output directory if it doesn't exist
os.makedirs(output, exist_ok=True)

for filename in os.listdir(input):
    if filename.endswith(".assoc.txt"):
        file_path = os.path.join(input, filename)
        apply_threshold(file_path, threshold_value, chrom)
