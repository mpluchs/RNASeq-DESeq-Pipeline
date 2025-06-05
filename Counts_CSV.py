import os
import pandas as pd
from pathlib import Path
import argparse

def convert_to_csv(file_path):
    """
    Convert a featureCounts output file from txt to csv, renaming the count column.
    """
    try:
        # Skip the first comment line and read the table
        df = pd.read_csv(file_path, sep='\t', skiprows=1)
        
        # Rename the last column to simply "Count"
        df.columns = list(df.columns[:-1]) + ['Count']
        
        # Create output filename
        output_path = file_path.with_suffix('.csv')
        
        # Save to CSV
        df.to_csv(output_path, index=False)
        print(f"Converted {file_path.name} -> {output_path.name}")
        
    except Exception as e:
        print(f"Error processing {file_path}: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='Convert featureCounts txt files to csv format.')
    parser.add_argument('input_dir', help='Directory containing featureCounts output .txt files')
    
    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    
    # Process all .txt files in the directory
    count_files = list(input_dir.glob('*_counts.txt'))
    
    if not count_files:
        print(f"No count files found in {input_dir}")
        return
    
    for file_path in count_files:
        convert_to_csv(file_path)

if __name__ == "__main__":
    main()
