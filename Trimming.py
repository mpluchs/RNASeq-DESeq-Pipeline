import os
import subprocess
import glob

def trim_and_fastqc(folder_path, output_dir):
    # Find all R1 and R2 files in the folder
    r1_files = glob.glob(os.path.join(folder_path, '*_R1_001_repair.fastq.gz'))
    r2_files = glob.glob(os.path.join(folder_path, '*_R2_001_repair.fastq.gz'))

    # Sort the files to ensure they are in the same order for pairing
    r1_files.sort()
    r2_files.sort()

    # Check if the number of R1 and R2 files match
    if len(r1_files) != len(r2_files):
        print("Error: Number of R1 and R2 files do not match.")
        return

    # Iterate over each pair of R1 and R2 files
    for r1, r2 in zip(r1_files, r2_files):
        # Define output file names for trim_galore
            # Define output file names for trim_galore
        base_name = os.path.basename(r1).replace('_R1_001_repair.fastq.gz', '')
        trimmed_r1_out = os.path.join(output_dir, f"{base_name}_R1_trimmed.fastq.gz")
        trimmed_r2_out = os.path.join(output_dir, f"{base_name}_R2_trimmed.fastq.gz")
        fastqc_report_out = os.path.join(output_dir, f"{base_name}_Fastqc_Report")

        # Run trim_galore
        trim_cmd = f"trim_galore --fastqc --paired --length 20 --clip_R2 15 --three_prime_clip_R1 15 -o {fastqc_report_out} {r1} {r2}"
        subprocess.run(trim_cmd, shell=True, check=True)

def check_trimming():
    folder_path = input("Enter the folder path containing your data: ").strip()
    output_dir = input("Enter the output directory for trimmed files: ").strip()
    trim_and_fastqc(folder_path, output_dir)  # Call your trim_and_fastqc function

# Call the function to start the process
check_trimming()
