import os
import subprocess

def automated_alignment(input_folder, output_folder):
    # Walk through the input folder and its subfolders
    for root, dirs, files in os.walk(input_folder):
        for filename in files:
            if filename.endswith('_R1_001_val_1.fq.gz'):
                # Extract sample name
                sample_name = filename.split('_R1_001_val_1.fq.gz')[0]

                # Check if corresponding R2 file exists
                r2_filename = sample_name + '_R2_001_val_2.fq.gz'
                r2_path = os.path.join(root, r2_filename)
                if os.path.exists(r2_path):
                    # Build command
                    command = f'hisat2 -t --rna-strandness RF --summary-file {output_folder}/{sample_name}.txt -p 24 -x /Users/miseq/Desktop/Luchs/Reference_Genome2.0/genome/genome_index -1 {os.path.join(root, filename)} -2 {r2_path} | samtools view -@ 24 -Shu - | samtools sort -@ 24 -n -o {output_folder}/{sample_name}.bam'

                    # Print the command before executing
                    print(f"Executing command: {command}")

                    # Execute command and display output
                    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                    for line in process.stdout:
                        print(line.decode().strip())

def check_alignment():
    folder_path = input("Enter the folder path containing your data: ").strip()
    output_dir = input("Enter the desired output directory: ").strip()
    automated_alignment(folder_path, output_dir)  # Call your trim_and_fastqc function

# Call the function to start the process
check_alignment()
