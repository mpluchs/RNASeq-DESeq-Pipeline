import os
import subprocess
import glob
import sys

def check_file_exists(filepath, description):
    """Verify if a file or directory exists"""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"{description} not found: {filepath}")

def ensure_directory_exists(directory):
    """Create directory if it doesn't exist"""
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
            print(f"Created output directory: {directory}")
        except Exception as e:
            raise RuntimeError(f"Failed to create directory {directory}: {str(e)}")

def run_feature_counts(gtf_file, input_folder, output_folder):
    # Verify input files and directories exist
    check_file_exists(gtf_file, "GTF file")
    check_file_exists(input_folder, "Input folder")
    ensure_directory_exists(output_folder)
    
    # Keep track of successes and failures
    successful_runs = []
    failed_runs = []
    
    # Find all BAM files
    bam_files = [f for f in os.listdir(input_folder) if f.endswith('.bam')]
    if not bam_files:
        raise RuntimeError(f"No BAM files found in {input_folder}")
    
    print(f"Found {len(bam_files)} BAM files to process")
    
    # Iterate through files in the input folder
    for filename in bam_files:
        # Extract sample name
        sample_name = filename.split('.bam')[0]
        # Construct full paths
        bam_file = os.path.join(input_folder, filename)
        output_file = os.path.join(output_folder, f"{sample_name}_counts.txt")
        
        # Verify BAM file exists and is not empty
        if os.path.getsize(bam_file) == 0:
            failed_runs.append((sample_name, "BAM file is empty"))
            print(f"Warning: {bam_file} is empty, skipping...")
            continue
        
        # Build command for featureCounts
        command = [
            'featureCounts', '-a', gtf_file, '-o', output_file,
            '-T', '2', '-s', '2', '-Q', '0', '-t', 'gene', '-g', 'gene_id',
            '--minOverlap', '1', '--fracOverlap', '0', '--fracOverlapFeature', '0',
            '-p', '-C', bam_file
        ]
        
        try:
            # Execute command and display output
            print(f"\nProcessing {sample_name}...")
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True  # Get string output instead of bytes
            )
            
            # Capture stdout and stderr
            stdout, stderr = process.communicate()
            
            # Check if process was successful
            if process.returncode == 0:
                successful_runs.append(sample_name)
                print(f"✓ Successfully processed {sample_name}")
                
                # Verify output file was created and has content
                if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                    print(f"  Output file created: {output_file}")
                else:
                    raise RuntimeError("Output file is empty or wasn't created")
                    
            else:
                failed_runs.append((sample_name, stderr))
                print(f"✗ Failed to process {sample_name}")
                print(f"Error: {stderr}")
                
        except Exception as e:
            failed_runs.append((sample_name, str(e)))
            print(f"✗ Error processing {sample_name}: {str(e)}")
    
    # Print summary
    print("\n=== Pipeline Summary ===")
    print(f"Total files processed: {len(bam_files)}")
    print(f"Successful: {len(successful_runs)}")
    print(f"Failed: {len(failed_runs)}")
    
    if failed_runs:
        print("\nFailed samples:")
        for sample, error in failed_runs:
            print(f"- {sample}: {error}")
    
    return successful_runs, failed_runs

def fcounts_pipeline():
    try:
        gtf_location = "/Volumes/CaChannel/Luchs/Reference_gtf/genomic.gtf"
        folder_path = input("Enter the folder path containing your data: ").strip()
        output_dir = input("Enter the desired output directory: ").strip()
        
        successful, failed = run_feature_counts(gtf_location, folder_path, output_dir)
        
        if failed:
            sys.exit(1)  # Exit with error code if any runs failed
            
    except KeyboardInterrupt:
        print("\nPipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nPipeline failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    fcounts_pipeline()
