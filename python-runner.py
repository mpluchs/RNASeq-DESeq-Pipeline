#!/usr/bin/env python3

"""
DESeq2 Comparison Runner

This script runs multiple DESeq2 comparisons by calling the R script for each comparison.
It manages the execution of all specified comparisons and tracks their progress.
"""

import os
import subprocess
import argparse
import time
import concurrent.futures
from datetime import datetime

# List of all comparisons to run
COMPARISONS = [
    # S18_11.5 comparisons
    {"name": "S18_11C_vs_S18_11R", "group1": "S18_11C", "group2": "S18_11R"},
    {"name": "S18_11C_vs_S18_Sib1", "group1": "S18_11C", "group2": "S18_Sib1"},
    {"name": "S18_11C_vs_S18_Sib2", "group1": "S18_11C", "group2": "S18_Sib2"},
    {"name": "S18_11C_vs_Selfies_11-18", "group1": "S18_11C", "group2": "Selfies_11-18"},
    {"name": "S18_11R_vs_S18_Sib1", "group1": "S18_11R", "group2": "S18_Sib1"},
    {"name": "S18_11R_vs_S18_Sib2", "group1": "S18_11R", "group2": "S18_Sib2"},
    {"name": "S18_11R_vs_Selfies_11-18", "group1": "S18_11R", "group2": "Selfies_11-18"},
    {"name": "S18_Sib1_vs_S18_Sib2", "group1": "S18_Sib1", "group2": "S18_Sib2"},
    {"name": "S18_Sib1_vs_Selfies_11-18", "group1": "S18_Sib1", "group2": "Selfies_11-18"},
    {"name": "S18_Sib2_vs_Selfies_11-18", "group1": "S18_Sib2", "group2": "Selfies_11-18"},
    
    # S18_12.5 comparisons
    {"name": "S18_12C_vs_S18_12R", "group1": "S18_12C", "group2": "S18_12R"},
    {"name": "S18_12C_vs_S18_Sib1", "group1": "S18_12C", "group2": "S18_Sib1"},
    {"name": "S18_12C_vs_S18_Sib2", "group1": "S18_12C", "group2": "S18_Sib2"},
    {"name": "S18_12C_vs_Selfies_12-18", "group1": "S18_12C", "group2": "Selfies_12-18"},
    {"name": "S18_12R_vs_S18_Sib1", "group1": "S18_12R", "group2": "S18_Sib1"},
    {"name": "S18_12R_vs_S18_Sib2", "group1": "S18_12R", "group2": "S18_Sib2"},
    {"name": "S18_12R_vs_Selfies_12-18", "group1": "S18_12R", "group2": "Selfies_12-18"},
    {"name": "S18_Sib1_vs_Selfies_12-18", "group1": "S18_Sib1", "group2": "Selfies_12-18"},
    {"name": "S18_Sib2_vs_Selfies_12-18", "group1": "S18_Sib2", "group2": "Selfies_12-18"},
    
    # S30_11.5 comparisons
    {"name": "S30_11C_vs_S30_11R", "group1": "S30_11C", "group2": "S30_11R"},
    {"name": "S30_11C_vs_S30_Sib1", "group1": "S30_11C", "group2": "S30_Sib1"},
    {"name": "S30_11C_vs_S30_Sib2", "group1": "S30_11C", "group2": "S30_Sib2"},
    {"name": "S30_11C_vs_Selfies_11-30", "group1": "S30_11C", "group2": "Selfies_11-30"},
    {"name": "S30_11R_vs_S30_Sib1", "group1": "S30_11R", "group2": "S30_Sib1"},
    {"name": "S30_11R_vs_S30_Sib2", "group1": "S30_11R", "group2": "S30_Sib2"},
    {"name": "S30_11R_vs_Selfies_11-30", "group1": "S30_11R", "group2": "Selfies_11-30"},
    {"name": "S30_Sib1_vs_S30_Sib2", "group1": "S30_Sib1", "group2": "S30_Sib2"},
    {"name": "S30_Sib1_vs_Selfies_11-30", "group1": "S30_Sib1", "group2": "Selfies_11-30"},
    {"name": "S30_Sib2_vs_Selfies_11-30", "group1": "S30_Sib2", "group2": "Selfies_11-30"},
    
    # S30_12.5 comparisons
    {"name": "S30_12C_vs_S30_12R", "group1": "S30_12C", "group2": "S30_12R"},
    {"name": "S30_12C_vs_S30_Sib1", "group1": "S30_12C", "group2": "S30_Sib1"},
    {"name": "S30_12C_vs_S30_Sib2", "group1": "S30_12C", "group2": "S30_Sib2"},
    {"name": "S30_12C_vs_Selfies_12-30", "group1": "S30_12C", "group2": "Selfies_12-30"},
    {"name": "S30_12R_vs_S30_Sib1", "group1": "S30_12R", "group2": "S30_Sib1"},
    {"name": "S30_12R_vs_S30_Sib2", "group1": "S30_12R", "group2": "S30_Sib2"},
    {"name": "S30_12R_vs_Selfies_12-30", "group1": "S30_12R", "group2": "Selfies_12-30"},
    {"name": "S30_Sib1_vs_Selfies_12-30", "group1": "S30_Sib1", "group2": "Selfies_12-30"},
    {"name": "S30_Sib2_vs_Selfies_12-30", "group1": "S30_Sib2", "group2": "Selfies_12-30"},
    
    # Treatment Stage Comparisons
    {"name": "S18_11C_vs_S18_12C", "group1": "S18_11C", "group2": "S18_12C"},
    {"name": "S18_11R_vs_S18_12R", "group1": "S18_11R", "group2": "S18_12R"},
    {"name": "Selfies_11-18_vs_Selfies_12-18", "group1": "Selfies_11-18", "group2": "Selfies_12-18"},
    {"name": "S30_11C_vs_S30_12C", "group1": "S30_11C", "group2": "S30_12C"},
    {"name": "S30_11R_vs_S30_12R", "group1": "S30_11R", "group2": "S30_12R"},
    {"name": "Selfies_11-30_vs_Selfies_12-30", "group1": "Selfies_11-30", "group2": "Selfies_12-30"},
    
    # Fixation Stage Comparisons
    {"name": "S18_11C_vs_S30_11C", "group1": "S18_11C", "group2": "S30_11C"},
    {"name": "S18_11R_vs_S30_11R", "group1": "S18_11R", "group2": "S30_11R"},
    {"name": "Selfies_11-18_vs_Selfies_11-30", "group1": "Selfies_11-18", "group2": "Selfies_11-30"},
    {"name": "S18_12C_vs_S30_12C", "group1": "S18_12C", "group2": "S30_12C"},
    {"name": "S18_12R_vs_S30_12R", "group1": "S18_12R", "group2": "S30_12R"},
    {"name": "Selfies_12-18_vs_Selfies_12-30", "group1": "Selfies_12-18", "group2": "Selfies_12-30"}
]

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Run multiple DESeq2 comparisons')
    parser.add_argument('--base-dir', required=True, help='Base directory containing all sample directories')
    parser.add_argument('--output-dir', required=True, help='Directory for output files')
    parser.add_argument('--r-script', required=True, help='Path to the R script')
    parser.add_argument('--parallel', type=int, default=1, help='Number of comparisons to run in parallel (default: 1)')
    parser.add_argument('--selected', nargs='+', help='Only run specified comparisons (by name)')
    parser.add_argument('--log-file', default='deseq2_runner.log', help='Log file path')
    return parser.parse_args()

def setup_logging(log_file):
    """Create a log file for the run"""
    os.makedirs(os.path.dirname(os.path.abspath(log_file)), exist_ok=True)
    with open(log_file, 'w') as f:
        f.write(f"DESeq2 Comparison Runner\n")
        f.write(f"Started at: {datetime.now()}\n")
        f.write(f"{'='*50}\n\n")

def log_message(log_file, message):
    """Write a message to the log file"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, 'a') as f:
        f.write(f"[{timestamp}] {message}\n")
    print(f"[{timestamp}] {message}")

def run_comparison(r_script, base_dir, output_dir, comparison, log_file):
    """Run a single DESeq2 comparison"""
    group1_dir = os.path.join(base_dir, comparison["group1"])
    group2_dir = os.path.join(base_dir, comparison["group2"])
    comp_output_dir = os.path.join(output_dir, comparison["name"])
    
    # Create output directory if it doesn't exist
    os.makedirs(comp_output_dir, exist_ok=True)
    
    # Construct command
    cmd = [
        "Rscript", 
        r_script, 
        group1_dir, 
        group2_dir, 
        comp_output_dir, 
        comparison["name"]
    ]
    
    # Log the command
    log_message(log_file, f"Running: {' '.join(cmd)}")
    
    # Run the command
    start_time = time.time()
    try:
        # Redirect stdout and stderr to log files in the output directory
        stdout_log = os.path.join(comp_output_dir, "stdout.log")
        stderr_log = os.path.join(comp_output_dir, "stderr.log")
        
        with open(stdout_log, 'w') as stdout_file, open(stderr_log, 'w') as stderr_file:
            process = subprocess.Popen(
                cmd, 
                stdout=stdout_file, 
                stderr=stderr_file
            )
            process.wait()
        
        # Check if process completed successfully
        if process.returncode == 0:
            elapsed_time = time.time() - start_time
            log_message(log_file, f"✓ Completed {comparison['name']} in {elapsed_time:.2f} seconds")
            return True
        else:
            with open(stderr_log, 'r') as f:
                error_message = f.read()
            log_message(log_file, f"✗ Failed {comparison['name']}: {error_message}")
            return False
            
    except Exception as e:
        elapsed_time = time.time() - start_time
        log_message(log_file, f"✗ Error in {comparison['name']}: {str(e)} after {elapsed_time:.2f} seconds")
        return False

def create_summary_report(output_dir, comparisons_run, comparisons_succeeded, log_file):
    """Create a summary report of all comparisons"""
    summary_file = os.path.join(output_dir, "deseq2_summary_report.html")
    
    # Gather statistics for successful comparisons
    comparison_stats = []
    for comp in comparisons_succeeded:
        comp_dir = os.path.join(output_dir, comp["name"])
        summary_path = os.path.join(comp_dir, f"{comp['name']}_summary.txt")
        
        # Extract key statistics from summary file
        stats = {"name": comp["name"], "group1": comp["group1"], "group2": comp["group2"]}
        
        if os.path.exists(summary_path):
            with open(summary_path, 'r') as f:
                for line in f:
                    if "Total genes analyzed:" in line:
                        stats["total_genes"] = line.split(":")[-1].strip()
                    elif "Significant genes (FDR < 0.05):" in line:
                        stats["sig_genes"] = line.split(":")[-1].strip()
                    elif "Up-regulated genes:" in line:
                        stats["up_genes"] = line.split(":")[-1].strip()
                    elif "Down-regulated genes:" in line:
                        stats["down_genes"] = line.split(":")[-1].strip()
        
        comparison_stats.append(stats)
    
    # Create HTML report
    with open(summary_file, 'w') as f:
        f.write(f"""<!DOCTYPE html>
<html>
<head>
    <title>DESeq2 Analysis Summary Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 20px;
        }}
        h1, h2, h3 {{
            color: #2c3e50;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin-bottom: 20px;
        }}
        th, td {{
            text-align: left;
            padding: 8px;
            border: 1px solid #ddd;
        }}
        th {{
            background-color: #f2f2f2;
        }}
        tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}
        .success {{
            color: green;
        }}
        .failure {{
            color: red;
        }}
        .summary {{
            margin: 20px 0;
            padding: 10px;
            background-color: #f8f9fa;
            border-radius: 5px;
        }}
    </style>
</head>
<body>
    <h1>DESeq2 Analysis Summary Report</h1>
    <p>Generated on: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
    
    <div class="summary">
        <h2>Run Summary</h2>
        <p>Total comparisons: {len(comparisons_run)}</p>
        <p>Successful: <span class="success">{len(comparisons_succeeded)}</span></p>
        <p>Failed: <span class="failure">{len(comparisons_run) - len(comparisons_succeeded)}</span></p>
    </div>
    
    <h2>Comparison Results</h2>
    <table>
        <tr>
            <th>Comparison</th>
            <th>Group 1</th>
            <th>Group 2</th>
            <th>Total Genes</th>
            <th>Significant Genes</th>
            <th>Up-regulated</th>
            <th>Down-regulated</th>
            <th>Status</th>
        </tr>
""")
        
        # Add rows for each comparison
        for comp in comparisons_run:
            # Find if this comparison succeeded
            succeeded = any(s["name"] == comp["name"] for s in comparisons_succeeded)
            
            # Find stats if available
            stats = next((s for s in comparison_stats if s["name"] == comp["name"]), 
                        {"total_genes": "N/A", "sig_genes": "N/A", "up_genes": "N/A", "down_genes": "N/A"})
            
            status_class = "success" if succeeded else "failure"
            status_text = "Success" if succeeded else "Failed"
            
            f.write(f"""
        <tr>
            <td>{comp["name"]}</td>
            <td>{comp["group1"]}</td>
            <td>{comp["group2"]}</td>
            <td>{stats.get("total_genes", "N/A")}</td>
            <td>{stats.get("sig_genes", "N/A")}</td>
            <td>{stats.get("up_genes", "N/A")}</td>
            <td>{stats.get("down_genes", "N/A")}</td>
            <td class="{status_class}">{status_text}</td>
        </tr>""")
        
        f.write("""
    </table>
    
    <h2>Links to Individual Reports</h2>
    <ul>
""")
        
        # Add links to individual comparison reports
        for comp in comparisons_succeeded:
            f.write(f"""
        <li><a href="{comp['name']}/{comp['name']}_DESeq2.csv">{comp['name']} Results</a> 
            (<a href="{comp['name']}/{comp['name']}_MA_plot.pdf">MA Plot</a> | 
             <a href="{comp['name']}/{comp['name']}_pvalue_hist.pdf">p-value Histogram</a> | 
             <a href="{comp['name']}/{comp['name']}_PCA_plot.pdf">PCA Plot</a>)</li>""")
        
        f.write("""
    </ul>
</body>
</html>""")
    
    log_message(log_file, f"Summary report created: {summary_file}")
    return summary_file

def main():
    """Main function to run multiple DESeq2 comparisons"""
    args = parse_arguments()
    
    # Set up paths
    base_dir = os.path.abspath(args.base_dir)
    output_dir = os.path.abspath(args.output_dir)
    r_script = os.path.abspath(args.r_script)
    log_file = os.path.abspath(args.log_file)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Set up logging
    setup_logging(log_file)
    log_message(log_file, f"Starting DESeq2 analysis with {args.parallel} parallel processes")
    log_message(log_file, f"Base directory: {base_dir}")
    log_message(log_file, f"Output directory: {output_dir}")
    log_message(log_file, f"R script: {r_script}")
    
    # Filter comparisons if selected
    comparisons_to_run = COMPARISONS
    if args.selected:
        selected_set = set(args.selected)
        comparisons_to_run = [comp for comp in COMPARISONS if comp['name'] in selected_set]
        log_message(log_file, f"Running {len(comparisons_to_run)} selected comparisons")
    else:
        log_message(log_file, f"Running all {len(comparisons_to_run)} comparisons")
    
    # Run comparisons in parallel
    comparisons_succeeded = []
    
    if args.parallel > 1:
        log_message(log_file, f"Running comparisons in parallel with {args.parallel} workers")
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.parallel) as executor:
            future_to_comp = {
                executor.submit(
                    run_comparison, 
                    r_script, 
                    base_dir, 
                    output_dir, 
                    comp, 
                    log_file
                ): comp for comp in comparisons_to_run
            }
            
            for future in concurrent.futures.as_completed(future_to_comp):
                comp = future_to_comp[future]
                try:
                    success = future.result()
                    if success:
                        comparisons_succeeded.append(comp)
                except Exception as e:
                    log_message(log_file, f"Exception in {comp['name']}: {str(e)}")
    else:
        log_message(log_file, "Running comparisons sequentially")
        for comp in comparisons_to_run:
            success = run_comparison(r_script, base_dir, output_dir, comp, log_file)
            if success:
                comparisons_succeeded.append(comp)
    
    # Create summary report
    log_message(log_file, f"Completed {len(comparisons_succeeded)} of {len(comparisons_to_run)} comparisons successfully")
    summary_file = create_summary_report(output_dir, comparisons_to_run, comparisons_succeeded, log_file)
    
    # Final message
    log_message(log_file, f"Analysis complete! Summary report: {summary_file}")

if __name__ == "__main__":
    main()
