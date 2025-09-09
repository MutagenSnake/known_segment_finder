import os
import subprocess

# Base folder containing subfolders with fasta files
base_folder = 'data/cluster_alignments/'

# Set this to True to run only one file (for testing)
# test_mode = True

for root, dirs, files in os.walk(base_folder):
    for f in files:
        if not f.endswith('.fasta'):
            continue
        
        file_path = os.path.join(root, f)
        treefile = os.path.join(root, f + ".treefile")
        
        if os.path.exists(treefile):
            print(f"Skipping {f}, treefile already exists.")
            continue
        else:
            print(f"Running IQ-TREE on {f}...")
            cmd = [
                "iqtree",
                "-s", file_path,
                "-m", "GTR+G",
                "-bb", "1000",
                "-alrt", "1000",
                "-nt", "AUTO"
            ]
            subprocess.run(cmd)
        
        # # If in test mode, run only on the first file
        # if test_mode:
        #     print("Test mode active: stopping after first file.")
        #     exit()