# Analyzing SLiMs Present in Viral nsP3 #

## First, put the following scripts into your 'modules' folder for 'ELM_motif_query' in 'my_modules': 

    ELM_motif_query_automation_filtering.py,
    filtering_results.py,
    motif_matrix.py,
    motif_conservation.py

## Running ELM Database motif query on your protein:

### Ensure that your nsP3 FASTA files are in the 'input_folder'.

### Then, run the following command:
    python "my_modules/ELM_motif_query/modules/ELM_motif_query_automation_filtering.py" --filter

#### This process may take time, but should be continually updating in the terminal.

#### After this is done, look at your results for each sequence and see how many motifs were found in either disordered or all regions! There are thousands for each!

## Filtering motif hits to unique motifs within disordered regions:

### Copy-paste your previous results into the 'input_folder'.

### Then, run the following command:

    python "my_modules/ELM_motif_query/modules/filtering_results.py"

#### Check your outputs after the process is completed in the terminal. There should be an additional file for each sequence that ends in "_unique.csv", indicating that each motif is only listed once per instance. Check how many hits are in individual sequences now.

## Combine your results:

### Clear out your 'input_folder', then copy-paste all of your results files in the 'output_folder' into the 'input_folder'. 

#### At this point, you can save a copy of all of this data in a separate location and clear out your 'output_folder' as well.

### Then, run the following command:

    python "my_modules/ELM_motif_query/modules/motif_matrix.py"

#### Look at the matrix that indicates whether each sequence (row) has a specific motif, which is listed as either absent (0) or present (1).
#### Save this matrix file, it is very important!

## Lastly, visualize your data to spot trends by eye:

### Clear the input and output folders after saving appropriate data files. Move your motif matrix file to the 'input_folder'.

### Run the following command to create bar charts to visualize your data:

    python "my_modules/ELM_motif_query/modules/motif_conservation.py"