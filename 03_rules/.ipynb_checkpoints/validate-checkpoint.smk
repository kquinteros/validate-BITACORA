###--- create output directories ---###
rule directories:
    input: 
        "05_output"
    params:
        s_start= config['seq_start'],
        s_end= config['seq_end'],
        p_start= config['pro_start'],
        p_end= config['pro_end']
    output:
        "05_output/rule-directories.out"
    message:
        "Creating directories for interproscan and validation steps"
    shell: 
        """
        ./04_scripts/make_directory.sh {input} {params.s_start} {params.s_end} {params.p_start} {params.p_end} > {output}
        """
        
###--- run interproscan with pfam ---###
rule run_interproscan:
    input:
        rule_directory= "05_output/rule-directories.out",
        fasta= "00_data/{sam}/{pro}/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered.fasta"
    output:
        tsv= "05_output/{sam}/{pro}/interproscan/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered.fasta.tsv"
    params:
        format= "TSV"
    log: 
        "05_output/{sam}/{pro}/interproscan/interproscan.log"
    threads: 10
    message:
        "Running InterProScan on predicted protein sequences using the provided Bash script."
    shell:
        """
        interproscan -app pfam -cpu {threads} -i {input.fasta} -f {params.format} --outfile {output.tsv} > {log}
        """

###--- Filter protein by protein domains ---###
rule run_filter_protein_domain:
    input:
        fasta= "00_data/{sam}/{pro}/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered.fasta",
        interproscan= "05_output/{sam}/{pro}/interproscan/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered.fasta.tsv"
    output:
        "05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_R1.fasta"
    params:
        acc_ids= lambda wildcards: '"{}"'.format(protein_pfam_input(wildcards)),
        script = "04_scripts/filter-protein-script.R"
    conda:
        config['env-r']
    message:
        "Using custom R script that takes InterProScan results to filter out protein sequences lacking the appropriate PFAM protein domain."
    shell:
        """
        Rscript {params.script} --input {input.fasta} --output {output} --interproscan {input.interproscan} --acc_ids {params.acc_ids}
        """

###--- Filter protein by protein sequence length size ---###
rule run_filter_size_selection:
    input:
        "05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_R1.fasta"
    output:
        "05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_R2.fasta"
    params:
        min = lambda wildcards: protein_size_min(wildcards)["min"],
        max = lambda wildcards: protein_size_min(wildcards)["max"],
        script = "04_scripts/filter-size-selection-script.R"
    conda:
        config['env-r']
    message:
        "Using custom R script that executes a function to filter FASTA sequences within a specified range."
    shell:
        """
        Rscript {params.script} --input {input} --output {output} --min {params.min} --max {params.max}
        """
###--- run TMbed to evaluate transmembrane domains ---###
rule run_TMbed:
    input:
        "05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_R2.fasta"
    output:
        embed = "05_output/{sam}/{pro}/TMbed/{pro}_embeddings.h5",
        pred = "05_output/{sam}/{pro}/TMbed/{pro}_transmembrane.pred"
    resources:
        nvidia_gpu=1
    threads: 20
    conda:
        config['env-python']
    message:
        "Run TMbed, this program predicts transmembrane protein domains"
    shell:
        """
        # Set CUDA_LAUNCH_BLOCKING
        export CUDA_LAUNCH_BLOCKING=1

        # Set CUBLAS_WORKSPACE_CONFIG
        export CUBLAS_WORKSPACE_CONFIG=:16:8

        # Check if CUDA device is available
        if command -v nvidia-smi &> /dev/null
        then
            echo "CUDA device found. Running with GPU."
            gpu_option="--use-gpu"
        else
            echo "No CUDA device found. Running without GPU."
            gpu_option="--no-use-gpu"
        fi

        # Generate embeddings for a set of protein sequences
        python -m tmbed embed -f {input} -e {output.embed} $gpu_option --cpu-fallback 

        # Predict transmembrane proteins and segments
        python -m tmbed predict -f {input} -e {output.embed} -p {output.pred} $gpu_option --cpu-fallback --out-format=0
        """
        
###--- run AWK script to count transmembrane domains from TMBED output---###
rule run_number_of_trans_domains:
    input:
        "05_output/{sam}/{pro}/TMbed/{pro}_transmembrane.pred"
    output:
        "05_output/{sam}/{pro}/TMbed/{pro}_transmembrane_counts.csv"
    message:
        "The script uses AWK to count occurrences of transmembrane domains in protein sequences, associating counts with sequence headers and outputting the results in CSV format. The input file is a .pred from TMbed"
    shell:
        """
        awk 'NR % 3 == 0 {{count = gsub(/[Hh]+/, "&", $0); printf "%s,%s\\n", header, count}} /^>/ {{header = $0}}' {input} > {output}
        """

###--- run R script that filters protein sequences based ---###
rule run_filter_trans_domains:
    input:
        fasta="05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_R2.fasta",
        preds="05_output/{sam}/{pro}/TMbed/{pro}_transmembrane_counts.csv"
    output:
        target="05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_R3.target.fasta"
    params:
        num=7,
        script="04_scripts/filter-trans-domains-script.R"
    conda:
        config['env-r']
    message:
        "Using custom R script that filters protein sequences based on the number of transmembrane domains"
    shell:
        """
        Rscript {params.script} --fasta {input.fasta} --pred {input.preds} --target {output.target} --number {params.num}
        """

###--- Exlude sequences in gff not found in fasta file ---###
rule exclude_sequences:
    input:
        fasta="05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_R3.target.fasta",
        gff="00_data/{sam}/{pro}/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered.gff3"
    output:
        "05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_target.gff3"
    params:
        script="bitacora/Scripts/Tools/exclude_sequences_ingff_notinfasta.pl"
    conda:
        config['env-genomics']
    message:
        "Create a GFF file that contains only our validate sequences."
    shell:
        """
        perl {params.script} {input.fasta} {input.gff} {output}
        """
###--- Ensure that fasta files and GFF have the same name for features and sequences  ---###
rule ensure_feature_names:
    input:
        fasta= "00_data/{sam}/{sam}.fasta",
        gff= "05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_excluded.gff3"
    output:
        cds = "05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_final.cds.fasta",
        pep = "05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_final.pep.fasta"
    params:
        script= "bitacora/Scripts/Tools/gff2fasta_v3.pl",
        prefix = "05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_final"
    conda:
        config['env-genomics']
    message:
        "Ensure that fasta files and GFF have the same name for features and sequences"
    shell:
        """
        perl {params.script} {input.fasta} {input.gff} {params.prefix}
        """
rule run_filter_overlapping:
    input:
        fasta = "05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_final.pep.fasta",
        gff = "05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_excluded.gff3"
    output:
        "05_output/{sam}/{pro}/validated/{pro}_genomic_and_annotated_proteins_trimmed_idseqsclustered_excluded_final.gff3"
    conda:
        config['env-genomics']
    message:
        "This script fixes and/or standardizes any GTF/GFF file into full sorted GTF/GFF file. It AGAT parser removes duplicate features, fixes duplicated IDs, adds missing ID and/or Parent attributes, deflates factorized attributes (attributes with several parents are duplicated with uniq ID), add missing features when possible (e.g. add exon if only CDS described, add UTR if CDS and exon described), fix feature locations"
    shell:
        """
        agat_convert_sp_gxf2gxf.pl --gff {input.gff} -o {output}
        """