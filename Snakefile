import pandas as pd

# ---- Config ----
SAMPLESHEET = "sample_sheet.csv"
OUTDIR = "data/SRA"
MAX_SRA = 20   # change this to limit total number of SRAs, None = all

# ---- Parse sample sheet ----
df = pd.read_csv(SAMPLESHEET)

sra_list = []
for _, row in df.iterrows():
    organism = row["Organism"].replace(" ", "_")
    sras = row["Values"].split(";")
    for sra in sras:
        sra_list.append((organism, sra))

# Do far there is a limmiter, since the memory is limited
if MAX_SRA is not None:
    sra_list = sra_list[:MAX_SRA]

# All final outputs
final_fastqs = []
for org, sra in sra_list:
    final_fastqs.append(f"data/SRA/{org}/{sra}_1.fastq")
    final_fastqs.append(f"data/SRA/{org}/{sra}_2.fastq")

final_contigs = [f"data/contigs/{org}/{sra}.contigs.fa" for org, sra in sra_list]

reference_genomes = {}
for folder in os.listdir('data/segments/'):
    folder_path = os.path.join('data/segments', folder)
    if not os.path.isdir(folder_path):
        continue
    reference_genomes[folder] = [
        os.path.join(folder_path, file)
        for file in os.listdir(folder_path)
        if file.endswith('.fasta')
    ]

contig_files = {
    folder: [
        os.path.join('data/contigs', folder, f)
        for f in os.listdir(os.path.join('data/contigs', folder))
        if f.endswith('.contigs.fa')
    ]
    for folder in os.listdir('data/contigs')
    if os.path.isdir(os.path.join('data/contigs', folder))
}

all_bams = []
for org in reference_genomes:
    for ref in reference_genomes[org]:
        ref_base = os.path.basename(ref).replace('.fasta','')
        for contig in contig_files.get(org, []):
            contig_base = os.path.basename(contig).replace('.contigs.fa','')
            all_bams.append(f"data/bam/{org}/{ref_base}__{contig_base}.sorted.bam")

# all_consensus = [
#     bam.replace("data/bam/", "data/consensus/") \
#        .replace(".sorted.bam", ".consensus.fa")
#     for bam in all_bams
# ]

all_consensus = [
    f"data/alignments/{org}/{os.path.basename(bam).replace('.sorted.bam','')}.consensus.fa"
    for bam in all_bams
    for org in [bam.split("/")[2]]   # extract organism folder from BAM path
]

rule all:
    input:
        final_contigs + all_bams

# Skipping SRA download for now, since we have contigs already.

# rule all:
#     input:
#         final_fastqs + final_contigs + all_bams

# rule download_sra:
#     output:
#         "data/SRA/{organism}/{sra}_1.fastq",
#         "data/SRA/{organism}/{sra}_2.fastq"
#     conda:
#         "envs/sra.yaml"
#     run:
#         import os
#         out_dir = os.path.dirname(output[0])
#         os.makedirs(out_dir, exist_ok=True)
#         shell(f"fasterq-dump {wildcards.sra} -O {out_dir} --threads 4")

rule download_logan_contigs:
    output:
        "data/contigs/{organism}/{sra}.contigs.fa"
    conda:
        "envs/aws.yaml"
    run:
        import os
        out_dir = os.path.dirname(output[0])
        os.makedirs(out_dir, exist_ok=True)
        zst_path = os.path.join(out_dir, f"{wildcards.sra}.contigs.fa.zst")
        # download
        shell(f"aws s3 cp s3://logan-pub/c/{wildcards.sra}/{wildcards.sra}.contigs.fa.zst {zst_path} --no-sign-request")
        # decompress
        shell(f"zstd -d {zst_path} -o {output[0]}")

rule align_contigs_to_refs:
    input:
        reference=lambda wildcards: reference_genomes[wildcards.organism],
        contig=lambda wildcards: contig_files[wildcards.organism]
    output:
        bam="data/bam/{organism}/{ref_basename}__{contig_basename}.sorted.bam",
        bai="data/bam/{organism}/{ref_basename}__{contig_basename}.sorted.bam.bai"
    conda:
        "envs/alignment.yaml"
    params:
        out_dir=lambda wildcards: f"data/bam/{wildcards.organism}"
    run:
        import os
        os.makedirs(params.out_dir, exist_ok=True)

        for ref in input.reference:
            ref_base = os.path.basename(ref).replace('.fasta','')
            for contig in input.contig:
                contig_base = os.path.basename(contig).replace('.contigs.fa','')
                out_bam = f"{params.out_dir}/{ref_base}__{contig_base}.sorted.bam"

                # concatenate reference (if multiple)
                tmp_ref = f"{out_bam}.ref.fa"
                shell(f"cat {ref} > {tmp_ref}")
                shell(f"samtools faidx {tmp_ref}")
                shell(f"bwa index {tmp_ref}")
                shell(f"minimap2 -a {tmp_ref} {contig} > {out_bam}.sam")
                shell(f"samtools view -bS {out_bam}.sam | samtools sort -o {out_bam}")
                shell(f"samtools index {out_bam}")
                shell(f"rm {out_bam}.sam {tmp_ref}*")

# rule consensus_from_bam:
#     input:
#         bam="data/bam/{organism}/{basename}.sorted.bam"
#     output:
#         consensus="data/alignments/{organism}/{basename}.consensus.fa",
#         bam_copy="data/alignments/{organism}/{basename}.sorted.bam"
#     conda:
#         "envs/alignment.yaml"
#     run:
#         import os
#         out_dir = os.path.dirname(output.consensus)
#         os.makedirs(out_dir, exist_ok=True)

#         # Copy BAM to alignments folder (optional)
#         if not os.path.exists(output.bam_copy):
#             shell(f"cp {input.bam} {output.bam_copy}")
        
#         # Generate consensus
#         shell(
#             f"samtools mpileup -uf {input.bam} | "
#             f"bcftools call -c | "
#             f"bcftools consensus -o {output.consensus}"
#         )

# rule consensus_from_bam:
#     input:
#         bam="data/bam/{organism}/{basename}.sorted.bam"
#     output:
#         consensus="data/alignments/{organism}/{basename}.consensus.fa",
#         bam_copy="data/alignments/{organism}/{basename}.sorted.bam"
#     conda:
#         "envs/alignment.yaml"
#     run:
#         import os
#         import glob
        
#         out_dir = os.path.dirname(output.consensus)
#         os.makedirs(out_dir, exist_ok=True)

#         # Copy BAM to alignments folder (optional)
#         if not os.path.exists(output.bam_copy):
#             shell(f"cp {input.bam} {output.bam_copy}")

#         # Make sure BAM index exists
#         shell(f"samtools index {input.bam} || true")

#         # Figure out reference FASTA dynamically
#         # Assume it's the one used for this BAM:
#         # match basename in data/segments/<organism>/
#         ref_dir = os.path.join("data/segments", wildcards.organism)
#         ref_files = glob.glob(os.path.join(ref_dir, "*.fasta"))
#         # pick the one that matches the prefix in BAM
#         ref = None
#         bam_base = wildcards.basename
#         for f in ref_files:
#             if os.path.basename(f).replace(".fasta","") in bam_base:
#                 ref = f
#                 break
#         if ref is None:
#             raise ValueError(f"No reference FASTA found for BAM {input.bam}")

#         # Generate consensus
#         shell(
#             f"samtools mpileup -uf {ref} {input.bam} | "
#             f"bcftools call -c | "
#             f"bcftools consensus -o {output.consensus}"
#         )

# rule consensus_from_bam:
#     input:
#         bam="data/bam/{organism}/{ref_basename}__{contig_basename}.sorted.bam",
#         ref="data/segments/{organism}/{ref_basename}.fasta"
#     output:
#         consensus="data/alignments/{organism}/{ref_basename}__{contig_basename}.consensus.fa",
#         bam_copy="data/alignments/{organism}/{ref_basename}__{contig_basename}.sorted.bam"
#     conda:
#         "envs/alignment.yaml"
#     run:
#         import os
#         os.makedirs(os.path.dirname(output.consensus), exist_ok=True)
#         if not os.path.exists(output.bam_copy):
#             shell(f"cp {input.bam} {output.bam_copy}")
#         # shell(f"samtools mpileup -uf {input.ref} {input.bam} | bcftools call -c | bcftools consensus -o {output.consensus}")
#         # inside align_contigs_to_refs
#         shell(f"samtools mpileup -uf {tmp_ref} {out_bam} | bcftools call -c | bcftools consensus -o data/alignments/{org}/{ref_base}__{contig_base}.consensus.fa")
rule consensus_from_bam:
    input:
        bam="data/bam/{organism}/{ref_basename}__{contig_basename}.sorted.bam",
        ref="data/segments/{organism}/{ref_basename}.fasta"
    output:
        consensus="data/alignments/{organism}/{ref_basename}__{contig_basename}.consensus.fa",
        bam_copy="data/alignments/{organism}/{ref_basename}__{contig_basename}.sorted.bam"
    conda:
        "envs/alignment.yaml"
    run:
        import os
        # make output dir
        os.makedirs(os.path.dirname(output.consensus), exist_ok=True)
        
        # copy BAM to alignments folder
        if not os.path.exists(output.bam_copy):
            shell(f"cp {input.bam} {output.bam_copy}")
        
        # generate consensus using the BAM and its correct reference
        shell(
            f"samtools mpileup -uf {input.ref} {input.bam} | "
            f"bcftools call -c | "
            f"bcftools consensus -o {output.consensus}"
        )





