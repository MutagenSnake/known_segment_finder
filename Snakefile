import pandas as pd

df = pd.read_csv("sample_sheet.csv")
df["Organism"] = df["Organism"].str.replace(" ", "_")

# Some SRAs failed to download from Logan's bucket; exclude them
with open("failing_contigs.txt") as f:
    failing_sras = [line.strip() for line in f]
df = df[~df["SRA"].isin(failing_sras)]



# limit SRA
MAX_SRA = 700
df = df.head(MAX_SRA)

rule all:
    input:
        # expand("data/contigs/{organism}/{sra}.contigs.fa",
        expand("data/contigs/{sra}.contigs.fa",
               zip, organism=df["Organism"], sra=df["SRA"]),
        expand("data/bam/{organism}/{ref_base}__{sra}.sorted.bam",
               zip, organism=df["Organism"], ref_base=df["Reference"].apply(lambda x: os.path.basename(x).replace(".fasta","")), sra=df["SRA"]),
        expand("data/alignments/{organism}/{ref_base}__{sra}.consensus.fa",
               zip, organism=df["Organism"], ref_base=df["Reference"].apply(lambda x: os.path.basename(x).replace(".fasta","")), sra=df["SRA"])



rule download_logan_contigs:
    output:
        # "data/contigs/{organism}/{sra}.contigs.fa"
        "data/contigs/{sra}.contigs.fa"
    conda:
        "envs/aws.yaml"
    run:
        import os
        out_dir = os.path.dirname(output[0])
        os.makedirs(out_dir, exist_ok=True)
        zst_path = os.path.join(out_dir, f"{wildcards.sra}.contigs.fa.zst")
        shell(f"aws s3 cp s3://logan-pub/c/{wildcards.sra}/{wildcards.sra}.contigs.fa.zst {zst_path} --no-sign-request")
        shell(f"zstd -d {zst_path} -o {output[0]}")


rule align_contigs_to_refs:
    input:
        ref="data/segments/{organism}/{ref_base}.fasta",
        # contig="data/contigs/{organism}/{sra}.contigs.fa"
        contig="data/contigs/{sra}.contigs.fa"
    output:
        bam="data/bam/{organism}/{ref_base}__{sra}.sorted.bam",
        bai="data/bam/{organism}/{ref_base}__{sra}.sorted.bam.bai"
    conda:
        "envs/alignment.yaml"
    run:
        import os
        os.makedirs(os.path.dirname(output.bam), exist_ok=True)
        shell(f"cat {input.ref} > {output.bam}.ref.fa")
        shell(f"samtools faidx {output.bam}.ref.fa")
        shell(f"bwa index {output.bam}.ref.fa")
        shell(f"minimap2 -a {output.bam}.ref.fa {input.contig} > {output.bam}.sam")
        shell(f"samtools view -bS {output.bam}.sam | samtools sort -o {output.bam}")
        shell(f"samtools index {output.bam}")
        shell(f"rm {output.bam}.sam {output.bam}.ref.fa*")


rule consensus_from_bam:
    input:
        bam="data/bam/{organism}/{ref_base}__{sra}.sorted.bam",
        ref="data/segments/{organism}/{ref_base}.fasta"
    output:
        "data/alignments/{organism}/{ref_base}__{sra}.consensus.fa"
    conda:
        "envs/alignment.yaml"
    run:
        import os
        out_dir = os.path.dirname(output[0])
        os.makedirs(out_dir, exist_ok=True)
        
        temp_bcf = os.path.join(out_dir, f"{wildcards.ref_base}__{wildcards.sra}.bcf")
        shell(f"bcftools mpileup -f {input.ref} {input.bam} | "
              f"bcftools call -c -Ob -o {temp_bcf}")
        shell(f"bcftools index {temp_bcf}")
        shell(f"bcftools consensus -f {input.ref} {temp_bcf} -o {output[0]}")
        shell(f"rm {temp_bcf} {temp_bcf}.csi")


