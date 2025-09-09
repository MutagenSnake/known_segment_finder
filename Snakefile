import pandas as pd

df = pd.read_csv("sample_sheet.csv")
df["Organism"] = df["Organism"].str.replace(" ", "_")

# Some SRAs failed to download from Logan's bucket; exclude them
with open("failing_contigs.txt") as f:
    failing_sras = [line.strip() for line in f]
df = df[~df["SRA"].isin(failing_sras)]

rule all:
    input:
        # All SRA contigs to be downloaded
        expand(
            "data/contigs/{sra}.contigs.fa",
            zip,
            organism=df["Organism"],
            sra=df["SRA"]
        ),
        # All consensus sequences to be generated
        expand(
            "data/alignments/{organism}/{ref_base}__{sra}.consensus.fa",
            zip,
            organism=df["Organism"],
            ref_base=df["Reference"].apply(lambda x: os.path.basename(x).replace(".fasta","")),
            sra=df["SRA"]
        ),
        # Run the finalizer code (So far it requires finalization.done to not be present. I think, not perfect.)
        "finalization.done"


rule download_logan_contigs:
    output:
        "data/contigs/{sra}.contigs.fa"
    conda:
        "envs/aws.yaml"
    run:
        # Simply download from Logan's public S3 bucket and decompress (Keeps the zst file, maybe unnecessary at this point.
        import os
        out_dir = os.path.dirname(output[0])
        os.makedirs(out_dir, exist_ok=True)
        zst_path = os.path.join(out_dir, f"{wildcards.sra}.contigs.fa.zst")
        shell(f"aws s3 cp s3://logan-pub/c/{wildcards.sra}/{wildcards.sra}.contigs.fa.zst {zst_path} --no-sign-request")
        shell(f"zstd -d {zst_path} -o {output[0]}")

rule consensus_from_contigs:
    input:
        ref="data/nucleotide_sequences/{organism}/{ref_base}.fasta",
        contig="data/contigs/{sra}.contigs.fa"
    output:
        "data/alignments/{organism}/{ref_base}__{sra}.consensus.fa"
    conda:
        "envs/alignment.yaml"
    threads: 4 # Can still use laptop for other things while this on is running (Fast anyway though)
    run:
        import os
        out_dir = os.path.dirname(output[0])
        os.makedirs(out_dir, exist_ok=True)

        bam   = os.path.join(out_dir, f"{wildcards.ref_base}__{wildcards.sra}.bam")
        vcf   = os.path.join(out_dir, f"{wildcards.ref_base}__{wildcards.sra}.vcf.gz")
        mask  = os.path.join(out_dir, f"{wildcards.ref_base}__{wildcards.sra}.nocov.bed")
        refix = os.path.join(out_dir, f"{wildcards.ref_base}__{wildcards.sra}.ref.fa")

        # Indexing, since other steps require indexed files
        shell("cp {input.ref} {refix}")
        shell("samtools faidx {refix}")
        shell("bwa index {refix}")

        # Generating BAM files (can be saved and inspected in IGV or Genious, but for now deleting to save space)
        shell("minimap2 -a {refix} {input.contig} | samtools sort -o {bam}")
        shell("samtools index {bam}")

        # Calling variants
        shell(
            "bcftools mpileup -f {input.ref} {bam} -Ou -a FORMAT/DP| "
            "bcftools call -mv -Oz -o {vcf}"
        )
        shell("bcftools index -f {vcf}")

        # If depth is 0, should be masked. Not my code though, but seems to work well. Even though molecular clock is not really working in the end.
        shell(
            "samtools depth -a {bam} | "
            "awk 'BEGIN{{OFS=\"\\t\"}} $3==0 {{ "
            "  if(c==$1 && p+1==$2) {{ p=$2 }} else {{ "
            "    if(c!=\"\") print c, s-1, p; "
            "    c=$1; s=$2; p=$2 "
            "  }} "
            "}} END{{ if(c!=\"\") print c, s-1, p }}' > {mask}"
        )

        # Generating final consensus alignment
        shell("bcftools consensus -f {input.ref} -m {mask} -H A {vcf} > {output[0]}")

        # Deleting intermediate files to save space
        shell("rm -f {bam} {bam}.bai {vcf} {vcf}.csi {mask} {refix}*")


rule finalize:
    output:
        touch("finalization.done")
    conda:
        "envs/summarizer.yaml"
    shell:
        """
        python scripts/summarizer.py
        # python scripts/create_iqtree_files.py # Not using this for now, takes too long, needs improvement before proper run.
        echo done > {output}
        """

