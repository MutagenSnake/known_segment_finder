import os
import pandas as pd
from Bio import SeqIO
import plotly.express as px

alignments_root = "data/alignments"
clusters_root = "data/clusters"
sra_file = "accessions_sra.csv"


df_sra = pd.read_csv(sra_file)
sra_dict = df_sra.set_index("Accession")[["Has_SRA","SRR_IDs"]].to_dict(orient="index")
rows = []

for organism in os.listdir(alignments_root):
    org_dir = os.path.join(alignments_root, organism)
    if not os.path.isdir(org_dir):
        continue

    # Load cluster info
    cluster_file = os.path.join(clusters_root, f"{organism}_clusters.txt")
    clusters = {}
    if os.path.isfile(cluster_file):
        with open(cluster_file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                cluster_name, seqs_str = line.split(":")
                cluster_name = cluster_name.strip()
                seq_ids = [s.strip() for s in seqs_str.split(",")]
                for seq_id in seq_ids:
                    clusters[seq_id] = cluster_name
    else:
        print(f"Warning: No cluster file for {organism}")

    # Iterate over consensus files
    for file in os.listdir(org_dir):
        if not file.endswith(".consensus.fa"):
            continue
        filepath = os.path.join(org_dir, file)

        # Extract header and SRA from filename
        base = os.path.basename(file).replace(".consensus.fa","")
        if "__" not in base:
            print(f"Warning: Unexpected filename format: {file}")
            continue
        header, sra = base.split("__", 1)

        # Get coverage, basically a percentage (not really using it at the moment)
        record = next(SeqIO.parse(filepath, "fasta"))
        seq = str(record.seq)
        total_len = len(seq)
        covered_len = sum(1 for b in seq if b.upper() != "N")
        coverage = (covered_len / total_len) * 100

        # Get cluster
        cluster = clusters.get(header, "Unknown")

        if cluster == "Unknown":
            print(f"Warning: Cluster not found for header '{header}' in organism '{organism}'")

        # Get SRA info if present
        sra_info = sra_dict.get(header, {"Has_SRA": False, "SRR_IDs": ""})

        rows.append({
            "Organism": organism,
            "Header": header,
            "SRA": sra,
            "Cluster": cluster,
            "Coverage": coverage,
            "Length": total_len,
            "Has_SRA": sra_info["Has_SRA"],
            "SRR_IDs": sra_info["SRR_IDs"]
        })

df = pd.DataFrame(rows)

df["Has_SRA_matched"] = df.apply(
    lambda row: str(row["SRA"]) in str(row["SRR_IDs"]).split(","),
    axis=1
)

srr_clusters = (
    df.dropna(subset=["SRR_IDs"])
      .groupby(["Organism", "Cluster"])["SRR_IDs"]
      .apply(list)
      .reset_index()
      .rename(columns={"SRR_IDs": "Cluster_SRR_IDs"})
)

df_merged = df.merge(srr_clusters, on=["Organism", "Cluster"], how="left")

# Ensure empty lists instead of NaN for Cluster_SRR_IDs
df_merged["Cluster_SRR_IDs"] = df_merged["Cluster_SRR_IDs"].apply(lambda x: x if isinstance(x, list) else [])

# Create Has_SRA_cluster_matched
df_merged["Has_SRA_cluster_matched"] = df_merged.apply(
    lambda row: row["SRA"] in row["Cluster_SRR_IDs"],
    axis=1
)

# Load all SRA metadata, important for dates and later maybe locations
df_all = pd.read_csv("all_sra_metadata.csv")



df_combined = df_merged.merge(
    df_all[['ID', 'releasedate', 'collection_date_sam']],
    left_on="SRA",   
    right_on="ID",   
    how="left"       # keep all rows from df_merged
).drop(columns=["ID"])

df_combined["best_date"] = df_combined["collection_date_sam"].fillna(
    df_combined["releasedate"]
)

# Convert to datetime, different formats present between dataframes
df_combined["best_date"] = pd.to_datetime(
    df_combined["best_date"], format="ISO8601", errors="coerce"
).dt.date

df_combined.to_csv("consensus_samples_with_dates.csv", index=False)

# Create inputs for IQ-Tree
alignment_root = "data/alignments"
output_root = "data/cluster_alignments"
# output_root = "output/cluster_alignments"

for (organism, cluster), group in df_combined.groupby(["Organism", "Cluster"]):
    out_dir = os.path.join(output_root, organism)
    os.makedirs(out_dir, exist_ok=True)

    out_file = os.path.join(out_dir, f"{cluster}.fasta")
    seen_headers = set()
    counter = 1

    with open(out_file, "w") as out_handle:
        for _, row in group.iterrows():
            header = row["Header"]
            sra = row["SRA"]
            date = row["best_date"]

            aln_path = os.path.join(alignment_root, organism, f"{header}__{sra}.consensus.fa")
            if not os.path.exists(aln_path):
                print(f"Warning: {aln_path} not found, skipping")
                continue

            for record in SeqIO.parse(aln_path, "fasta"):
                # base new header
                new_header = f"{header}__{sra}|{date}"

                # if duplicate, append a counter
                while new_header in seen_headers:
                    counter += 1
                    new_header = f"{header}__{sra}|{date}_{counter}"

                seen_headers.add(new_header)
                record.id = new_header
                record.description = ""
                SeqIO.write(record, out_handle, "fasta")

# Sunbursy HTML preparation
df_merged["ColorCategory"] = df_merged["Has_SRA_cluster_matched"].map({True: "Matched", False: "Not Matched"})

# Aggregate counts for the sunburst
agg = df_merged.groupby(
    ["Organism", "Cluster", "Header", "SRA", "ColorCategory"],
    dropna=False
).size().reset_index(name="Count")

# Build sunburst
fig = px.sunburst(
    agg,
    path=["Organism", "Cluster", "Header", "SRA"],
    values="Count",
    color="ColorCategory",
    color_discrete_map={"Matched": "red", "Not Matched": "green"},
    title="Consensus Samples Sunburst (SRA in Cluster SRR_IDs)"
)

fig.update_layout(margin=dict(t=40, l=0, r=0, b=0))

# Save to HTML
fig.write_html("output/samples_sunburst_cluster_matched.html")
