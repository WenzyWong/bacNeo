#################################################################################
# Extracting phylogenetic tree, OTU table, and tree-plot from Kraken2 report
## Inspired by a GitHub issue: https://github.com/jenniferlu717/KrakenTools/issues/46#issuecomment-2387744942
## A reference tool: https://github.com/jonasoh/gracken.git
## The output tree and OTU table can be imported into R environment for better illustration
## The generated picture displays the top species based on the overall abundance across samples,
### with top 10 species highlighted in bold font

# Author: Yunzhe WANG, yunzhewang24@m.fudan.edu.cn
# Updated: 2025-04-22
#################################################################################
# Other prerequisite: six, numpy, pandas, sip, PyQt5, opencv-python-headless
from ete3 import NCBITaxa, TreeStyle, NodeStyle, TextFace
import pandas as pd
from pandas.errors import EmptyDataError
import math, os, sys, glob, argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Creates a phylogenetic tree and OTU table from Kraken2 reports by pruning NCBI trees."
    )
    parser.add_argument(
        "--input_dir",
        "-i",
        required=True,
        help="Directory containing Kraken2 report files",
    )
    parser.add_argument(
        "--out_prefix",
        "-o",
        default="output",
        help="prefix for output files (default: output). Creates <prefix>.tree and <prefix>.otu.csv",
    )
    parser.add_argument(
        "--top_number",
        "-t",
        default=30,
        help="top ranked number based on species abundance across all samples",
    )
    args = parser.parse_args()
    return args

def read_kraken2_report(file_path):
    """Extract species abundances from Kraken2 report."""
    try:
        df = pd.read_csv(file_path, sep="\t", header=None)
    except EmptyDataError:
        # Return an empty DataFrame with the correct columns
        return pd.DataFrame(columns=["abundance", "taxid", "species"])
    if df.empty:
        return pd.DataFrame(columns=["abundance", "taxid", "species"])

    # Determine the format based on the number of columns
    num_cols = df.shape[1]

    # Different Kraken2 formats:
    # %coverage, #reads, #reads direct, rank code, NCBI ID, name
    if num_cols == 6:
        rank_col = 3
        taxid_col = 4
        name_col = 5
    # %coverage, #reads, #reads direct, ?, ?, rank code, NCBI ID, name
    elif num_cols == 8:
        rank_col = 5
        taxid_col = 6
        name_col = 7
    else:
        raise ValueError(f"Unsupported Kraken2 report format: {num_cols} columns")

    filtered = df[(df[rank_col] == "S") & (df[1] > 0)]

    result = pd.DataFrame()
    result["abundance"] = filtered[1]
    result["taxid"] = filtered[taxid_col]
    result["species"] = filtered[name_col].str.strip()

    return result

def get_sample_cols(cols):
    tax_cols = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    return [c for c in cols if c not in tax_cols]

def filter_bacteria_taxids(ncbi, taxids):
    """Filter TaxIDs to retain only those under Bacteria (TaxID: 2)."""
    bacteria_taxid = 2
    bacteria_taxids = []
    for taxid in taxids:
        try:
            lineage = ncbi.get_lineage(int(taxid))
            if bacteria_taxid in lineage:
                bacteria_taxids.append(taxid)
        except:
            continue  # Skip invalid TaxIDs
    return bacteria_taxids

def build_lineage_mapping(ncbi, otu_table):
    """
    Build a mapping from species to taxonomy info using NCBI lineage.
    """
    species_taxid = otu_table.groupby("species")["taxid"].first().to_dict()
    lineage_mapping = {}
    for sp, tid in species_taxid.items():
        lineage = ncbi.get_lineage(int(tid))
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)
        tax_cols = ["domain", "phylum", "class", "order", "family", "genus", "species"]
        taxonomy = {c: "" for c in tax_cols}
        for rank in taxonomy:
            if rank == "species":
                continue
            for t in lineage:
                actual_rank = ranks.get(t)
                if rank == "domain" and actual_rank == "superkingdom":
                    taxonomy[rank] = names.get(t, "")
                    break
                elif actual_rank == rank:
                    taxonomy[rank] = names.get(t, "")
                    break
        lineage_mapping[sp] = taxonomy
    return lineage_mapping

def save_tree(tree, output_file):
    tree.write(format = 1, outfile = output_file)
    tree.render(output_file.replace(".nwk", ".png"), w = 183, units = "mm")

def get_color_gradient(value, min_val, max_val):
    """Map a value to a color gradient inspired by RColorBrewer's RdYlBu (Red-Yellow-Blue)."""
    if max_val == min_val:
        ratio = 0
    else:
        ratio = (value - min_val) / (max_val - min_val)

    # Copy RdYlBu colors from R package-RcolorBrewer
    red = (215, 25, 28)    # #D7191C
    blue = (44, 123, 182)  # #2C7BB6

    # Interpolate directly from Blue to Red
    r = int(blue[0] + (red[0] - blue[0]) * ratio)
    g = int(blue[1] + (red[1] - blue[1]) * ratio)
    b = int(blue[2] + (red[2] - blue[2]) * ratio)

    # Convert to hex color code
    return f"#{r:02x}{g:02x}{b:02x}"

def save_tree_image(tree, output_file, species_abundance, top_species, height):
    """Render and save the phylogenetic tree with enhanced styling for abundance."""
    abundances = list(species_abundance.values())
    max_abundance = max(abundances) if abundances else 1
    min_abundance = min(abundances) if abundances else 0

    # Style nodes based on abundance
    def abundance_layout(node):
        nstyle = NodeStyle()
        if node.is_leaf():
            species = node.name
            abundance = species_abundance.get(species, 0)
            # Logarithmic thickness scaling for better contrast (1 to 15 pixels)
            if abundance > 0 and max_abundance > 0:
                log_abundance = math.log(abundance)
                log_max = math.log(max_abundance)
                log_min = math.log(min_abundance) if min_abundance > 0 else 0
                thickness = 1 + 14 * ((log_abundance - log_min) / (log_max - log_min)) if log_max > log_min else 1
            else:
                thickness = 1
            thickness = int(max(1, min(15, thickness)))
            nstyle["hz_line_width"] = thickness
            nstyle["vt_line_width"] = thickness
            # Apply colours
            color = get_color_gradient(abundance, min_abundance, max_abundance)
            nstyle["hz_line_color"] = color
            nstyle["vt_line_color"] = color
            # Add labels
            is_top_species = species in top_species
            label = TextFace(species, fsize=8, bold=is_top_species)
            node.add_face(label, column=0, position="branch-right")
        else:
            nstyle["hz_line_width"] = 1
            nstyle["vt_line_width"] = 1
            nstyle["hz_line_color"] = "#000000"
            nstyle["vt_line_color"] = "#000000"
        node.set_style(nstyle)

    # Apply the layout to TreeStyle
    ts = TreeStyle()
    ts.layout_fn = abundance_layout
    ts.show_leaf_name = False  # Disable default leaf name rendering
    ts.show_scale = False # Disable default scale

    ts.title.add_face(TextFace("Top Bacterial Species by Abundance", fsize=12, bold=True), column=0)

    ts.legend.add_face(TextFace("Colour(abundance): ", fgcolor="black"), column=0)
    ts.legend.add_face(TextFace("Blue (Low)", fgcolor="#2C7BB6"), column=1)
    ts.legend.add_face(TextFace(" -> ", fgcolor="black"), column=2)
    ts.legend.add_face(TextFace("Red (High)", fgcolor="#D7191C"), column=3)

    ts.legend.add_face(TextFace("Thickness(abundance): ", fgcolor="black"), column=0)
    ts.legend.add_face(TextFace(" Thin (Low)"), column=1)
    ts.legend.add_face(TextFace(" -> "), column=2)
    ts.legend.add_face(TextFace("Thick (High)"), column=3)
    ts.legend_position = 2

    ts.margin_right = 30

    tree.render(output_file, w=800, h=height, units="px", tree_style=ts)

# Generate OTU tables for de-contamination processes
def sanitize_name(name: str) -> str:
    """Create a filesystem- and R-friendly OTU ID from a species string."""
    if pd.isna(name):
        name = ""
    s = re.sub(r'[^A-Za-z0-9_.-]', '_', str(name))
    s = re.sub(r'_+', '_', s).strip('_')
    if s == "":
        s = "OTU"
    return s

def make_unique(ids: List[str]) -> List[str]:
    seen = {}
    out = []
    for id_ in ids:
        base = id_
        i = 1
        while id_ in seen:
            id_ = f"{base}_{i}"
            i += 1
        seen[id_] = True
        out.append(id_)
    return out

def write_decontam_outputs(out_prefix: str, wide_otu_table: pd.DataFrame, tax_cols: List[str]):
    """Write counts (OTU x samples) and taxonomy (OTU x tax_cols) TSVs suitable for decontam/phyloseq."""
    if wide_otu_table.empty:
        # Write empty files with headers so downstream code won't fail
        counts_path = f"{out_prefix}.counts.tsv"
        taxonomy_path = f"{out_prefix}.taxonomy.tsv"
        # counts: just header with OTU column (no samples)
        with open(counts_path, "w") as fh:
            fh.write("OTU_ID\n")
        # taxonomy: headers
        pd.DataFrame(columns=["OTU_ID"] + tax_cols).to_csv(taxonomy_path, sep="\t", index=False)
        return

    sample_cols = get_sample_cols(wide_otu_table.columns)
    counts = wide_otu_table.set_index("species")[sample_cols].fillna(0)
    try:
        counts = counts.astype(int)
    except Exception:
        counts = counts.fillna(0)
    # Sanitize OTU IDs
    original_species = list(counts.index.astype(str))
    otu_ids = [sanitize_name(s) for s in original_species]
    otu_ids = make_unique(otu_ids)
    counts.index = otu_ids

    counts_path = f"{out_prefix}.counts.tsv"
    counts.to_csv(counts_path, sep="\t", index=True, index_label="OTU_ID")

    # Taxonomy table
    tax_df = wide_otu_table.set_index("species")[tax_cols].copy()
    tax_df.index.name = "species"
    tax_df = tax_df.reindex(original_species)
    tax_df["OTU_ID"] = otu_ids
    taxonomy_path = f"{out_prefix}.taxonomy.tsv"
    tax_df.reset_index(drop=False).set_index("OTU_ID").to_csv(taxonomy_path, sep="\t", index=True)

def main():
    args = parse_args()

    otu_table = pd.DataFrame()
    file_pattern = "*.standard"
    matching_files = glob.glob(os.path.join(args.input_dir, file_pattern))
    if not matching_files:
        print(f"Error: No files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    for f in matching_files:
        name = os.path.splitext(os.path.basename(f))[0]
        sample_otu = read_kraken2_report(f)

        if sample_otu.empty or not all(
            col in sample_otu.columns
            for col in ["taxid", "abundance"]
        ):
            print(
                f"Warning: Skipping {f} because it's empty or malformed.",
                file=sys.stderr,
            )
            continue

        sample_otu["sample"] = name
        otu_table = pd.concat([otu_table, sample_otu], ignore_index=True)
    if otu_table.empty:
        print("Warning: No taxa detected in any samples. Writing empty outputs and exiting.", file=sys.stderr)
        write_decontam_outputs(args.out_prefix, pd.DataFrame(), ["domain", "phylum", "class", "order", "family", "genus", "species"])
        pd.DataFrame(columns=["domain","phylum","class","order","family","genus","species"]).to_csv(f"{args.out_prefix}.otu.csv", index=False)
        sys.exit(0)

    ncbi = NCBITaxa()
    taxid_list = list(otu_table["taxid"].unique())
    bacteria_taxids = filter_bacteria_taxids(ncbi, taxid_list)
    if not bacteria_taxids:
        print("Error: No Bacteria TaxIDs found.", file=sys.stderr)
        write_decontam_outputs(args.out_prefix, pd.DataFrame(), ["domain", "phylum", "class", "order", "family", "genus", "species"])
        sys.exit(1)

    otu_table = otu_table[otu_table["taxid"].isin(bacteria_taxids)]

    translator = ncbi.get_taxid_translator(bacteria_taxids)
    try:
        translator = ncbi.get_taxid_translator(bacteria_taxids)
    except Exception:
        translator = {}

    otu_table["species"] = otu_table["taxid"].map(lambda tid: translator.get(int(tid), str(tid)))

    wide_otu_table = otu_table.pivot_table(index="species", columns="sample", values="abundance", aggfunc="sum").reset_index().reset_index()

    mapping = build_lineage_mapping(ncbi, otu_table)

    tax_cols = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    for col in tax_cols[:-1]:
        wide_otu_table[col] = wide_otu_table["species"].map(lambda sp: mapping.get(sp, {}).get(col, ""))

    sample_cols = get_sample_cols(wide_otu_table.columns)
    wide_otu_table = wide_otu_table[tax_cols + sample_cols]

    wide_otu_table.to_csv(f"{args.out_prefix}.otu.csv", sep = ",", index = False)
    write_decontam_outputs(args.out_prefix, wide_otu_table, tax_cols)

    # Tree building: guard against empty taxid list and exceptions
    taxid_list_nonempty = [int(t) for t in otu_table["taxid"].unique() if str(t).strip() != ""]
    if not taxid_list_nonempty:
        print("Warning: No valid taxids to build trees.", file=sys.stderr)
        return

    try:
        tree = ncbi.get_topology(taxid_list_nonempty)
        for leaf in tree.get_leaves():
            try:
                leaf.name = translator.get(int(leaf.name), leaf.name)
            except Exception:
                pass
        tree.write(outfile=f"{args.out_prefix}.tree", format = 1, quoted_node_names = False)
        save_tree(tree, f"{args.out_prefix}_tree.nwk")
    except Exception as e:
        print(f"Warning: Failed to build full tree: {e}", file=sys.stderr)

    # Calculate total abundance per species across all samples
    species_abundance = otu_table.groupby("species")["abundance"].sum().to_dict()
    # Top species
    abundance_df = pd.DataFrame.from_dict(species_abundance, orient="index", columns=["abundance"])

    # Top 100
    top_100_species = abundance_df.sort_values(by="abundance", ascending=False).head(100).index.tolist()
    otu_table_100 = otu_table[otu_table["species"].isin(top_100_species)]
    top_100_taxids = otu_table_100["taxid"].unique().tolist()
    try:
        tree_100 = ncbi.get_topology([int(t) for t in top_100_taxids if str(t).strip() != ""])
        for leaf in tree_100.get_leaves():
            try:
                leaf.name = translator.get(int(leaf.name), leaf.name)
            except Exception:
                pass
        tree_100.write(outfile=f"{args.out_prefix}_100.tree", format = 1, quoted_node_names = False)
        save_tree(tree_100, f"{args.out_prefix}_tree_100.nwk")
    except Exception as e:
        print(f"Warning: Failed to build top-100 tree: {e}", file=sys.stderr)

    # Top 500
    top_500_species = abundance_df.sort_values(by="abundance", ascending=False).head(500).index.tolist()
    otu_table_500 = otu_table[otu_table["species"].isin(top_500_species)]
    top_500_taxids = otu_table_500["taxid"].unique().tolist()
    try:
        tree_500 = ncbi.get_topology([int(t) for t in top_500_taxids if str(t).strip() != ""])
        for leaf in tree_500.get_leaves():
            try:
                leaf.name = translator.get(int(leaf.name), leaf.name)
            except Exception:
                pass
        tree_500.write(outfile=f"{args.out_prefix}_500.tree", format = 1, quoted_node_names = False)
        save_tree(tree_500, f"{args.out_prefix}_tree_500.nwk")
    except Exception as e:
        print(f"Warning: Failed to build top-500 tree: {e}", file=sys.stderr)

    # Top N (user specified)
    top_n = int(args.top_number)
    top_n_species = abundance_df.sort_values(by="abundance", ascending=False).head(top_n).index.tolist()
    otu_table_n = otu_table[otu_table["species"].isin(top_n_species)]
    top_n_taxids = [t for t in otu_table_n["taxid"].unique().tolist() if str(t).strip() != ""]
    top_10_species = abundance_df.sort_values(by="abundance", ascending=False).head(10).index.tolist()
    try:
        tree_n = ncbi.get_topology([int(t) for t in top_n_taxids])
        for leaf in tree_n.get_leaves():
            try:
                leaf.name = translator.get(int(leaf.name), leaf.name)
            except Exception:
                pass
        save_tree_image(tree_n, f"{args.out_prefix}_tree_{top_n}.png", species_abundance, top_10_species, top_n * 30)
    except Exception as e:
        print(f"Warning: Failed to build top-{top_n} tree/image: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
