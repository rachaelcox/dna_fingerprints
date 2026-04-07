#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import numpy as np
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import ListedColormap, to_hex
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
from scipy.spatial.distance import squareform

# read in the files
fasta_path = "consensus/combined_consensus_sequences.fasta"
output_file = "figures/images/sequences_consensus_dendrogram.png"

# read in the fasta
records = list(SeqIO.parse(fasta_path, "fasta"))
df = pd.DataFrame({
    "id": [r.id for r in records],
    "len": [len(r.seq) for r in records],
    "seq": [str(r.seq).upper() for r in records],
})

# double check the file is reading correctly (628 sequences)
print("n sequences:", len(df))
print("\nLength summary:")
print(df["len"].describe())

# check there are only expected letters (nucleotides)
alphabet = set("".join(df["seq"].tolist()))
print("\nAlphabet observed:", "".join(sorted(alphabet)))

# define the logic to make the cosine distance matrix
def kmer_cosine_distance_matrix(seqs, k=5):
    vectorizer = CountVectorizer(analyzer="char", ngram_range=(k, k), lowercase=False)
    X = vectorizer.fit_transform(seqs)             
    Xn = normalize(X, norm="l2", axis=1) # normalize the kmer count vectors           

    S = (Xn @ Xn.T).toarray() # create array for how pairwise cosine similarity using the kmers
    D = 1.0 - S # calculate the distance between thes sequences

    D[D < 0] = 0 # distances are never less than 0
    D[D > 1] = 1 # distances are never greater than 1
    np.fill_diagonal(D, 0.0) # diagnol distance, self : self is always 0
    return D, X, vectorizer

# flavor meaning substrate and take everything before _k
def flavor_from_id(s: str) -> str:
    return s.split("_k", 1)[0]

# create the distance matrix
seqs = df["seq"].tolist()
D, X, vec = kmer_cosine_distance_matrix(seqs, k=5)

print("Distance matrix shape:", D.shape)
print("k-mer vocab size:", len(vec.get_feature_names_out()))

# Linkage
condensed = squareform(D, checks=False) # condense the distance matrix
Z = linkage(condensed, method="average") # hierarchial clustering using the average linkage

labels = df["id"].tolist()
flavors = [flavor_from_id(x) for x in labels]

# Color Mapping
uniq_flavors = list(dict.fromkeys(flavors))

# Using your specific palette
pal = [
    "#ff81c0",  # bubblegum pink
    "#c90076",  # raspberry pink
    "#b266ff",  # lavender purple
    "#7a1fa2",  # royal purple
    "#1f8f6b",  # deep jade green
    "#6ea8fe",  # periwinkle blue
    "#4c72b0",  # denim blue
    "#353540"   # charcoal
]

# Map unique flavors found in script to the palette colors
flavor2hex = {f: pal[i % len(pal)] for i, f in enumerate(uniq_flavors)}
flavor2rgb = {f: matplotlib.colors.to_rgb(flavor2hex[f]) for f in uniq_flavors}
bar_cmap = ListedColormap([flavor2hex[f] for f in uniq_flavors])
flavor2code = {f: i for i, f in enumerate(uniq_flavors)}

# Tree and Clade Logic
root = to_tree(Z, rd=False) # build the tree using the Z linkage
leaf_id_to_flavor = {i: flavors[i] for i in range(len(labels))}
node_id_to_flavor_set = {}

def collect_flavors(node):
    if node.is_leaf():
        fset = {leaf_id_to_flavor[node.id]}
    else:
        fset = collect_flavors(node.left) | collect_flavors(node.right)
    node_id_to_flavor_set[node.id] = fset
    return fset

collect_flavors(root)

GREY = "#404040"

def link_color_func(node_id: int):
    fset = node_id_to_flavor_set.get(node_id, set())
    if len(fset) == 1:
        return flavor2hex[next(iter(fset))]
    return GREY

# plot dendrogram
fig = plt.figure(figsize=(14, 6))
gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[12, 0.6], hspace=0.02)

ax_d = fig.add_subplot(gs[0, 0])
ax_s = fig.add_subplot(gs[1, 0], sharex=ax_d)

ddata = dendrogram(
    Z,
    no_labels=True,
    link_color_func=link_color_func,
    above_threshold_color=GREY,
    ax=ax_d,
)

ax_d.set_title("Hierarchical clustering of contiguous sequences")
ax_d.set_ylabel("Distance")

# Leaf strip
leaf_order = ddata["leaves"]
codes = np.array([flavor2code[flavors[i]] for i in leaf_order], dtype=int)
n = len(leaf_order)
xmin, xmax = 0, 10 * n

ax_s.imshow(
    codes.reshape(1, -1),
    aspect="auto",
    interpolation="nearest",
    cmap=bar_cmap,
    extent=[xmin, xmax, 0, 1],
)
ax_s.set_axis_off()

# Legend
handles = [mpatches.Patch(color=flavor2hex[f], label=f) for f in uniq_flavors]
ax_d.legend(handles=handles, title="Substrate", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)

# Save and Finish
plt.savefig(output_file, format="png", bbox_inches="tight")
print(f"\nDone! Output saved to {output_file}")