#!/usr/bin/env python3
"""Root the 1KP ASTRAL Newick tree with Chromista as outgroup.

[Generated via AI, tested by Cecilia]

The input has support and branch-length annotations on internal branches, but
no branch lengths on leaves.  ETE's topology-only output format (format 9)
drops the former, while its standard output format invents ``:1`` lengths for
the latter.  The small serializer below therefore writes the annotations held
by ETE only for internal branches.
"""

from pathlib import Path

from ete3 import Tree


# Input and output files
PROJECT = Path(__file__).resolve().parents[1]
TREE_DIR = PROJECT / "source_data/3.phylogenetic_tree/1kp_trees"
INPUT_TREE = TREE_DIR / "astral_trees_33_percent-FAA_estimated_species_tree.tree"
OUTPUT_TREE = TREE_DIR / "astral_trees_33_percent-FAA_estimated_species_tree.rooted.tree"


# The 33 taxa classified as "Outgroup" in the 1KP paper's
# supplementary table 1 listing the species (file 41586_2019_1693_MOESM3_ESM.xslx)
CHROMISTA = [
    "APTP", "ASZK", "BAJW", "BAKF", "BOGT", "DBYD", "EBWI", "FIDQ",
    "FIKG", "FOMH", "FSQE", "HFIK", "IAYV", "IRZA", "JCXF", "JGGD",
    "LDRY", "LIRF", "LLEN", "LXRN", "NMAK", "QLMZ", "RAPY", "RFAD",
    "ROZZ", "RWXW", "SRSQ", "ULXR", "VJED", "VKVG", "VRGZ", "VYER",
    "YRMA",
]


# 1. Read the unrooted Newick tree.
tree = Tree(INPUT_TREE.read_text().strip(), format=0)

# 2. Find the node containing all Chromista taxa.
chromista_node = tree.get_common_ancestor(*CHROMISTA)

# Check that this node contains Chromista only.
observed_chromista = set(chromista_node.get_leaf_names())
if observed_chromista != set(CHROMISTA):
    raise ValueError("The Chromista taxa do not form a single clade in this tree")

# 3. Place the root on the branch leading to the Chromista clade.
tree.set_outgroup(chromista_node)

# Check that Chromista is now one of the two groups descending from the root.
root_groups = [set(child.get_leaf_names()) for child in tree.children]
if len(root_groups) != 2 or set(CHROMISTA) not in root_groups:
    raise RuntimeError("The tree was not rooted correctly")


def annotated_newick(node):
    """Return an annotated Newick representation of an ETE tree node.

    ETE represents the annotation on the branch leading to a node using two
    attributes on that node: ``support`` is the internal-node label before the
    colon, and ``dist`` is the branch length after it.  Thus an internal node
    with two leaves is written as, for example::

        (UZWG,NHUA)100.0:1.50643867296

    Calling ``tree.write(format=9)`` cannot be used here because format 9 is
    topology-only and discards both values.  Conversely, ETE assigns its
    default distance of 1.0 to leaves whose input branches have no lengths, so
    its annotated writers would add ``:1`` to every leaf in this particular
    tree.  That would create data which were not present in the ASTRAL input.

    This recursive serializer matches the input's annotation scheme:

    * leaves contribute only their taxon names;
    * internal nodes contribute their parenthesized children followed by the
      support and distance stored by ETE; and
    * the root contributes only its children, because no branch exists above
      the root on which a support value or length could be defined.

    ``tree.set_outgroup`` has already transferred the branch annotations to
    the rerooted topology.  In particular, ETE splits the length of the branch
    selected for rooting equally between the two new root branches.  Writing
    the node attributes here therefore preserves ETE's rerooted branch data
    without introducing default leaf lengths.
    """
    if node.is_leaf():
        # The input ASTRAL tree has names, but no lengths, on terminal branches.
        return node.name

    descendants = ",".join(annotated_newick(child) for child in node.children)
    if node.is_root():
        # A root has no parent branch and consequently no branch annotation.
        return f"({descendants})"

    # For every non-root internal node, restore both parts of the annotation on
    # the branch connecting it to its parent.
    return f"({descendants}){node.support}:{node.dist}"


# 4. Preserve support values and branch lengths in the rooted Newick output.
OUTPUT_TREE.write_text(annotated_newick(tree) + ";\n")

print(f"Rooted tree written to: {OUTPUT_TREE}")
print(f"Root groups contain {len(root_groups[0])} and {len(root_groups[1])} taxa")
