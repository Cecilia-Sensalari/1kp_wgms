"""Shared helpers for the ortholog-Ks / tree-path-distance calibration.

Everything the calibration scripts need that is not specific to one of them: the two trees,
tip-to-tip path distances, the species metadata (1KP code, latin name, FASTA), the paranome
bin heights and the coarse taxonomic groups.

Paths are absolute on purpose: these scripts are run by hand, from anywhere, under
/software/shared/apps/x86_64/python/3.11.4/bin/python3 (stdlib + numpy/scipy/matplotlib).
"""
import os, re, csv, glob, json, statistics

ROOT     = "/group/esb/cesen/1kp"
DATA_DIR = os.path.join(ROOT, "ks_analysis/1kp/calibration")

CONCAT_NEXUS  = os.path.join(ROOT, "source_data/3.phylogenetic_tree/1kp_trees/SupFigs1-3-FigTree.tre")
ANALYSIS_TREE = os.path.join(ROOT, "source_data/3.phylogenetic_tree/1kp_trees/"
                                   "astral_trees_33_percent-FAA_estimated_species_tree.rooted."
                                   "wgm_suptab3_mrca.nhx.manual_fixes.tree")
TAXONOMY_CSV  = os.path.join(ROOT, "source_data/1.species_dataset/"
                                   "1kp_paper_2019_suptab1_species_ncbi_taxonomy.csv")
SPECIES_TSV   = os.path.join(ROOT, "source_data/1.species_dataset/1kp_paper_2019_suptab1_species.tsv")
PARALOG_DIR   = os.path.join(ROOT, "ks_analysis/1kp/paralog_runs/paralog_distributions")
DATASET_TSV   = os.path.join(ROOT, "ks_analysis/1kp/focal_branch_to_dataset.tsv")
CONFIG_GLOB   = os.path.join(ROOT, "ks_analysis/1kp/datasets/*/configs/config_*.txt")

BIN_HEIGHTS_TSV      = os.path.join(DATA_DIR, "paralog_ks_bin_heights.tsv")
CALIBRATION_PAIRS    = os.path.join(DATA_DIR, "calibration_pairs.tsv")
CONCAT_NEWICK        = os.path.join(DATA_DIR, "concat_tree.nwk")

# the 20 bins of the median-bin-height indicator
BIN_LO, BIN_HI, BIN_WIDTH = 0.5, 2.5, 0.1
N_BINS = int(round((BIN_HI - BIN_LO) / BIN_WIDTH))


# --------------------------------------------------------------------------- trees
_COMMENT = re.compile(r"\[[^\]]*\]")
_PUNCT = set(":,();")


class Node:
    __slots__ = ("children", "parent", "name", "blen", "_leaves")

    def __init__(self):
        self.children = []
        self.parent = None
        self.name = None
        self.blen = 0.0
        self._leaves = None

    def leaves(self):
        if self._leaves is None:
            self._leaves = [self.name] if not self.children else \
                           [x for c in self.children for x in c.leaves()]
        return self._leaves


def parse_newick(text):
    """Parse a newick string into a Node tree. Tolerates NHX comments and quotes."""
    t = _COMMENT.sub("", text).replace("'", "")
    root = Node(); cur = root; i = 0; n = len(t)

    def number(i):
        start = i
        while i < n and t[i] in "0123456789.eE+-":
            i += 1
        return (float(t[start:i]) if i > start else 0.0), i

    while i < n:
        ch = t[i]
        if ch == "(":
            child = Node(); child.parent = cur; cur.children.append(child); cur = child; i += 1
        elif ch == ")":
            i += 1
            while i < n and t[i] not in _PUNCT:      # internal label / support value
                i += 1
            cur = cur.parent
            if i < n and t[i] == ":":
                cur.blen, i = number(i + 1)
        elif ch == ",":
            sib = Node(); sib.parent = cur.parent; cur.parent.children.append(sib); cur = sib; i += 1
        elif ch == ";":
            i += 1
        else:
            start = i
            while i < n and t[i] not in _PUNCT:
                i += 1
            cur.name = t[start:i].strip()
            if i < n and t[i] == ":":
                cur.blen, i = number(i + 1)
    while len(root.children) == 1 and not root.name:
        root = root.children[0]
        root.parent = None
    return root


def extract_concat_newick(write_cache=True):
    """Pull the `tree concat` line out of the NEXUS file.

    The ASTRAL tree in the same file cannot be used for path distances: its internal branches
    are in coalescent units and every terminal branch is a placeholder, so this refuses to
    return a tree whose terminal lengths are all identical.
    """
    text = None
    with open(CONCAT_NEXUS) as fh:
        for line in fh:
            s = line.strip()
            if s.lower().startswith("tree concat"):
                text = s[s.index("("):]
                break
    if text is None:
        raise SystemExit(f"no 'tree concat' line in {CONCAT_NEXUS}")
    root = parse_newick(text)
    terminal = [n.blen for n in _walk(root) if not n.children]
    if len(set(terminal)) <= 1:
        raise SystemExit("all terminal branch lengths are identical - this is not the concat tree")
    if write_cache:
        os.makedirs(DATA_DIR, exist_ok=True)
        with open(CONCAT_NEWICK, "w") as fh:
            fh.write(text if text.endswith(";") else text + ";")
    return root


def _walk(node):
    yield node
    for c in node.children:
        yield from _walk(c)


def tip_index(root):
    return {n.name: n for n in _walk(root) if n.name and not n.children}


class PathDistances:
    """Tip-to-tip path distance on a tree with branch lengths, with memoised lookups."""

    def __init__(self, root):
        self.index = tip_index(root)
        self._to_root = {name: self._path(node) for name, node in self.index.items()}
        self._cache = {}

    @staticmethod
    def _path(node):
        out = []
        while node is not None:
            out.append(node)
            node = node.parent
        return out

    def __contains__(self, name):
        return name in self.index

    def __call__(self, a, b):
        if a == b:
            return 0.0
        key = (a, b) if a < b else (b, a)
        if key in self._cache:
            return self._cache[key]
        if a not in self._to_root or b not in self._to_root:
            self._cache[key] = None
            return None
        pa, pb = self._to_root[a], self._to_root[b]
        seen = {id(x): i for i, x in enumerate(pb)}
        result = None
        for i, node in enumerate(pa):
            if id(node) in seen:
                j = seen[id(node)]
                result = sum(x.blen for x in pa[:i]) + sum(x.blen for x in pb[:j])
                break
        self._cache[key] = result
        return result


def load_concat_tree():
    """Return (root, PathDistances) for the concatenation tree."""
    if os.path.exists(CONCAT_NEWICK):
        root = parse_newick(open(CONCAT_NEWICK).read().strip())
    else:
        root = extract_concat_newick()
    return root, PathDistances(root)


def load_analysis_tree():
    return parse_newick(open(ANALYSIS_TREE).read().strip())


def branch_nodes(root):
    """Every node that represents a branch, in the traversal order used by the dataset table."""
    out = []

    def walk(n):
        if n.parent is not None:
            out.append(n)
        for c in n.children:
            walk(c)

    walk(root)
    return out


# ------------------------------------------------------------------- species metadata
# Published genomes are absent from the 1KP supplementary table; the tree labels them with a
# 5-letter prefix plus an assembly version, and the per-dataset ksrates configs carry the
# short code, the latin name and the CDS path for every species, genomes included.
def species_metadata(cache=os.path.join(DATA_DIR, "species_metadata.json")):
    """tree tip label -> {"code", "latin", "fasta"}, harvested from the generated configs."""
    if os.path.exists(cache):
        return json.load(open(cache))
    latin, fasta = {}, {}
    for path in glob.glob(CONFIG_GLOB):
        text = open(path).read()
        m = re.search(r"latin_names\s*=\s*(.*?)\n\s*\n", text, re.S)
        f = re.search(r"fasta_filenames\s*=\s*(.*?)\n\s*\n", text, re.S)
        if not m or not f:
            continue
        for entry in m.group(1).split(","):
            if ":" in entry:
                code, value = entry.split(":", 1)
                latin.setdefault(code.strip(), value.strip())
        for entry in f.group(1).split(","):
            if ":" in entry:
                code, value = entry.split(":", 1)
                fasta.setdefault(code.strip(), value.strip())
    meta = {value.split()[-1]: {"code": code, "latin": value, "fasta": fasta.get(code)}
            for code, value in latin.items()}
    os.makedirs(DATA_DIR, exist_ok=True)
    json.dump(meta, open(cache, "w"), indent=0, sort_keys=True)
    return meta


def resolve_tip(name, meta):
    """Map a latin-name suffix to a tree tip label.

    The master_1_2026 database stores Vitis as "Vitis vinifera Vitvi" while the tree calls it
    "Vitvi_Genoscope.12X", so a bare 5-letter genome code has to be matched on the prefix.
    """
    if name in meta:
        return name
    for tip in meta:
        if tip.split("_")[0] == name:
            return tip
    return None


def lineages():
    """1KP sample id -> NCBI lineage string (genomes are not in this table)."""
    out = {}
    with open(TAXONOMY_CSV, newline="") as fh:
        for row in csv.DictReader(fh):
            out[row["sample_id"]] = row.get("lineage_names", "")
    return out


def busco():
    out = {}
    with open(SPECIES_TSV) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                out[row["1KP Index ID"]] = float(row["% BUSCOs (complete+fragmented)"])
            except (TypeError, ValueError, KeyError):
                pass
    return out


# Published genomes, which carry no 1KP lineage: assigned by hand, verified against NCBI for
# every species that does appear in the 1KP table (826 of 847 in the lycophyte-sister clade,
# zero disagreements).
GENOME_GROUP = {
    "Chlre_v5.5": "alga", "Chlva_v1.0": "alga", "Micpu_v3.0": "alga", "Ostlu_v2.0": "alga",
    "Cyame_v1.0": "alga", "Chabr_v0.1": "alga", "Klefl_v1.0": "alga",
    "Phypa_v3.0": "bryophyte", "Selmo_v1.0": "lycophyte", "Pinta_v2.0": "gymnosperm",
}
GROUP_ORDER = ["alga", "bryophyte", "lycophyte", "fern", "gymnosperm", "angiosperm"]
_GROUP_KEYS = (("Lycopodiopsida", "lycophyte"), ("Polypodiopsida", "fern"),
               ("Acrogymnospermae", "gymnosperm"), ("Magnoliopsida", "angiosperm"),
               ("Liliopsida", "angiosperm"), ("Spermatophyta", "gymnosperm"),
               ("Bryophyta", "bryophyte"), ("Marchantiophyta", "bryophyte"),
               ("Anthocerotophyta", "bryophyte"))


def groups(meta=None, lin=None):
    """tree tip label -> one of GROUP_ORDER."""
    meta = meta or species_metadata()
    lin = lin if lin is not None else lineages()
    out = {}
    for tip, info in meta.items():
        if tip in GENOME_GROUP:
            out[tip] = GENOME_GROUP[tip]
            continue
        lineage = lin.get(info["code"], "")
        out[tip] = next((g for key, g in _GROUP_KEYS if key in lineage),
                        "angiosperm" if tip not in GENOME_GROUP and not lineage else "alga")
    return out


# ------------------------------------------------------------------- paranome bin heights
def bin_heights(path=BIN_HEIGHTS_TSV):
    """tree tip label -> list of N_BINS weighted-node bar heights over Ks 0.5-2.5."""
    if not os.path.exists(path):
        raise SystemExit(f"{path} not found - run paralog_ks_bin_heights.py first")
    out = {}
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < N_BINS + 1:
                continue
            out[f[0]] = [float(x) for x in f[1:N_BINS + 1]]
    return out


def median_bin_heights(path=BIN_HEIGHTS_TSV):
    return {k: statistics.median(v) for k, v in bin_heights(path).items()}


# ------------------------------------------------------------------- subtree construction
# These mirror what make_branch_subtrees.py does, so that the figures describe the same
# clades the survey actually analyses.
A_PER_CHILD, A_MAX = 2, 4


def pick_A(node, median_bin):
    """The A representatives of a test branch: the best two species of each direct child."""
    children = node.children if node.children else [node]
    out = []
    for child in children:
        out.extend(sorted(child.leaves(), key=lambda s: -median_bin.get(s, 0))[:A_PER_CHILD])
    return out[:A_MAX]


def sister_clades(node, how_many=3):
    """The B, C and D clades of a test branch: the siblings met walking up to the root."""
    out, cur = [], node
    while cur.parent is not None and len(out) < how_many:
        out.extend(k for k in cur.parent.children if k is not cur)
        cur = cur.parent
    return out[:how_many]


def clade_containing(nodes, members, min_fraction=0.85, outsider_penalty=3):
    """The node that best corresponds to a taxonomic group: most members, fewest outsiders."""
    best, best_score = None, None
    for node in nodes:
        leaves = set(node.leaves())
        inside = len(leaves & members)
        if inside < min_fraction * len(members):
            continue
        score = inside - outsider_penalty * len(leaves - members)
        if best_score is None or score > best_score:
            best, best_score = node, score
    return best


# ------------------------------------------------------------------- path -> expected Ks
# Fitted on the 124 calibration pairs with a path distance of at most RELIABLE_PATH (see
# ortholog_ks_path_distance_calibration.md, 3.3). A power law, not a single conversion factor:
# the relationship is concave, so any constant factor fitted on the shallow pairs
# UNDERESTIMATES the divergence of the deep ones - the one direction the selection rule must
# not err in, since an underestimate lets a saturated candidate through the budget.
#
# The fit is restricted to the pairs whose ortholog Ks distribution actually HAS a peak. Beyond
# path ~1.0 the distributions flatten into a smear against the max_ks_orthologs cutoff and the
# "peak" the KDE returns is an artefact of that cutoff, not a measurement (3.5). The earlier
# version of this fit used all 205 pairs and so was fitted partly on those artefacts. Over the
# pairs that are real the two fits agree closely (median |error| 7.5% against 8.1%); the
# difference is in the extrapolation, which is exactly what we no longer quote.
POWER_A, POWER_EXPONENT = 0.4531, 1.266
# p90 of observed/predicted over the fitted pairs: applying it stops ~90% of pairs from being
# underestimated, at the cost of a systematic overestimate. For the saturation flag, not for
# the reported value.
CONSERVATIVE_MARGIN = 1.19
# SD of log(observed/predicted) over the fitted pairs: the calibration error is multiplicative,
# so the ranking noise is a fraction of the predicted value, not a fixed number of Ks units.
LOG_NOISE_SD = 0.141

# ------------------------------------------------------------------- reliability of a path
# Measured on the 86 ortholog Ks distributions available on disk, as the fraction of Ks values
# within +-0.5 of the distribution mode (a flat smear over the 0-10 window gives 10%):
#
#     path   0-0.6  0.6-0.8  0.8-1.0  1.0-1.2  1.2-1.5  1.5-1.8  1.8-2.1  2.1-2.4  2.4-3.6
#     conc     68%      55%      41%      21%      20%      17%      16%      16%      16%
#
# A real mode survives to about path 1.0 and is gone by 1.5. These two numbers therefore split
# the range into three, and they are thresholds on PATH rather than on Ks because the quantity
# that degrades is the measurability of the Ks peak, not its value.
RELIABLE_PATH, BORDERLINE_PATH = 1.0, 1.5
RELIABILITY_LEVELS = ("reliable", "borderline", "unreliable")


def path_reliability(path):
    """How far an ortholog Ks peak measured over this path distance can be trusted.

    'reliable'   - the ortholog Ks distribution has a clear mode; the peak is a measurement.
    'borderline' - the mode is broad and a large share of the Ks values fall outside the
                   window; usable for ordering candidates, to be reported with a warning.
    'unreliable' - the distribution is indistinguishable from a flat smear against the
                   max_ks_orthologs cutoff. The branch is still analysed (its paranome Ks
                   peak is unaffected), but the rate adjustment built on this pair cannot
                   position the divergence, so the BRANCH ASSIGNMENT is what is in doubt -
                   not the WGM call itself.
    """
    if path is None:
        return None
    if path <= RELIABLE_PATH:
        return "reliable"
    return "borderline" if path <= BORDERLINE_PATH else "unreliable"


def path_to_ks(path, conservative=False):
    """Expected ortholog Ks peak for a concat-tree path distance.

    Always returns a number, including past RELIABLE_PATH: the selection rule needs an
    ORDERING of the candidates, and path distance stays monotone in divergence even where the
    Ks peak itself has saturated. Use report_ks() for anything shown to a reader, which
    refuses to quote an extrapolated value as if it were a measurement.
    """
    ks = (path / POWER_A) ** POWER_EXPONENT
    return ks * CONSERVATIVE_MARGIN if conservative else ks


def report_ks(ks, path=None, verbose=False):
    """Format an expected Ks for reporting, refusing to quote a value off the calibration.

    Beyond RELIABLE_PATH there is no measured Ks peak to compare a prediction against, so the
    number is reported as a bound rather than as a value. The short form is the default
    because this string goes into figures and tables, where a sentence does not fit; pass
    verbose=True for the one that names the path.
    """
    if ks is None:
        return "NA"
    level = path_reliability(path) if path is not None else "reliable"
    if level == "reliable":
        return f"{ks:.2f}"
    limit = path_to_ks(RELIABLE_PATH)
    if verbose:
        return f">{limit:.1f} ({level}: path {path:.2f} is beyond the calibration)"
    return f">{limit:.1f} ({level})"


def expected_ks(candidate, a_species, dist, conservative=False):
    """Expected ortholog Ks of a candidate, from the median path distance to the A species.

    The median over A, not the mean: the A representatives are a small, deliberately
    heterogeneous sample and one unusually fast A species should not shift every candidate.
    """
    path = median_path(candidate, a_species, dist)
    return None if path is None else path_to_ks(path, conservative)


def median_path(candidate, a_species, dist):
    """Median concat-tree path distance from the A species to a candidate."""
    paths = [dist(a, candidate) for a in a_species]
    paths = [p for p in paths if p is not None]
    return statistics.median(paths) if paths else None
