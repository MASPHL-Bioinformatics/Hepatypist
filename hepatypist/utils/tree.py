#!/usr/bin/env python
import base64
from io import BytesIO
import pandas as pd
import seaborn as sns
import baltic as bt
from ete3 import Tree
from Bio import Phylo

import matplotlib as mpl

mpl.use("agg")  # noninteractive backend, only allows writing to file

from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D


##############################################################################
##############################################################################

pd.set_option("display.max_columns", None)
# typeface='Helvetica Neue' ## set default matplotlib font and font size
FSIZE = 7
LEGFSIZE = 12

##############################################################################
##############################################################################

mpl.rcParams["font.weight"] = 500
# mpl.rcParams["font.weight"] = "heavy"
mpl.rcParams["axes.labelweight"] = 300
# mpl.rcParams['font.family']=typeface
mpl.rcParams["font.size"] = FSIZE
mpl.rcParams["lines.linewidth"] = 1


def prep_tree(treepath):
    """
    Reads in treefile as ete Tree object.
    Midpoint roots, adds genotype as node feature
    Returns the prepped Tree object
    """
    pt = Tree(treepath)  # quoted_node_names=True)#, quoted_node_names=True)
    pt.standardize()
    for leaf in pt.get_leaves():
        if "|" in leaf.name:
            leaf.add_features(genotype=leaf.name.split("|")[-1])
        else:
            leaf.add_features(genotype="QUERY")
    return pt


def count_nodes_to_root(tr, tipname):
    tip_node = next(iter(tr.getExternal(lambda x: x.name == tipname)), None)
    count = 0
    current_node = tip_node
    while current_node.parent:  # Traverse up to the root
        count += 1
        current_node = current_node.parent
    return count


def dist_to_root(tr, tipname):
    tip_node = next(iter(tr.getExternal(lambda x: x.name == tipname)), None)
    distance = 0.0
    current_node = tip_node
    while current_node.parent:  # Traverse up to the root
        distance += current_node.length  # Add branch length
        current_node = current_node.parent
    return distance


def tip_to_tip_dist(tr, tip1, tip2):
    # note: counts the tips themselves as external nodes; MIGHT also count the root...
    tip1node = next(iter(tr.getExternal(lambda x: x.name == tip1)), None)
    tip2node = next(iter(tr.getExternal(lambda x: x.name == tip2)), None)

    # Find the most recent common ancestor (MRCA)
    mrca = tr.commonAncestor([tip1node, tip2node])

    # If the MRCA is the root, avoid overestimation
    if mrca.parent is None:
        return dist_to_root(tr, tip1) + dist_to_root(tr, tip2)

    # Compute the distances from both tips to the MRCA
    distance1 = 0.0
    current_node = tip1node
    while current_node != mrca:
        distance1 += current_node.length
        current_node = current_node.parent

    distance2 = 0.0
    current_node = tip2node
    while current_node != mrca:
        distance2 += current_node.length
        current_node = current_node.parent

    # The distance between the two tips is the sum of both distances to the MRCA
    return distance1 + distance2


def tip_to_tip_nodes(tr, tip1, tip2):
    tip1node = next(iter(tr.getExternal(lambda x: x.name == tip1)), None)
    tip2node = next(iter(tr.getExternal(lambda x: x.name == tip2)), None)

    # Find the most recent common ancestor (MRCA)
    mrca = tr.commonAncestor([tip1node, tip2node])

    # If the MRCA is the root, avoid overestimation
    if mrca.parent is None:
        return count_nodes_to_root(tr, tip1) + count_nodes_to_root(tr, tip2)

    count1 = 0
    current_node = tip1node
    while current_node != mrca:  # Traverse up to the mrca
        count1 += 1
        current_node = current_node.parent

    count2 = 0
    current_node = tip2node
    while current_node != mrca:  # Traverse up to the mrca
        count2 += 1
        current_node = current_node.parent

    return count1 + count2


def get_query_dists_baltic(mtt, queryname):
    nonqueries = [
        leaf.name for leaf in mtt.getExternal() if leaf.traits["gt"] != "QUERY"
    ]
    distdict = {}
    for leaf in nonqueries:
        dist = tip_to_tip_dist(mtt, queryname, leaf)
        distdict[leaf] = dist
    return distdict


# def check_mono(queryname, querygt, allqueries, pt):
#     """
#     Collects a list of tree leaves with given genotype.
#     Checks whether all of them are monophyletic with ete3 check_monophyly function, returns boolean True if so.
#     Finds MRCA (most recent common ancestor) and checks whether all its descendants are the specified genotype, returns boolean True if so.
#     Also returns mrca node (to get support value)
#     """
#     # non-queries only
#     gt_leaves = [
#         leaf.name
#         for leaf in pt.get_leaves()
#         if (leaf.name == queryname or leaf.genotype == querygt)
#     ]
#     mrca = pt.get_common_ancestor(gt_leaves)
#     # mrca without queries...
#     if len(mrca) == len(gt_leaves):
#         mrca_descend = True
#     else:  # check if mrca contains other genotypes or just queries
#         mrca_descend = False

#     results = pt.check_monophyly(values=gt_leaves, target_attr="name")

#     if results[0]:
#         monophyletic = True
#     else:
#         monophyletic = False

#     return monophyletic, mrca_descend, mrca


def check_mono_baltic(queryname, allquerynames, gt, mtt):
    gt_leaf_names = [
        k.name
        for k in mtt.getExternal()
        if "|" in k.name and k.name.split("|")[-1] == gt
    ]
    ancestor = mtt.commonAncestor(
        mtt.getExternal(lambda k: k.name in gt_leaf_names or k.name == queryname)
    )  # [queryname, 'CY012432_NewZealand_2000.81643836', 'CY011960_NewZealand_2000.6630137', 'CY009404_NewZealand_2001.50410959', 'CY008139_NewZealand_2000.65479452'])) ## identify common ancestor node of two (or more) tips

    # if mrca is root....then what?
    # check descendant nodes:
    mrcadescendants = [w.name for w in ancestor.children if w.branchType == "leaf"]
    allpos = gt_leaf_names + allquerynames
    if "label" in ancestor.traits:
        mrcasup = float(ancestor.traits["label"])
    else:
        mrcasup = 0
    if all(item in allpos for item in mrcadescendants):
        return True, mrcasup
    else:
        return False, mrcasup


def load_prep_baltic_tree(rooted_treepath):
    mtt = bt.loadNewick(rooted_treepath, absoluteTime=False)
    # sorts children of each internal node in the tree according to specified sorting function and order (default is alphabetical - but can also sort by branch length, for example) - and redraws tree afterwards (required to update x and y positions of each branch) - useful after (for example) collapsing branches to make a less complex tree
    mtt.sortBranches()

    for k in mtt.getExternal():
        if "|" not in k.name:
            k.traits["gt"] = "QUERY"
        else:
            k.traits["gt"] = str(k.name.split("|")[-1])
    return mtt


def check_placement(queryname, blth, pt, supth, mtt):
    """
    Finds nearest neighbor to query and assigns putative genotype based on it.
    Checks monophyly of all leaves in this genotype including query.
    """
    if supth > 1.0:
        support_threshold = round(supth / 100, 2)
    else:
        support_threshold = supth

    # find query tip in tree
    querysearch = pt.search_nodes(name=queryname)

    if len(querysearch) == 0:
        return ValueError("Query ID not found in tree")

    querynode = querysearch[0]
    # finds nearest neighbor and gets its genotype
    dd = get_query_dists_baltic(mtt, querynode.name)
    nn_name = min(dd, key=dd.get)
    nn_dist = min(dd)
    nn = pt.search_nodes(name=nn_name)[0]

    numgen = len(
        pt.search_nodes(genotype=nn.genotype)
    )  # numgen non-query leaves in Tree (also excludes current query)
    query_parent = querynode.up  # gets bootstrap support for node up from query

    if len(blth) > 0:
        blast_gt = blth["match_gt"].iloc[0]
        blastth = blth["matchid"].iloc[0]
    else:
        blast_gt = "not_assignable_by_blast"

    # check_monophyly:   # returns (1) boolean if leaf names form monophyletic clade, (2) monophyletic or polyphyletic relationship, (3) set of leaf names breaking monophyly
    # monophyletic: Are all leaves in this genotype, including the query of interest but not any other queries, monophyletic?
    # mrca_descend: Are all descendants of the MRCA the same genotype? Includes this query but not other queries
    # mtt = bt.loadNewick(rooted_treepath, absoluteTime=False)
    # monophyletic, mrca_descend, mrca = check_mono(queryname, nn.genotype, pt)
    allquerynames = [k.name for k in mtt.getExternal() if "|" not in k.name]
    monophyletic, mrca_support = check_mono_baltic(
        queryname, allquerynames, nn.genotype, mtt
    )

    warnings = []
    if not monophyletic:
        warnings.append("NNGENOTYPE_NOT_MONOPHYLETIC")

    nodedict = {"querynode": querynode, "nn": nn}
    for nodename, node in nodedict.items():
        if node.support < support_threshold:
            warnings.append(f"{nodename}_SUPPORT_LOW")

    # checks whether NN genotype matches BLAST top hit genotype
    if nn.genotype != blast_gt:
        warnings.append(f"MISMATCH_GENOTYPE_NN_{nn.genotype}_BLAST_TOPHIT_{blast_gt}")

    if len(warnings) == 0:
        placement_gt = nn.genotype
    else:
        placement_gt = "UNASSIGNABLE"

    if "IC" in [nn.genotype, blast_gt]:
        placement_gt = "IC"
        warnings.append(
            "CANDIDATE_GENOTYPE_IC_PUTATIVE_IA_IC_NOT_FORMALLY_RECOGNIZED_BY_ICTV"
        )

    # should be added as a warning, but not failure-worthy
    if mrca_support < support_threshold:
        warnings.append(f"{nn.genotype}_mrca_SUPPORT_LOW")

    reportcols = [
        "queryname",
        "querynode_support",
        "queryparent_support",
        "nn_name",
        "nn_gt",
        "nn_dist",
        "nn_support",
        "numgen_nngt",
        "monophyletic",
        "mrca_support",
        "blast_tophit",
        "blast_gt",
        "placement_gt",
        "warnings",
    ]

    report = [
        queryname,
        querynode.support,
        query_parent.support,
        nn.name,
        nn.genotype,
        nn_dist,
        nn.support,
        numgen,
        monophyletic,
        mrca_support,
        blastth,
        blast_gt,
        placement_gt,
        "|".join(warnings),
    ]

    repdf = pd.DataFrame([report])
    repdf.columns = reportcols
    return repdf


def get_single_query_placement(
    rooted_treepath, queryname, blastparsed, output_tsv, supth
):
    """
    Loads a per-query reference tree (one query + full reference) and extracts
    placement statistics for that single query.  Called once per query from the
    place_single_query Snakemake rule.
    """
    if supth > 1.0:
        support_threshold = round(supth / 100, 2)
    else:
        support_threshold = supth

    pt = prep_tree(rooted_treepath)
    mtt = load_prep_baltic_tree(rooted_treepath)

    br = pd.read_csv(blastparsed, sep="\t")
    blth = br[br["queryid"] == queryname]

    repdf = check_placement(queryname, blth, pt, support_threshold, mtt)
    repdf.to_csv(output_tsv, sep="\t", index=False)
    return None


def root_tree(unrooted, rooted_treepath):
    """
    Uses biopython Phylo to midpoint root tree and output to new file
    """
    tr = Phylo.read(unrooted, "newick")
    tr.root_at_midpoint()
    Phylo.write(tr, rooted_treepath, "newick")
    return None


def get_placement(rooted_treepath, blastres, placement_hits, supth):
    if supth > 1.0:
        support_threshold = round(supth / 100, 2)
    else:
        support_threshold = supth

    pt = prep_tree(rooted_treepath)

    mtt = load_prep_baltic_tree(rooted_treepath)
    querylist = [
        leaf.name for leaf in mtt.getExternal() if leaf.traits["gt"] == "QUERY"
    ]

    br = pd.read_csv(blastres, sep="\t")

    dflist = []
    for q in querylist:
        blth = br[br["queryid"] == q]
        repdf = check_placement(q, blth, pt, support_threshold, mtt)
        if len(repdf) > 0:
            dflist.append(repdf)

    combreport = pd.concat(dflist)
    combreport.to_csv(placement_hits, sep="\t", index=False)
    return None


def get_branch_col(node, high_support_threshold, low_support_threshold):
    """
    Sets up criteria for branch coloration
    """
    if float(node.traits["label"]) <= 1.0:
        bs = 100 * float(node.traits["label"])
    else:
        bs = float(node.traits["label"])
    if bs >= high_support_threshold:
        return "darkgreen"
    elif high_support_threshold > bs > low_support_threshold:
        return "darkorange"
    else:
        return "crimson"


def make_treeviz(
    rooted_treepath, png_path, html_path, high_support_threshold, low_support_threshold
):
    btr = load_prep_baltic_tree(rooted_treepath)

    gtset = sorted(
        set([w.traits["gt"] for w in btr.getExternal() if w.traits["gt"] != "QUERY"])
    )

    if len(gtset) > 7:
        pal = sns.color_palette("Set2")  # , n_colors=len(gtset))
        pal.as_hex()

    else:
        # an ordered pull of sns Set2: aqua, coral, periwinkle, pink, spring green, yellow, tan, grey

        # mediumblue, goldenrod, royalgreen, red, lightpurple, browngray, lilacpink, graygreen, springgreen, cerulean

        # darker, more muted versions of above

        # Zack's suggestion, from https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40
        # minus hotpink because I don't think it stood out
        colorblind_plus_grey = [
            "#D25F89",
            "#14558E",
            "#99C791",
            "#43B19A",
            "#E5971E",
            "#4D0006",
            "#9C0AAB",
            "#b3b3b3",
        ]  # ,"#FF0776",]

        # To increase the saturation of a color, you can pass a value greater than 0 for the s parameter, while keeping the hue and lightness values the same.
        [sns.set_hls_values(ocol, h=None, l=None, s=1) for ocol in colorblind_plus_grey]

        pal = colorblind_plus_grey  # set2

    pad_nodes = {}

    cd = {}
    for n in range(0, len(gtset)):
        cd[gtset[n]] = pal[n]

    def c_func(k):
        return (
            cd[k.traits["gt"]]
            if k.is_leaf() and "QUERY" not in k.traits["gt"]
            else "black"
        )

    maxX = max(list([k.x for k in btr.Objects])) + 0.05

    figY = int(len(btr.getExternal()) / 8)
    fig, ax = plt.subplots(figsize=(10, figY), facecolor="w")
    plt.margins(x=0.01, y=0.01)  # reduce margins to 1% of the data range

    nodepadding = len(btr.getExternal()) * 0.4 / FSIZE
    regnodepadding = nodepadding * 0.5
    x_min, x_max = ax.get_xlim()
    ptoffset = (x_max - x_min) * 0.03

    querynodes = [n for n in btr.getExternal() if "QUERY" in n.traits["gt"]]

    # very minimal padding for ALL tips: increase spacing a bit so font label is more clearly visible
    for k in btr.getExternal():
        pad_nodes[k] = regnodepadding

    # pad nodes for the mrca for each genotype, or just the singular node if there's only one tip in the genotype
    for gt in gtset:
        gtnodes = [
            node
            for node in btr.getExternal()
            if "QUERY" not in node.traits["gt"] and node.traits["gt"] == gt
        ]
        if len(gtnodes) == 1:
            pad_nodes[gtnodes[0]] = nodepadding
        else:
            mrca = btr.commonAncestor(gtnodes)
            pad_nodes[mrca] = nodepadding

    # now pad queries...
    for qn in querynodes:
        pad_nodes[qn] = nodepadding

    btr.drawTree(pad_nodes=pad_nodes)

    # plot tree
    def x_attr(k):
        return k.x

    def pts_target(k):
        return k.is_leaf()  ## target list of which branches will be annotated

    def pts_x_attr(k):
        return float(k.x) + ptoffset

    def tree_c_func(k):
        return (
            get_branch_col(k, high_support_threshold, low_support_threshold)
            if "label" in k.traits
            else "darkgrey"
        )

    btr.plotTree(ax, x_attr=x_attr, colour=tree_c_func)
    btr.plotPoints(
        ax,
        target=pts_target,
        x_attr=pts_x_attr,
        size=20,
        outline=False,
        colour=c_func,
        zorder=100,
    )  ## tips

    for k in btr.getExternal():
        if (
            k.traits["gt"] != "QUERY"
        ):  ## genotype-colored lines from every tip's date to where we will draw text
            ax.plot(
                [k.x + ptoffset, maxX],
                [k.y, k.y],
                ls="--",
                lw=1,
                color=cd[k.traits["gt"]],
            )
            ax.text(
                maxX,
                k.y,
                k.name,
                color=cd[k.traits["gt"]],
                size=FSIZE,
                ha="left",
                va="center",
                zorder=102,
            )  # **kwargs)#style="oblique", ha="left", va="center")
        else:  # tag each query with a yellow flag box
            # Draw a gold dashed leader line to maxX (same pattern as reference
            # tips) so the label box is anchored at the right edge and cannot
            # overlap reference dashed lines or text.
            ax.plot(
                [k.x + ptoffset, maxX],
                [k.y, k.y],
                ls="--",
                lw=1,
                color="goldenrod",
                zorder=101,
            )
            ax.text(
                maxX,
                k.y,
                f"{k.name}  QUERY",
                fontsize=FSIZE - 1,
                ha="left",
                va="center",
                bbox=dict(
                    fc="#fffacd",  # light lemon-yellow, softer than pure yellow
                    edgecolor="goldenrod",
                    linewidth=0.6,
                    pad=3,
                    alpha=0.92,
                ),
                zorder=103,
            )

    # add legend for genotypes:
    square_verts = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
    custom_elements = []
    for k in sorted(cd.keys()):
        custom_elements.append(
            Polygon(
                square_verts, closed=True, facecolor=cd[k], edgecolor=None, label=f"{k}"
            )
        )
        # Line2D([0], [0], color=cd[k], marker="o", label=k))
    if len(querynodes) > 0:
        vertices = [(0, 0), (0, 1), (1, 1), (1, 0)]
        custom_elements.append(
            Polygon(
                vertices,
                closed=True,
                facecolor="#fffacd",
                edgecolor="goldenrod",
                linewidth=0.6,
                label="QUERY",
            )
        )
    legend1 = ax.legend(
        title="Genotype color key",
        title_fontsize=LEGFSIZE,
        fontsize=LEGFSIZE,
        handles=custom_elements,
        loc="lower left",
        bbox_to_anchor=(0.1, 0.1),
    )
    ax.add_artist(legend1)

    # add legend for branches ("darkgreen", "darkorange", "crimson")
    support_elements = [
        Line2D(
            [0], [0], color="crimson", marker=",", label=f"<{low_support_threshold}"
        ),
        Line2D(
            [0],
            [0],
            color="darkorange",
            marker=",",
            label=f"between {low_support_threshold} and {high_support_threshold}",
        ),
        Line2D(
            [0], [0], color="darkgreen", marker=",", label=f">={high_support_threshold}"
        ),
    ]
    legend2 = ax.legend(
        title="Branch support values",
        title_fontsize=LEGFSIZE,
        fontsize=LEGFSIZE,
        handles=support_elements,
        loc="lower left",
        bbox_to_anchor=(0.1, 0.3),
    )
    ax.add_artist(legend2)

    if len(querynodes) == 0:
        title = "Reference tree (no queries placed — all failed preprocessing or BLAST)"
    elif len(querynodes) == 1:
        title = "Placement of query (n = 1) in reference tree"
    else:
        title = f"Placement of queries (n = {len(querynodes)}) in reference tree"
    ax.set_title(
        title,
        fontsize=LEGFSIZE,
    )

    ax.tick_params(axis="x", size=0)  ## no labels
    ax.tick_params(axis="y", size=0)  ## no labels
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_yticklabels([])
    [
        ax.spines[loc].set_visible(False) for loc in ["top", "right", "left", "bottom"]
    ]  ## no axes

    tmpimg = BytesIO()
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    fig.savefig(tmpimg, dpi=200, bbox_inches="tight")
    plt.close()

    tmpimg.seek(0)
    plot_url = base64.b64encode(tmpimg.getvalue()).decode("utf8")

    htmls = f"""
        <!doctype html>
        <title>Query-ref tree</title>
        <section>
        <img src="data:image/png;base64, {plot_url}">
        </section>
        """

    with open(html_path, "w") as f:
        f.write(htmls)

    return None
