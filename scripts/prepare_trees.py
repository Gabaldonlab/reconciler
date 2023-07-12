# Usage
import argparse
from ete3 import PhyloTree, Tree
import os

# Parse data
parser = argparse.ArgumentParser(
    description="Remove multifurcations from newick tree"
)

parser.add_argument(
    "-i",
    "--input",
    dest="tree",
    action="store",
    default="None",
    help="Input trees file",
    required=True
)
parser.add_argument(
    "-o",
    "--output",
    dest="out",
    action="store",
    default="None",
    help="output file",
    required=True
)
parser.add_argument(
    "-a",
    "--aln",
    dest="aln",
    action="store",
    default="None",
    help="alignment folder for generax",
    required=False
)
parser.add_argument(
    "-m",
    "--mode",
    choices=['asteroid', 'apro', 'duptree', 'stride', 'speciesrax', 'generax'],
    dest="mode",
    action="store",
    default="None",
    help="Prepare files for which software",
    required=True
)

args = parser.parse_args()


def obtain_duptree_file(treeFile, duptreeFile, midpoint=False, weighted=False):
    """Prepare files for duptree.
    Parameters
    ----------
    treeFile : str
        best trees file from PhylomeDB
    duptreeFile : str
        File where the tress will be written
    midpoint: bool
        Midpoint root gene trees?
    Returns
    -------
    type
        File that will be the input of duptree
    """

    outfile = open(duptreeFile, "w")
    for line in open(treeFile):
        line = line.strip()
        t = PhyloTree(line)
        # for leaf in t.iter_leaves():
        #     leaf.name = leaf.species
        string = ""
        if weighted:
            nms = [
                node.support
                for node in t.traverse()
                if not node.is_leaf() and node.support != 1.0
            ]
            if nms != []:
                weight = round(sum(nms) / (len(nms)), 3)
                if weight == 0:
                    weight = 0.01
                string = "[&U][&WEIGHT=" + str(weight) + "]"
            if nms == []:
                string = "[&U][&WEIGHT=0.01]"

        # this is very important! before each tree was rooted with midpoint rooting I don't know why.
        t.resolve_polytomy()
        # Root tree for new duptree
        # if root:
        if midpoint:
            t.set_outgroup(t.get_midpoint_outgroup())
            # if spe2age:
        outfile.write(string + t.write(format=9) + "\n")

    outfile.close()

def get_disco_data(gene_trees, out_file):
    """Obtain data to run astral-pro.
    Parameters
    ----------
    gene_trees : str
        Best trees file from PhylomeDB.
    out_file : str
        outputfile.
    Returns
    -------
    file
        Write outputfile in specified directory.
    """

    with open(out_file, "w") as o:
        with open(gene_trees) as t:
            for line in t:
                tree = Tree(line.strip())
                for leaf in tree.iter_leaves():
                    leaf.name = '_'.join(leaf.name.split("_")[::-1])
                string = tree.write()
                o.write(string + "\n")

def get_asteroid_data(gene_trees, out_file):
    """Obtain mapping file to run asteroid multicopy.
    Parameters
    ----------
    gene_trees : str
        Best trees file from PhylomeDB.
    out_file : str
        output mapping file.
    Returns
    -------
    file
        Write outputfile in specified directory.
    """

    with open(out_file, "w") as o:
        sp_dict = {}
        with open(gene_trees) as t:
            for line in t:
                tree = Tree(line.strip())
                for leaf in tree.iter_leaves():
                    sp = leaf.name.split("_")[1]
                    if sp not in sp_dict:
                        sp_dict[sp] = []
                    sp_dict[sp].append(leaf.name)
        for sp in sp_dict.keys():
            string = sp + ':' + ';'.join(sp_dict[sp])
            o.write(string + "\n")

def get_stride_data(gene_trees, out_dir):
    """Obtain data to run stride.
    Parameters
    ----------
    gene_trees : str
        Best trees file from PhylomeDB.
    out_file : str
        outputfile.
    Returns
    -------
    file
        Write outputfile in specified directory.
    """

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    with open(gene_trees, "r") as g:
        for line in g:
            line = line.strip().split()
            id = line[0]
            tree = Tree(line[3])
            out_file = out_dir + "/" + id + "_stride.nwk"

            for leaf in tree.iter_leaves():
                leaf.name = '_'.join(leaf.name.split("_")[::-1])

            tree.write(outfile=out_file)


def get_speciesrax_data(gene_trees, out_dir, keep_model=True):
    """Prepare data for speciesrax.
    Parameters
    ----------
    gene_trees : str
        Best trees file from phylomeDB.
    out_dir : str
        output directory.
    keep_model : bool
        Keep the best model found in phylomeDB
    Returns
    -------
    type
        A directory with the structure and files needed for SpeciesRax to run.
    """

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    if not os.path.isdir(out_dir + "/trees/"):
        os.mkdir(out_dir + "/trees/")
    if not os.path.isdir(out_dir + "/mapping/"):
        os.mkdir(out_dir + "/mapping/")

    out_map = out_dir + "/mapping/mapping"
    out_family = out_dir + "/family.txt"
    
    curpath = os.getcwd() + '/'

    with open(out_family, "w") as f:
        f.write("[FAMILIES]\n")
        with open(gene_trees, "r") as g:
            for line in g:

                line = line.strip().split()
                id = line[0]
                tree = Tree(line[3])

                leaves = [leaf.name for leaf in tree.iter_leaves()]
                species = list(set(leaf.split("_")[1] for leaf in leaves))

                sp_gene_dict = {}

                for s in species:
                    sp_genes = []
                    for gene in leaves:
                        nm = gene.split("_")
                        if nm[1] == s:
                            sp_genes.append(gene)
                    sp_gene_dict[s] = sp_genes

                out_mapfile = out_map + "_" + id + ".txt"
                with open(out_mapfile, "w") as m:
                    for k in sp_gene_dict.keys():
                        if len(sp_gene_dict[k]) > 0:
                            m.write(str(k) + ":" + ";".join(sp_gene_dict[k]) + "\n")

                f.write("- " + id + "\n")
                id_file = curpath + out_dir + "/trees/" + id + ".gene.newick"
                tree.write(outfile=id_file)
                f.write("starting_gene_tree = " + id_file + "\n")
                if keep_model:
                    model = line[1]
                    f.write("subst_model = " + model + "\n")
                f.write("mapping = " + curpath + out_mapfile + "\n")


def get_generax_data(gene_trees, out_dir, aln_dir, keep_model=True):
    """Prepare data for generax.
    Parameters
    ----------
    gene_trees : str
        Best trees file from phylomeDB.
    aln_dir: str
        folder with clean alignments
    out_dir : str
        output directory.
    keep_model : bool
        Keep the best model found in phylomeDB
    Returns
    -------
    type
        A directory with the structure and files needed for SpeciesRax to run.
    """

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    if not os.path.isdir(out_dir + "/trees/"):
        os.mkdir(out_dir + "/trees/")
    if not os.path.isdir(out_dir + "/mapping/"):
        os.mkdir(out_dir + "/mapping/")

    out_map = out_dir + "/mapping/mapping"
    out_family = out_dir + "/family.txt"

    curpath = os.getcwd() + '/'

    with open(out_family, "w") as f:
        f.write("[FAMILIES]\n")
        with open(gene_trees, "r") as g:
            for line in g:
                line = line.strip().split()
                id = line[0]

                aln_file = curpath + aln_dir + "/" + id + ".clean.fasta"

                if not os.path.exists(aln_file):
                    continue

                tree = Tree(line[3])

                leaves = [leaf.name for leaf in tree.iter_leaves()]
                species = list(set(leaf.split("_")[1] for leaf in leaves))

                sp_gene_dict = {}

                for s in species:
                    sp_genes = []
                    for gene in leaves:
                        nm = gene.split("_")
                        if nm[1] == s:
                            sp_genes.append(gene)
                    sp_gene_dict[s] = sp_genes

                out_mapfile = out_map + "_" + id + ".txt"
                with open(out_mapfile, "w") as m:
                    for k in sp_gene_dict.keys():
                        if len(sp_gene_dict[k]) > 0:
                            m.write(str(k) + ":" + ";".join(sp_gene_dict[k]) + "\n")

                f.write("- " + id + "\n")
                id_file = curpath + out_dir + "/trees/" + id + ".gene.newick"
                tree.write(outfile=id_file)
                f.write("starting_gene_tree = " + id_file + "\n")
                f.write("alignment = " + aln_file + "\n")
                f.write("mapping = " + curpath + out_mapfile + "\n")
                if keep_model:
                    model = line[1]
                    if "JTTDCMut" in model:
                        model = model.replace("JTTDCMut","JTT-DCMut")
                    f.write("subst_model = " + model + "\n")

if __name__=="__main__":
    if args.mode=="duptree":
        obtain_duptree_file(args.tree, args.out, weighted=True)
    if args.mode=="asteroid":
        get_asteroid_data(args.tree, args.out)
    if args.mode=="apro":
        get_disco_data(args.tree, args.out)
    if args.mode=="stride":
        get_stride_data(args.tree, args.out)
    if args.mode=="speciesrax":
        get_speciesrax_data(args.tree, args.out, keep_model=False)
    if args.mode=="generax":
        get_generax_data(args.tree, args.out, args.aln, keep_model=True)