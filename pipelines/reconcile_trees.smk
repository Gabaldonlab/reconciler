configfile: "data/configs/reconcile_pdb.yml"

outdir=config["outdir"]
sptrees=config["sptrees"]
datadir=config["datadir"]

# SAMPLES WITH PHYID
# read phyids and get 000-code
with open(config['ids']) as input_ids:
    phyids = ["{:04d}".format(int(code)) for code in input_ids]

models = ["DL", "DTL"]

rule all:
    input:
        # expand(outdir+"phylome_{code}/generax_input", code=phyids),
        expand(outdir+"phylome_{code}/generax_{model}", code=phyids, model=models),
        expand(outdir+"phylome_{code}/results/phylome_{code}_{model}_rates.tsv", code=phyids, model=models),
        expand(outdir+"phylome_{code}/results/phylome_{code}_HGT.svg", code=phyids)
        # expand(outdir+"data/all_trees/rooted_mv_best_trees_{code}.nwk", code=phyids),


rule prepare_GeneRax:
    input:
        trees=datadir+"txts/best_trees_{code}.txt",
        aln=datadir+"aln/aln_{code}.tar.gz"
    output:
        temp(directory(outdir+"phylome_{code}/input"))
    shell:
        """
mkdir -p {output}/alns
tar tf {input.aln} | grep clean > {output}/{wildcards.code}_ids.txt
tar -C {output}/alns/ --strip-components 1 -xzf {input.aln} --files-from {output}/{wildcards.code}_ids.txt
python scripts/check_alns.py -i {output}/alns/
python scripts/prepare_trees.py -i {input.trees} -a {output}/alns/ -o {output} -m generax
"""

rule run_GeneRax:
    input:
        folder=outdir+"phylome_{code}/input",
        st=sptrees+"result_{code}_consensus.nwk"
    output:
        directory(outdir+"phylome_{code}/generax_{model}")
    params:
        model=lambda wcs: wcs.model
    benchmark:
        outdir+"benchmarks/{code}_{model}.txt"
    threads:
        12
    shell:
        """
mpiexec --oversubscribe -np {threads} generax --families {input.folder}/family.txt --species-tree {input.st} \
--rec-model Undated{params.model} --prefix {output} --per-family-rates --strategy RECONCILE
"""
# This corrects gene trees but its veeeeeeery slow
# --strategy SPR --max-spr-radius 1

# RECONCILE, undocumented option, does not correct gene trees but estimates DTL parameters,
# this may change in the future and EVAL should have this behaviout
# the phylogenetic likelihood is always 1!

# Statistical support for horizontal gene transfer
# https://groups.google.com/g/generaxusers/c/e_pmNI_aYBg

# to get statistical support either use --reconciliation-samples 100 to "bootstrap" reconciliations
# or compute aic difference between DTL (3*number of families) and DL (2*number of families)
# if the difference between the log likelihood is greater than one, then the AIC test says that there were transfers in this family.

rule explore_GeneRax:
    input:
        outdir+"phylome_{code}/generax_{model}"
    output:
        outdir+"phylome_{code}/results/phylome_{code}_{model}_rates.tsv"
    params:
        code = lambda wcs: wcs.code
    shell:
        """
Rscript scripts/explore_generax.R -i {input} -o {output} -c {params.code}
"""

rule mark_HGT:
    input:
        dl=outdir+"phylome_{code}/results/phylome_{code}_DL_rates.tsv",
        dtl=outdir+"phylome_{code}/results/phylome_{code}_DTL_rates.tsv"
    output:
        df=outdir+"phylome_{code}/results/phylome_{code}_HGT.tsv",
        # ids=outdir+"phylome_{code}/results/phylome_{code}_HGT.ids",
        files=outdir+"phylome_{code}/results/phylome_{code}_HGT.paths",
        plot=outdir+"phylome_{code}/results/phylome_{code}_HGT.svg"
    params:
        alpha=config['alpha'],
        t=config['thirdkind_collapse'],
        tk_config=config['thirdkind_config']
    shell:
        """
recdir=$(dirname {input.dl} | sed 's/results/generax_DTL\/reconciliations/')
Rscript scripts/mark_HGT.R -d {input.dl} -t {input.dtl} -o {output.df} -r $recdir -f {output.files} -s {params.alpha}
thirdkind -f {output.files} -m -t {params.t} -x -T 0 -o {output.plot} -c {params.tk_config}
"""
# -L to landscape mode
# convert -units PixelsPerInch  -density 300 test.svg test.png


# TODO
# --per-family-rates or --per-species-rates?

# clean or raw alignments??
# remove U in alignments!

# root gene trees: 3 different ways, min variance, midpoint and s2a dictionary
# outgroup rooting of gene trees???
# compare min variance and midpoint roting of gene trees
# root gene tree with QR or astralpro
# https://github.com/simonpenel/thirdkind/


# phyparts for reconciliation analysis and new trees!

# phylome 486 had a wrong tar alignment file