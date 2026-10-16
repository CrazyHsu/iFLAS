"""GO enrichment plotting stage backed by the packaged Rscript bridge."""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterable

from ..rscript import run_r_function
from ..runtime import ToolRunner


PLOT_TARGET_GENES_GO_ENRICHMENT_R = r'''
    library(GO.db)
    library(stringr)
    library(plyr)
    library(ggplot2)
    library(AnnotationDbi)
    library(clusterProfiler)

    plotTargetGenesGoEnrichment <- function(targetGeneFiles, sampleNames, gene2goFile, outName, cutoff, filterBy, showCategory){
        targetGeneFiles <- unlist(strsplit(targetGeneFiles, ","))
        sampleNames <- unlist(strsplit(sampleNames, ","))
        gene2go <- read.delim(gene2goFile, header=T, sep="\t")
        go2term <- unique(toTable(GOTERM)[,c(1,3)])
        go2gene <- gene2go[,c(2,1)]
        names(targetGeneFiles) <- sampleNames
        outPdf <- paste0(outName, ".goEnrichResults.pdf")   # legacy passed width/height here, which paste0 appended to the name (BUG-170)
        if(length(sampleNames)>1){
            r_bind <- data.frame()
            for (i in seq(length(targetGeneFiles))){
                sample <- names(targetGeneFiles[i])
                genes <- read.csv(targetGeneFiles[i], header=FALSE)
                genes <- as.vector(genes[,1])
                goRes <- enricher(genes, TERM2GENE = go2gene, TERM2NAME = go2term, pvalueCutoff = 1, qvalueCutoff=1, maxGSSize = 100000)
                goResult <- as.data.frame(goRes)
                goResult <- goResult[which(goResult[filterBy]<=cutoff),]
                goResult$Ontology <- Ontology(as.vector(goResult$ID))
                goResult <- goResult[c(1, ncol(goResult), 2:(ncol(goResult)-1))]
                goResult$Ontology <- revalue(goResult$Ontology, c("BP"="Biological process", "MF"="Molecular function", "CC"="Cellular component"))
                goResult <- goResult[order(goResult$Ontology, goResult$pvalue),]
                goResult$Description <- factor(goResult$Description, levels = as.vector(goResult$Description))
                goResult$Cluster <- sample
                goResult$group <- sample
                goResult$level <- sample
                r_bind <- rbind(r_bind, goResult)
            }

            r_bind <- na.omit(r_bind)
            # one table for every sample (legacy rewrote it inside the loop, keeping only the last sample: BUG-172)
            write.table(r_bind[, setdiff(names(r_bind), c("group", "level"))], file=paste0(outName, ".goEnrichResults.txt"),
                        sep = "\t", row.names = FALSE, col.names = TRUE, quote = F)

            r_bind$level <- factor(r_bind$level, levels=sampleNames)

            r_bind_new <- new("compareClusterResult", compareClusterResult = r_bind)
            pdf(outPdf, width=10, height=8)
            # x as a string: current enrichplot rejects the legacy formula x=~group (BUG-171)
            p <- dotplot(r_bind_new, showCategory=showCategory, x="group") + scale_color_continuous(low='purple', high='green') +
                scale_y_discrete(labels=function(x) str_wrap(x, width=100)) +
                scale_size(range=c(0, 5)) + ggplot2::facet_grid(~level)
            print(p)
            dev.off()
        }else{
            sample <- sampleNames
            genes <- read.csv(targetGeneFiles[[1]][1], header=FALSE)
            genes <- as.vector(genes[,1])
            goRes <- enricher(genes, TERM2GENE = go2gene, TERM2NAME = go2term, pvalueCutoff = 1, qvalueCutoff=0.05, maxGSSize = 100000)
            goResult <- as.data.frame(goRes)
            goResult <- goResult[which(goResult[filterBy]<=cutoff),]
            goResult$Ontology <- Ontology(as.vector(goResult$ID))
            goResult <- goResult[c(1, ncol(goResult), 2:(ncol(goResult)-1))]
            goResult$Ontology <- revalue(goResult$Ontology, c("BP"="Biological process", "MF"="Molecular function", "CC"="Cellular component"))
            goResult <- goResult[order(goResult$Ontology, goResult$pvalue),]
            goResult$Description <- factor(goResult$Description, levels = as.vector(goResult$Description))
            outFile <- paste0(outName, ".goEnrichResults.txt")
            write.table(goResult, file=outFile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = F)
            pdf(outPdf, width=10, height=8)
            p <- dotplot(goRes, showCategory=showCategory) + scale_color_continuous(low='purple', high='green') +
                scale_y_discrete(labels=function(x) str_wrap(x, width=100)) + scale_size(range=c(0, 5))
            print(p)
            dev.off()
        }

    }

'''


def _split_csv(value: str | Iterable[str], name: str) -> list[str]:
    if isinstance(value, str):
        items = value.split(",")
    else:
        items = [str(item) for item in value]
    items = [item.strip() for item in items if item and item.strip()]
    if not items:
        raise ValueError(f"{name} must not be empty")
    return items


def _resolve_existing_file(path: str | Path, label: str) -> str:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    return str(resolved)


def run_go_enrichment(
    target_gene_files: str | Iterable[str],
    sample_names: str | Iterable[str],
    gene2go_file: str | Path,
    *,
    out_name: str = "goEnrichment",
    cutoff: float = 0.05,
    filter_by: str = "p.adjust",
    show_category: int = 20,
    out_dir: str | Path = ".",
    runner: ToolRunner | None = None,
    r_function_runner: Callable | None = None,
) -> dict[str, str]:
    """Run legacy GO enrichment plotting through ``Rscript``.

    The R function source is copied verbatim from ``plotRscriptStrs.py``; this
    wrapper only validates paths and builds the argument contract.
    """
    if filter_by not in {"pvalue", "p.adjust"}:
        raise ValueError("filter_by must be one of: pvalue, p.adjust")

    targets = [
        _resolve_existing_file(path, "target gene file")
        for path in _split_csv(target_gene_files, "target_gene_files")
    ]
    samples = _split_csv(sample_names, "sample_names")
    if len(targets) != len(samples):
        raise ValueError("target gene file count must match sample name count")
    gene2go = _resolve_existing_file(gene2go_file, "gene2go file")

    go_dir = Path(out_dir).expanduser().resolve() / "GO"
    go_dir.mkdir(parents=True, exist_ok=True)
    r_runner = r_function_runner or run_r_function
    kwargs = {"cwd": str(go_dir)}
    if runner is not None:
        kwargs["runner"] = runner
    r_runner(
        PLOT_TARGET_GENES_GO_ENRICHMENT_R,
        "plotTargetGenesGoEnrichment",
        [
            ",".join(targets),
            ",".join(samples),
            gene2go,
            out_name,
            float(cutoff),
            filter_by,
            int(show_category),
        ],
        **kwargs,
    )

    return {
        "go_dir": str(go_dir),
        "result_table": str(go_dir / f"{out_name}.goEnrichResults.txt"),
        "result_pdf": str(go_dir / f"{out_name}.goEnrichResults.pdf"),
    }
