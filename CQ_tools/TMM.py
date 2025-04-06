import argparse
import os
import subprocess
import pandas as pd


def define_rscript(matrix_file, output_file):
    r_script = f"""
library(edgeR)

MakeMatrix <- read.table("{matrix_file}", header=TRUE, row.names=1, check.names=FALSE)
MakeMatrix <- as.matrix(MakeMatrix)
MakeMatrix <- round(MakeMatrix)
Matrix_exp <- DGEList(counts=MakeMatrix, group=factor(colnames(MakeMatrix)))
Matrix_exp <- calcNormFactors(Matrix_exp)
Matrix_exp$samples$eff.lib.size <- Matrix_exp$samples$lib.size * Matrix_exp$samples$norm.factors
write.table(Matrix_exp$samples, file="{output_file}", quote=FALSE, sep="\\t", row.names=FALSE)
"""
    return r_script


def temp_rscript(matrix_file):
    tmm_info_txt = f"{matrix_file}.TMM_info.txt"
    r_script_run = define_rscript(matrix_file, tmm_info_txt)
    r_script_file = "__tmp_runTMM.R"

    with open(r_script_file, "w") as f:
        f.write(r_script_run)
    try:
        subprocess.run(["Rscript", r_script_file], check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error running R script: {e}")
    finally:
        os.remove(r_script_file)

    return tmm_info_txt

def normalize_matrix(matrix_file,tmm_info_txt):
    tmm_read = pd.read_csv(tmm_info_txt, sep="\t")
    norm_factors = dict(zip(tmm_read['group'], tmm_read['norm.factors']))
    matrix_table = pd.read_csv(matrix_file, sep="\t", index_col=0)
    normalized_matrix = matrix_table.apply(lambda norm: norm / norm_factors[norm.name], axis=0)
    normalized_outpath = f"{matrix_file}_TMM.txt"
    normalized_matrix.to_csv(normalized_outpath, sep="\t", header=True, index=True)


def process_tmm(matrix_file):
    """
    Perform TMM normalization on a given matrix file.
    :param matrix_file: Path to the matrix file.
    """
    if not os.path.exists(matrix_file):
        raise FileNotFoundError(f"Matrix file not found: {matrix_file}")
    tmm_info = temp_rscript(matrix_file)
    normalize_matrix(matrix_file, tmm_info)


def man():
    parser = argparse.ArgumentParser(description="Perform TMM normalization on RNA-Seq count data.")
    parser.add_argument("--matrix", required=True, help="Matrix of raw read counts (not normalized).")
    args = parser.parse_args()
    process_tmm(args.matrix)



if __name__ == "__main__":
    man()