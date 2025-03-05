import subprocess
from pathlib import Path


def draw_plot(sample, count, expression, de, output, map=None):
    Rscript = Path(__file__).parent.absolute() / "draw_plot_R" / "huatu.R"
    if map :
        cmd = f"Rscript {Rscript} -s {sample} -c {count} -e {expression} -d {de} -o {output} -m {map}"
    else:
        cmd = f"Rscript {Rscript} -s {sample} -c {count} -e {expression} -d {de} -o {output}"

    subprocess.run(cmd, shell=True, check=True)

def go_plot(sample, output):
    Rscript = Path(__file__).parent.absolute() / "draw_plot_R" / "GO_KEGG.R"
    cmd = f"Rscript {Rscript} -s {sample} -o {output}"
    subprocess.run(cmd, shell=True, check=True)

def run_plot(Rscript, sample, count, expression, de, output, map=None):
    if Rscript == "go":
        go_plot(sample, output)
    if Rscript == "heatmap":
        draw_plot(sample, count, expression, de, output, map)