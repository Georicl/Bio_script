import os
import subprocess


##################
# function about align re-sequence to genome
##################
def index_bwa(reference_genome, prefix="A"):
    idx_bwa = ["bwa", "index", "-p", prefix, reference_genome]
    try:
        subprocess.run(idx_bwa, check=True, text=True)
        print("BWA index finished")
    except subprocess.CalledProcessError as e:
        print("BWA index's error log path:", e.stderr)
        raise


def align_bwa(reference_prefix, sample_1, sample_2, bwa_output, cpu,):
    run_BWA = ["bwa", "mem", "-t", str(cpu), "-o", bwa_output, reference_prefix, sample_1, sample_2]
    try:
        subprocess.run(run_BWA, check=True, text=True)
        print("BWA finished output path:", bwa_output)
    except subprocess.CalledProcessError as e:
        print("BWA's error log path:", e.stderr)
        raise


def run_bwa(reference_genome, sample_1, sample_2, bwa_output, cpu, prefix="A"):
    index_bwa(reference_genome, prefix)
    align_bwa(prefix, sample_1, sample_2, bwa_output, cpu)


##################
# function about use samtools to convert sam to bam
##################
def samtools(bwa_output, output_bam, output_sort, output_tsv, open_primaryalignment, output_secondaryalignment):
    try:
        print(f"Converting {bwa_output} to BAM: {output_bam}")
        subprocess.run(["samtools", "view", "-bS", "-o", output_bam, bwa_output], check=True)

        print(f"Sorting BAM: {output_sort}")
        subprocess.run(["samtools", "sort", "-o", output_sort, output_bam], check=True)

        print("Indexing BAM")
        subprocess.run(["samtools", "index", output_sort], check=True)

        print("Generating flagstat report")
        with open(output_tsv, 'w') as f:
            subprocess.run(["samtools", "flagstat", output_sort], stdout=f, check=True)
        print("Samtools processing finished successfully!")

    except subprocess.CalledProcessError as e:
        print(f"Samtools error: {e}")
        raise


def run_bwa_samtools(fasta, sample_1, sample_2, output_dir, cpu):
    bwa_output = os.path.join(output_dir, "bwa_output.sam")
    output_bam = os.path.join(output_dir, "output.bam")
    output_sort = os.path.join(output_dir, "output.sort.bam")
    output_tsv = os.path.join(output_dir, "output.tsv")
    output_primaryalignment = os.path.join(output_dir, "output.primaryalignment.bam")
    output_secondaryalignment = os.path.join(output_dir, "output.secondaryalignment.bam")

    # Run BWA alignment and indexing
    run_bwa(fasta, sample_1, sample_2, bwa_output, cpu)

    # Convert SAM to BAM and sort/index using samtools
    samtools(bwa_output, output_bam, output_sort, output_tsv, output_primaryalignment, output_secondaryalignment)


if __name__ == "__main__":
    # Example usage (replace these with actual paths and parameters)
    bwa_output = "bwa_output.sam"
    output_bam = "test.bam"
    output_sort = "test.sort.bam"
    output_tsv = "test.tsv"
    samtools(bwa_output, output_bam, output_sort, output_tsv)
