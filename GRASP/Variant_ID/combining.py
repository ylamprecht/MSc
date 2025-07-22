import csv
import os

class Combiner:
    def __init__(self, logger):
        self.logger = logger

    def combine_cohorts(self, classifications, cohorts):
        """
        Combines variant count data across cohorts into single CSVs per classification.
        Creates summary tables showing variant distribution across cohorts.

        Parameters:
            classifications (list): List of variant classifications to process.
            cohorts (list): List of cohort names (e.g., ['Cohort_1', ..., 'Cohort_5']).

        Returns:
            None: Generates combined cohort CSV files per classification.
        """

        def read_csv_to_dict(filepath):
            data = {}
            with open(filepath, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    data[row['Gene_Nucleotide']] = int(row['Sample_Count'])
            return data

        out_dir = "variant_classifications/all_cohorts"
        os.makedirs(out_dir, exist_ok=True)

        for classification in classifications:
            merged = {}
            # Read and merge data from each cohort
            for cohort in cohorts:
                fpath = os.path.join("variant_classifications", cohort, f"{cohort}_{classification}_gene_nucleotide.csv")
                if os.path.exists(fpath):
                    cohort_data = read_csv_to_dict(fpath)
                    for gn, count in cohort_data.items():
                        if gn not in merged:
                            merged[gn] = {c: 0 for c in cohorts}
                        merged[gn][cohort] = count
                else:
                    self.logger.warning(f"File missing: {fpath}")

            # Write combined variant counts
            combined_file = os.path.join(out_dir, f"combined_cohorts_{classification}_gene_nucleotide.csv")
            with open(combined_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=['Gene_Nucleotide'] + cohorts)
                writer.writeheader()
                for gn, counts in merged.items():
                    writer.writerow({'Gene_Nucleotide': gn, **counts})

            self.logger.info(f"Combined data written for {classification}.")
