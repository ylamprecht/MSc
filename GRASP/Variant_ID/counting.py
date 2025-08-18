import csv
import os
from collections import defaultdict

class Counter:
    def __init__(self, logger):
        self.logger = logger

    def count_variants(self, path_to_csv_files):
        """
        Counts and categorises variants by classification for each cohort.
        For each classification, groups data by Gene_Nucleotide, creating summary CSV files
        for further analysis.

        Parameters:
            path_to_csv_files (str): Path to the directory containing CSV files for the cohort samples.

        Returns:
            None: Generates CSV files summarising variant counts by Gene_Nucleotide.
        """

        cohort_name = os.path.basename(path_to_csv_files)

        # Variant classifications of interest
        classifications = [
            'PATHOGENIC', 'LIKELY_PATHOGENIC', 'POSSIBLY_PATHOGENIC_MODERATE', 'POSSIBLY_PATHOGENIC_LOW', 'UNCERTAIN_SIGNIFICANCE'
        ]

        gene_nuc_data = {cls: defaultdict(list) for cls in classifications}

        # Process each sample CSV in cohort directory
        for file in os.listdir(path_to_csv_files):
            if file.endswith('.csv'):
                sample_name = file.replace('.csv', '')
                filepath = os.path.join(path_to_csv_files, file)
                with open(filepath, 'r') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        gene, nuc = row['Gene'], row['Nucleotide']
                        classification = row['Genoox Classification']
                        gene_nuc = f"{gene}_{nuc}"

                        if classification in gene_nuc_data:
                            gene_nuc_data[classification][gene_nuc].append(sample_name)

        # Output variant counts to CSV
        for classification, data in gene_nuc_data.items():
            output = [{
                'Gene_Nucleotide': gn,
                'Sample_Count': len(samples),
                'Samples': ', '.join(samples)
            } for gn, samples in data.items()]
            output.sort(key=lambda x: x['Sample_Count'], reverse=True)

            out_dir = os.path.join('variant_classifications', cohort_name)
            os.makedirs(out_dir, exist_ok=True)
            out_file = os.path.join(out_dir, f"{cohort_name}_{classification}_gene_nucleotide.csv")
            with open(out_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=['Gene_Nucleotide', 'Sample_Count', 'Samples'])
                writer.writeheader()
                writer.writerows(output)

        self.logger.info(f"Variant counting complete for {cohort_name}.")
