import csv
import os
from scipy.stats import fisher_exact

class CandidateIdentifier:
    def __init__(self, logger):
        self.logger = logger

    def identify_candidates(self, classifications, cohort_sizes):
        """
        Identifies candidate variants with proportional difference > 0.3
        and statistically significant Fisher p-values (p < 0.05).

        Parameters:
            classifications (list): List of variant classifications to process.
            cohort_sizes (dict): A dictionary mapping cohort names (str) to their sizes (int).

        Returns:
            None: Generates a candidate_variants.csv file with prioritised variants.
        """

        # Threshold for selecting potential candidate variants
        threshold = 0.3

        candidates = []
        for classification in classifications:
            combined_file = os.path.join("variant_classifications/all_cohorts", f"combined_cohorts_{classification}_gene_nucleotide.csv")
            if not os.path.exists(combined_file):
                continue

            with open(combined_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    try:
                        c1 = int(row['Cohort_1'])
                        c3 = int(row['Cohort_3'])
                        c5 = int(row['Cohort_5'])
                    except ValueError:
                        continue

                    total1_3 = cohort_sizes['Cohort_1'] + cohort_sizes['Cohort_3']
                    prop_5 = c5 / cohort_sizes['Cohort_5']
                    prop_1_3 = (c1 + c3) / total1_3
                    prop_diff = prop_5 - prop_1_3

                    if prop_diff > threshold:
                        # Create 2x2 table for Fisher's exact test
                        a = c5
                        b = cohort_sizes['Cohort_5'] - c5
                        c = c1 + c3
                        d = total1_3 - (c1 + c3)
                        _, p_value = fisher_exact([[a, b], [c, d]], alternative='two-sided')

                        # Filter by p-value < 0.05
                        if p_value < 0.05:
                            candidates.append({
                                'Gene_Nucleotide': row['Gene_Nucleotide'],
                                'Classification': classification,
                                'Cohort_1': c1,
                                'Cohort_2': row['Cohort_2'],
                                'Cohort_3': c3,
                                'Cohort_4': row['Cohort_4'],
                                'Cohort_5': c5,
                                'Proportional_Difference': round(prop_diff, 4),
                                'p_value': round(p_value, 5)
                            })

        # Sort candidate variants by classification, then descending proportional difference
        classification_order = {
            'PATHOGENIC': 1,
            'LIKELY_PATHOGENIC': 2,
            'POSSIBLY_PATHOGENIC_MODERATE': 3,
            'POSSIBLY_PATHOGENIC_LOW': 4,
            'UNCERTAIN_SIGNIFICANCE': 5,
            'POSSIBLY_BENIGN': 6,
            'LIKELY_BENIGN': 7,
            'BENIGN': 8
        }
        candidates.sort(
            key=lambda x: (classification_order.get(x['Classification'], 99), -x['Proportional_Difference'])
        )

        out_file = "candidate_variants.csv"
        with open(out_file, 'w', newline='') as f:
            fieldnames = [
                'Gene_Nucleotide', 'Classification', 'Cohort_1', 'Cohort_2', 'Cohort_3',
                'Cohort_4', 'Cohort_5', 'Proportional_Difference', 'p_value'
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(candidates)

        self.logger.info(f"Candidate variants identified and saved to {out_file}.")
