"""
grasp.py

GRASP (Genetic Risk factors Associated with Severe Phenotypes) pipeline.
Main script for running variant identification across cohorts.
"""

import os
import logging
from cleaning import Cleaner
from counting import Counter
from combining import Combiner
from candidate import CandidateIdentifier

if __name__ == "__main__":
    # Set up logging: logs to both a file (grasp.log) and standard output.
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler("grasp.log")
        ]
    )
    logger = logging.getLogger("GRASP")

    # Define cohort raw input folders inside Frank   lin_output
    cohort_folders = [os.path.join('Franklin_output', f'Cohort_{i}_raw') for i in range(1, 6)]
    # Define cleaned output folders inside Franklin_output
    cohort_outputs = [os.path.join('Franklin_output', f'Cohort_{i}') for i in range(1, 6)]


    # Initialise cleaning and counting classes with logger
    cleaner, counter = Cleaner(logger), Counter(logger)

    # Loop through each cohort's raw data folder
    for infolder, outfolder in zip(cohort_folders, cohort_outputs):
        files = os.listdir(infolder)
        os.makedirs(outfolder, exist_ok=True)

        # Skip empty cohort folders
        if not files:
            logger.warning(f"Skipping {infolder}: folder empty.")
            continue

        # Find unique sample names from filenames (before _single_snp_variants)
        samples = {f.split('_single_snp_variants')[0] for f in files if 'single_snp_variants' in f}
        if not samples:
            logger.warning(f"No valid samples in {infolder}.")
            continue

        # Clean and merge both input CSVs for each sample
        for sample in samples:
            f_default = os.path.join(infolder, f"{sample}_single_snp_variants.csv")
            f_utr = os.path.join(infolder, f"{sample}_single_snp_variants (1).csv")
            f_out = os.path.join(outfolder, f"{sample}.csv")
            cleaner.clean_franklin(f_default, f_utr, f_out)

        # Count variants after cleaning, grouped by classification
        counter.count_variants(outfolder)

    # Define variant classifications of interest
    classifications = [
        "BENIGN", "LIKELY_BENIGN", "LIKELY_PATHOGENIC", "PATHOGENIC",
        "POSSIBLY_BENIGN", "POSSIBLY_PATHOGENIC_LOW", "POSSIBLY_PATHOGENIC_MODERATE", "UNCERTAIN_SIGNIFICANCE"
    ]
    # Cohort names for combining data
    cohorts = [f'Cohort_{i}' for i in range(1, 6)]
    # Total sample counts for each cohort (used to calculate proportions)
    cohort_sizes = {'Cohort_1': 18, 'Cohort_2': 15, 'Cohort_3': 15, 'Cohort_4': 13, 'Cohort_5': 23}

    # Combine cohort variant counts into unified tables per classification
    combiner = Combiner(logger)
    combiner.combine_cohorts(classifications, cohorts)

    # Identify candidate variants with prop_diff > 0.3 and significant Fisher p-values
    candidate = CandidateIdentifier(logger)
    candidate.identify_candidates(classifications, cohort_sizes)

    logger.info("GRASP pipeline complete.")
