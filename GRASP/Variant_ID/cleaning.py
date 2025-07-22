import csv

class Cleaner:
    def __init__(self, logger):
        self.logger = logger

    def clean_franklin(self, input_default, input_utr, output_file):
        """
        Cleans and merges two Franklin input CSV files by removing unnecessary columns
        and duplicate rows. Writes the cleaned data to an output CSV file.

        Parameters:
            input_default (str): Path to the first input CSV file (default variants).
            input_utr (str): Path to the second input CSV file (UTR variants).
            output_file (str): Path to save the cleaned and merged output CSV file.

        Returns:
            None: Writes the cleaned data to the specified output file.
        """

        # Columns to keep in the output file for easier analysis
        columns_to_keep = ["Gene", "Nucleotide", "Genoox_Classification", "Zygosity", "Inheritance_Model"]
        merged_rows = []

        # Read and merge rows from both default and UTR variant files
        with open(input_default, 'r', encoding='utf-8') as infile1, open(input_utr, 'r', encoding='utf-8') as infile2:
            reader1, reader2 = csv.DictReader(infile1), csv.DictReader(infile2)
            merged_rows.extend(reader1)
            merged_rows.extend(reader2)

        # Write cleaned, merged rows to output CSV
        with open(output_file, 'w', newline='', encoding='utf-8') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=columns_to_keep)
            writer.writeheader()

            seen_entries, duplicates = set(), []

            # Filter duplicates by unique (Gene, Nucleotide) combination
            for row in merged_rows:
                uid = (row["Gene"], row["Nucleotide"])
                if uid in seen_entries:
                    duplicates.append(row)
                else:
                    seen_entries.add(uid)
                    cleaned_row = {k: row[k].replace('""', '').strip('"') for k in columns_to_keep}
                    writer.writerow(cleaned_row)

        # Log duplicate entries (if any)
        if duplicates:
            self.logger.warning(f"Duplicate entries found in {output_file}. Skipping duplicates...")
        else:
            self.logger.info(f"{output_file}: No duplicates found.")

        self.logger.info(f"Cleaning complete for {output_file}.")
