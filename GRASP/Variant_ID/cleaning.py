import csv

class Cleaner:
    def __init__(self, logger):
        self.logger = logger

    def clean_franklin(self, input_file, output_file):
        """
        Cleans a Franklin variant CSV file by keeping only relevant columns
        and removing unnecessary quotation marks. Writes the cleaned data to a new CSV.

        Parameters:
            input_file (str): Path to the input CSV file.
            output_file (str): Path to save the cleaned output CSV file.

        Returns:
            None: Writes the cleaned data to the specified output file.
        """

        # Columns to keep in the output file for easier analysis
        columns_to_keep = ["Gene", "Nucleotide", "Genoox Classification", "Zygosity"]

        cleaned_rows = []

        with open(input_file, 'r', encoding='utf-8') as infile:
            reader = csv.DictReader(infile)
            for row in reader:
                cleaned_row = {}
                for col in columns_to_keep:
                    if col in row:
                        raw_value = row[col]
                        value = raw_value.replace('""', '').strip('"') if raw_value else ''
                        cleaned_row[col] = value
                cleaned_rows.append(cleaned_row)

        with open(output_file, 'w', newline='', encoding='utf-8') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=columns_to_keep)
            writer.writeheader()
            writer.writerows(cleaned_rows)

        self.logger.info(f"Cleaning complete for {output_file}.")
