# Define expected ground truth TaxIDs
# You have to get the IDs from /data/qix007/newDBkraken2/seqid2taxid.map, or wherever you install your database.
# In the future, I will work on a script that automates this process and extracts the IDs.
GROUND_TRUTH_TAXIDS = {
    "NC_000913.3": "511145",  # E. coli K-12 substr. MG1655
    "AF013254.1": "11250",    # genome name
    "ON811098.1": "2697049",  # genome name
    "OQ557947.1": "10244"     # genome name
}

# Function to load Kraken classification report from a file
def load_kraken_report(filename):
    taxon_counts = {}  # Dictionary to store read counts assigned to each TaxID
    total_reads = 0  # Variable to store the total number of classified reads

    # Open the Kraken report file for reading
    with open(filename, "r") as f:
        for line in f:  # Iterate through each line in the file
            parts = line.strip().split("\t")  # Split the line into tab-separated fields
            if len(parts) < 6:  # Ensure the line has at least 6 fields (expected Kraken output format)
                continue  # Skip malformed lines
            
            # Extract relevant fields from the Kraken output
            percent, total_clade, assigned_reads, rank, tax_id, tax_name = parts
            tax_id = tax_id.strip()  # Remove any leading or trailing whitespace from TaxID
            assigned_reads = int(assigned_reads.strip())  # Convert assigned read count to integer
            
            # Ensure the tax_id consists of digits (valid numeric TaxID)
            if tax_id.isdigit():  
                taxon_counts[tax_id] = assigned_reads  # Store the assigned read count for this TaxID
                total_reads += assigned_reads  # Increment the total read count

    return taxon_counts, total_reads  # Return the dictionary of counts and total reads

# Function to compute classification accuracy
def compute_accuracy(kraken_counts, total_reads, expected_taxids):
    # Sum up the counts of correctly classified reads by checking expected TaxIDs
    correct_reads = sum(kraken_counts.get(tax_id, 0) for tax_id in expected_taxids.values())
    
    # Compute the number of classified reads (excluding unclassified reads with tax_id = "0")
    classified_reads = total_reads - kraken_counts.get("0", 0)
    
    # Compute accuracy as the fraction of correctly classified reads over all classified reads
    accuracy = correct_reads / classified_reads if classified_reads > 0 else 0
    return accuracy  # Return the computed accuracy

# Load Kraken classification report from the specified file
kraken_counts, total_reads = load_kraken_report("kraken_report1.txt")

# Compute the classification accuracy using the expected TaxIDs
accuracy = compute_accuracy(kraken_counts, total_reads, GROUND_TRUTH_TAXIDS)

# Print the computed accuracy, formatted to 4 decimal places
print(f"Classification Accuracy: {accuracy:.4f}")
