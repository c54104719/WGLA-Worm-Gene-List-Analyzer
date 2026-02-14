# C:\Users\yukan\Downloads\COSBI_Lab\enrich_project\HW4PreprocessRawData\f5_add_go_terms.py
import csv
from pathlib import Path


def parse_obo_file(obo_file):
    """
    Parse the OBO file and extract GO ID -> GO term name mapping.
    
    Args:
        obo_file: Path to the gene_ontology.WS298.obo file
    
    Returns:
        Dictionary mapping GO ID to GO term name
    """
    go_terms = {}
    current_go_id = None
    current_name = None
    
    with open(obo_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            
            # Check if we're starting a new term section
            if line == '[Term]':
                # Save previous term if we had one
                if current_go_id is not None and current_name is not None:
                    go_terms[current_go_id] = current_name
                current_go_id = None
                current_name = None
            
            # Extract GO ID
            elif line.startswith('id: '):
                current_go_id = line[4:]  # Remove 'id: ' prefix
            
            # Extract GO term name
            elif line.startswith('name: '):
                current_name = line[6:]  # Remove 'name: ' prefix
        
        # Don't forget the last term
        if current_go_id is not None and current_name is not None:
            go_terms[current_go_id] = current_name
    
    return go_terms


def add_go_terms_to_csv(f4_csv, obo_file, f5_csv):
    """
    Add GO term column to the CSV file.
    
    Args:
        f4_csv: Input CSV file (f4_go_annotation_processed.csv)
        obo_file: OBO file with GO terms (gene_ontology.WS298.obo)
        f5_csv: Output CSV file (f5_...)
    """
    print("Parsing OBO file...")
    go_terms = parse_obo_file(obo_file)
    print(f"Loaded {len(go_terms)} GO terms from OBO file")
    
    print(f"Processing {f4_csv}...")
    with open(f4_csv, 'r', encoding='utf-8') as infile, \
         open(f5_csv, 'w', newline='', encoding='utf-8') as outfile:
        reader = csv.DictReader(infile)
        
        # Add GO_term to fieldnames
        fieldnames = reader.fieldnames + ['GO_term']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        
        # Write header
        writer.writeheader()
        
        # Process each row
        rows_processed = 0
        rows_with_term = 0
        
        for row in reader:
            go_id = row['GO ID']
            
            # Look up the GO term
            if go_id in go_terms:
                row['GO_term'] = go_terms[go_id]
                rows_with_term += 1
            else:
                row['GO_term'] = 'Unknown'
            
            writer.writerow(row)
            rows_processed += 1
            
            if rows_processed % 10000 == 0:
                print(f"  Processed {rows_processed} rows...")
        
        print(f"Total rows processed: {rows_processed}")
        print(f"Rows with GO term found: {rows_with_term}")
        print(f"Rows without GO term (marked as 'Unknown'): {rows_processed - rows_with_term}")
    
    print(f"Output saved to {f5_csv}")


if __name__ == '__main__':
    # File paths
    current_dir = Path(__file__).parent
    f4_csv = current_dir / 'f4_go_annotation_processed.csv'
    obo_file = current_dir / 'gene_ontology.WS298.obo'
    f5_csv = current_dir / 'f5_go_annotation_with_terms.csv'
    
    # Validate input files
    if not f4_csv.exists():
        print(f"Error: {f4_csv} not found!")
        exit(1)
    if not obo_file.exists():
        print(f"Error: {obo_file} not found!")
        exit(1)
    
    # Process
    add_go_terms_to_csv(str(f4_csv), str(obo_file), str(f5_csv))
