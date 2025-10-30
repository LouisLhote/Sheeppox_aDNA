#!/usr/bin/env python3
"""
fetch_ncbi_metadata_full.py

Input:
  - accessions.txt : one GenBank accession per line

Output:
  - metadata.csv : columns: accession,collection_date,geo_loc_name,warning,note

Logic:
  1. Use /collection_date if present (full date preserved)
  2. Fallback: parse /note for year or decade
  3. Location: /geo_loc_name or /country, fallback to parsing note
  4. Add warnings for inferred info
"""

from Bio import Entrez, SeqIO
import csv
import re

Entrez.email = "your_email@example.com"  # REQUIRED

input_file = "accessions.txt"
output_file = "metadata.csv"

# simple list of countries for note parsing (expand if needed)
COUNTRIES = [
    "Turkey","Egypt","China","Kenya","Israel","Iran","Russia","India","Sudan","Uganda","Bulgaria","Greece"
]

def infer_location(note):
    for country in COUNTRIES:
        if re.search(r"\b" + re.escape(country) + r"\b", note, re.IGNORECASE):
            return country
    return "NA"

with open(input_file) as f:
    accessions = [line.strip() for line in f if line.strip()]

with open(output_file, "w", newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["accession","collection_date","geo_loc_name","warning","note"])
    
    for acc in accessions:
        collection_date = "NA"
        geo_loc = "NA"
        warning = ""
        note_text = ""
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            for f in record.features:
                if f.type == "source":
                    note_text = f.qualifiers.get("note", [""])[0]
                    
                    # COLLECTION DATE
                    cd = f.qualifiers.get("collection_date", ["NA"])[0]
                    if cd != "NA":
                        collection_date = cd  # keep full date
                    else:
                        # fallback: parse note
                        m = re.search(r"(\d{4})", note_text)
                        if m:
                            collection_date = m.group(1)
                            warning += "Year inferred from note; "
                        else:
                            m2 = re.search(r"(early|mid|late)?\s*(\d{4})'s", note_text, re.IGNORECASE)
                            if m2:
                                collection_date = m2.group(2)
                                warning += "Year inferred from note (decade); "
                            else:
                                warning += "No date found; "
                    
                    # LOCATION
                    geo_loc = f.qualifiers.get("geo_loc_name", f.qualifiers.get("country", ["NA"]))[0]
                    if geo_loc == "NA":
                        inferred_loc = infer_location(note_text)
                        if inferred_loc != "NA":
                            geo_loc = inferred_loc
                            warning += "Location inferred from note; "
                    
                    break  # only first source feature
            
            writer.writerow([acc, collection_date, geo_loc, warning.strip(), note_text])
            print(f"{acc}\t{collection_date}\t{geo_loc}\t{warning.strip()}")
            
        except Exception as e:
            print(f"Error fetching {acc}: {e}")
            writer.writerow([acc,"ERROR","ERROR","ERROR",""])
