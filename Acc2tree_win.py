import csv
import time
import os
import warnings
import subprocess
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO
from Bio.SeqFeature import BiopythonParserWarning
import tkinter as tk
from tkinter import simpledialog, ttk

def perform_blast_search(accession_id, max_target_seqs, email):
    Entrez.email = email
    print("Performing BLAST search...")
    result_handle = NCBIWWW.qblast(
        "blastp",
        "nr",
        accession_id,
        hitlist_size=max_target_seqs,
        expect=10,
        matrix_name="BLOSUM62",
    )
    print("BLAST search completed.")
    return result_handle

def retrieve_protein_sequence(accession):
    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text", timeout=30)
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        print(f"Error retrieving protein sequence for {accession}: {e}")
        return None

def perform_RAxML(directory_path):
    print("\nRAxML is making the tree...")
    try:
        subprocess.run(['Raxml', "-s", f'{directory_path}\\aligned_output_file', "-w", f'{directory_path}', "-n", "output", "-m", "PROTGAMMAAUTO", "-f", "a", "-x", "1000", "-p", "1988", "-k", "-N", "autoMRE"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        print('\r' + 'Making the tree completed.', end="", flush=True)
    except subprocess.CalledProcessError:
        print("Making the tree failed.")

def perform_IQtree(directory_path):
    print("\nIQtree is making the tree...")
    try:
        subprocess.run(['iqtree', "-s", f'{directory_path}\\aligned_output_file', "-m", "MFP", "-nt", "AUTO", "-bb", "1000", "-alrt", "1000", "-quiet"])
        print('\r' + 'Making the tree completed.', end="", flush=True)
    except subprocess.CalledProcessError:
        print("Making the tree failed.")

def run_analysis():
    email = email_var.get()
    accession_id = accession_var.get()
    max_target_seqs = max_target_seqs_var.get()
    preferred_tree = tree_algorithm_var.get()

    result_handle = perform_blast_search(accession_id, max_target_seqs, email)

    blast_records = NCBIXML.parse(result_handle)
    selected_records = []
    fasta_sequence_dict = {}

    for blast_record in blast_records:
        try:
            taxon_records = {}
            for alignment in blast_record.alignments:
                accession = alignment.accession
                e_value = alignment.hsps[0].expect
                definition = alignment.hit_def
                scientific_name = definition.split("[")[1].split("]")[0]
                taxon = scientific_name.lower()

                if accession.startswith("WP_"):
                    if taxon not in taxon_records:
                        taxon_records[taxon] = [(accession, e_value)]
                    else:
                        taxon_records[taxon].append((accession, e_value))
                else:
                    if taxon not in taxon_records:
                        taxon_records[taxon] = [(accession, e_value)]
                        selected_records.append((taxon, accession, e_value))
                    else:
                        if e_value < taxon_records[taxon][0][1]:
                            taxon_records[taxon][0] = (accession, e_value)
                            selected_records = [rec for rec in selected_records if rec[0] != taxon]

            for taxon, records in taxon_records.items():
                if len(records) == 1:
                    selected_records.append((taxon, records[0][0], records[0][1]))

        except ValueError as e:
            print(f"Error processing record: {e}")
            continue

    if selected_records:
        selected_records.sort(key=lambda x: (x[1].startswith("WP_"), x[2]))

        fasta_records = {}

        for i, record in enumerate(selected_records):
            protein_record = retrieve_protein_sequence(record[1])
            if protein_record:
                fasta_records[record[1]] = protein_record.format("fasta")
                fasta_sequence_dict[record[1]] = protein_record.seq
                organism = protein_record.annotations.get('organism', '')

            print(f"\rRecord {i + 1}/{len(selected_records)} processed", end="", flush=True)

        directory_path = f'C:\\ACC2TREE\\{accession_id}'
        if not os.path.exists(directory_path):
            os.makedirs(directory_path)

        with open(f'{directory_path}\\combined_fasta_sequences.fasta', "w") as fasta_file:
            for accession, fasta_record in fasta_records.items():
                fasta_file.write(fasta_record)

        print("\nCombined FASTA sequences saved to combined_fasta_sequences.fasta")

        print("Writing results to CSV...")
        with open(f'{directory_path}\\blast_results.csv', "w", newline="") as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(["DEFINITION", "Accession ID", "e value", "FASTA Sequence", "Link"])
            for i in range(len(selected_records)):
                taxon, accession, e_value = selected_records[i]
                fasta_sequence = fasta_sequence_dict.get(accession, "N/A")
                link = f"https://www.ncbi.nlm.nih.gov/protein/{accession}"
                csv_writer.writerow([taxon, accession, e_value, fasta_sequence, link])
                print(f"\rRecord {i + 1}/{len(selected_records)} written to CSV", end="", flush=True)
                time.sleep(1)

        try:
            print("\nMuscle is Aligning the sequences...")
            subprocess.run(['muscle', '-align', f'{directory_path}\\combined_fasta_sequences.fasta', '-output',
                            f'{directory_path}\\aligned_output_file'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            print('\r' + 'Alignment completed.', end="", flush=True)
        except subprocess.CalledProcessError:
            print("Alignment failed.")

    if preferred_tree == "RAxML":
        perform_RAxML(directory_path)
    else:
        perform_IQtree(directory_path)

    print("\nProcess completed.")

# Create main Tkinter window
root = tk.Tk()
root.title("Protein Analysis Tool")

# Create and place widgets in the window
email_label = tk.Label(root, text="Email:")
email_var = tk.StringVar()
email_entry = tk.Entry(root, textvariable=email_var)

accession_label = tk.Label(root, text="Accession ID:")
accession_var = tk.StringVar()
accession_entry = tk.Entry(root, textvariable=accession_var)

max_target_seqs_label = tk.Label(root, text="Maximum Target Sequences:")
max_target_seqs_var = tk.IntVar()
max_target_seqs_entry = tk.Entry(root, textvariable=max_target_seqs_var, validate="key", validatecommand=(root.register(lambda char: char.isdigit()), '%S'))

tree_algorithm_label = tk.Label(root, text="Tree Algorithm:")
tree_algorithm_var = tk.StringVar()
tree_algorithm_combobox = ttk.Combobox(root, textvariable=tree_algorithm_var, values=["RAxML", "IQtree"])
tree_algorithm_combobox.set("RAxML")

run_button = tk.Button(root, text="Run Analysis", command=run_analysis)

email_label.grid(row=0, column=0, padx=5, pady=5, sticky="w")
email_entry.grid(row=0, column=1, padx=5, pady=5, sticky="w")

accession_label.grid(row=1, column=0, padx=5, pady=5, sticky="w")
accession_entry.grid(row=1, column=1, padx=5, pady=5, sticky="w")

max_target_seqs_label.grid(row=2, column=0, padx=5, pady=5, sticky="w")
max_target_seqs_entry.grid(row=2, column=1, padx=5, pady=5, sticky="w")

tree_algorithm_label.grid(row=3, column=0, padx=5, pady=5, sticky="w")
tree_algorithm_combobox.grid(row=3, column=1, padx=5, pady=5, sticky="w")

run_button.grid(row=4, column=0, columnspan=2, pady=10)

# Start Tkinter main loop
root.mainloop()
