import csv
import time
import os
import warnings
import subprocess
import threading
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO
from Bio.SeqFeature import BiopythonParserWarning
import tkinter as tk
from tkinter import *
from tkinter import simpledialog, ttk, messagebox

def perform_blast_search(email, accession_id, max_target_seqs, Display):
    Entrez.email = email
    Display.insert(tk.END, "Performing BLAST search...\n")
    result_handle = NCBIWWW.qblast(
        "blastp",
        "nr",
        accession_id,
        hitlist_size=max_target_seqs,
        expect=10,
        matrix_name="BLOSUM62",
    )
    Display.insert(tk.END, "BLAST search completed.\n")
    return result_handle

def retrieve_protein_sequence(accession):
    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text", timeout=1)
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        Display.insert(tk.END, f"Error retrieving protein sequence for {accession}: {e} \n")
        return None
    
def perform_RAxML(directory_path, Display):            
    Display.insert(tk.END,"\nRAxML is making the tree... \n")
    try:
        subprocess.run(['Raxml', "-s", f'{directory_path}\\aligned_output_file', "-w", f'{directory_path}', "-n", "output", "-m", "PROTGAMMAAUTO", "-f", "a", "-x", "1000", "-p", "1988", "-k", "-N", "autoMRE"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True) 
        Display.insert(tk.END, 'Making the tree completed.')
    except subprocess.CalledProcessError as e:
        Display.insert(tk.END, f"Making the tree failed. Error: {e}\n")

        print(f"RAxML Error Output: {e.stderr}")

def main(Display):
    email = email_var.get()
    accession_id = accession_var.get()
    max_target_seqs = max_target_seqs_var.get()
    preferred_tree = tree_algorithm_entry.get()
    if not accession_var.get():
        messagebox.showerror("Error", "Accession ID cannot be empty.")
        return
    warnings.filterwarnings("ignore", category=BiopythonParserWarning, message="Dropping bond qualifier in feature location")
    Display.config(state=tk.NORMAL)
    result_handle = perform_blast_search(email, accession_id, max_target_seqs, Display)
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
            Display.insert(tk.END, f"Error processing record: {e}\n")
            continue

    if selected_records:
        selected_records.sort(key=lambda x: (x[1].startswith("WP_"), x[2]))

        fasta_sequences = []
        fasta_records = {}  
                
        for i, record in enumerate(selected_records):
            protein_record = retrieve_protein_sequence(record[1])
            if protein_record:
                fasta_sequences.append(protein_record.format("fasta"))  
                fasta_records[record[1]] = protein_record.format("fasta")  
                fasta_sequence_dict[record[1]] = protein_record.seq  
                organism = protein_record.annotations.get('organism', '')  

            Display.delete(1.0, tk.END)
            Display.insert(tk.END, f"Record {i+1}/{len(selected_records)} processed")
            Display.update_idletasks()
        
        Display.insert(tk.END, "\n")
        
        directory_path = f'C:\\ACC2TREE\\{accession_id}'
        if not os.path.exists(directory_path):
            os.makedirs(directory_path)
        with open(f'{directory_path}\\combined_fasta_sequences.fasta', "w") as fasta_file:
            for accession, fasta_record in fasta_records.items():
                fasta_file.write(fasta_record)

        Display.insert(tk.END, "\nCombined FASTA sequences saved to combined_fasta_sequences.fasta \n")
        
        Display.insert(tk.END, "Writing results to CSV...\n")
        with open(f'{directory_path}\\blast_results.csv', "w", newline="") as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(["DEFINITION", "Accession ID", "e value", "FASTA Sequence", "Link"])
            for i in range(len(selected_records)):
                taxon, accession, e_value = selected_records[i]
                fasta_sequence = fasta_sequence_dict.get(accession, "N/A") 
                link = f"https://www.ncbi.nlm.nih.gov/protein/{accession}"  
                csv_writer.writerow([taxon, accession, e_value, fasta_sequence, link])
                Display.delete(1.0, tk.END)
                Display.insert(tk.END, f"Record {i+1}/{len(selected_records)} written to CSV")
                Display.update_idletasks()
                time.sleep(1)

        Display.insert(tk.END, "\n")

        Display.insert(tk.END, "\nMuscle is Aligning the sequences... \n") 
        try:
            subprocess.run(['muscle', '-align', f'{directory_path}\\combined_fasta_sequences.fasta', '-output', f'{directory_path}\\aligned_output_file'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            Display.insert(tk.END, 'Alignment completed.')
        except subprocess.CalledProcessError:
            Display.insert(tk.END, "Alignment failed. \n")
            
    if tree_algorithm_var.get() == "RAxML":
        perform_RAxML(directory_path, Display)
    elif tree_algorithm_var.get() == "IQtree":
        Display.insert(tk.END, "\nIQtree is making the tree... \n")
        try:
            subprocess.run(['iqtree', "-s", f'{directory_path}\\aligned_output_file', "-m", "MFP", "-nt", "AUTO", "-bb", "1000", "-alrt", "1000", "-quiet"])
            Display.insert(tk.END, 'Making the tree completed.')
        except subprocess.CalledProcessError:
            Display.insert(tk.END, "Making the tree failed. Error: {e}\n")
    
    Display.insert(tk.END, "\nProcess completed. \n")        
    

def update_Display(Display, main_func):
    main_func(Display)
    Display.update_idletasks()

def main_with_update(Display):
    threading.Thread(target=main, args=(Display,)).start()
    Display.after(100, update_Display, Display, main)
    
# Tkinter Window
win = tk.Tk()
win.geometry = ("100x100")
win.resizable(width=False, height=False)
win.title("Acc2Tree")
frame = tk.Frame(win)
frame.pack()

# LabelFrames
first_frame_info = tk.LabelFrame(frame, text="Email")
first_frame_info.grid(row=0, column=0, sticky="news", padx=20, pady=(5, 2))
second_frame_info = tk.LabelFrame(frame, text="Blast Entries")
second_frame_info.grid(row=1, column=0, sticky="news", padx=20, pady=2)
third_frame_info = tk.LabelFrame(frame, text="Tree Algorithm")
third_frame_info.grid(row=2, column=0, sticky="news", padx=20, pady=2)
fourth_frame_info = tk.LabelFrame(frame, relief=FLAT)
fourth_frame_info.grid(row=4, column=0, sticky="news", padx=20, pady=(2, 5))

# Entry Labels
email_label = tk.Label(first_frame_info, text="Your Email:        ")
email_label.grid(row=0, column=0, sticky="w")
accession_label = tk.Label(second_frame_info, text="Accession Id.")
accession_label.grid(row=1, column=0, sticky="w")
max_target_seqs_label = tk.Label(second_frame_info, text="Max Target Seq.")
max_target_seqs_label.grid(row=2, column=0, sticky="w")
tree_algorithm_label = tk.Label(third_frame_info, text="Choose One:")
tree_algorithm_label.grid(row=3, column=0, sticky="w")

# Entries and Variables
email_var = tk.StringVar()
email_entry = tk.Entry(first_frame_info, textvariable=email_var,)

accession_var = tk.StringVar()
accession_entry = tk.Entry(second_frame_info, textvariable=accession_var)

max_target_seqs_var = tk.IntVar(value=100)
max_target_seqs_entry = tk.Entry(second_frame_info, textvariable=max_target_seqs_var, validate="key", validatecommand=(win.register(lambda char: char.isdigit()), '%S'))

tree_algorithm_var = tk.StringVar()
tree_algorithm_entry = ttk.Combobox(third_frame_info, textvariable=tree_algorithm_var, values=["RAxML", "IQtree"])
tree_algorithm_entry.set("RAxML")

email_entry.grid(row=0, column=1, sticky="w")
accession_entry.grid(row=1, column=1)
max_target_seqs_entry.grid(row=2, column=1)
tree_algorithm_entry.grid(row=3, column=1,sticky="w")

# Display Area
Display = tk.Text(fourth_frame_info, height=10, width=31, state=tk.NORMAL)
scrollbar = tk.Scrollbar(fourth_frame_info, command=Display.yview)
Display.config(yscrollcommand=scrollbar.set)
Display.grid(row=4, column=0, padx=(0, 0), pady=10, sticky="w")
scrollbar.grid(row=4, column=1, sticky="wns", padx=(0, 20), pady=10)

# Button
start_button = tk.Button(frame, text="Start", command=lambda: main_with_update(Display))
start_button.grid(row=3, column=0, padx=20, pady=5, sticky="news")

# Configure padding for child widgets
for widget in first_frame_info.winfo_children():
    widget.grid_configure(padx=10, pady=5)
for widget in second_frame_info.winfo_children():
    widget.grid_configure(padx=10, pady=5)
for widget in third_frame_info.winfo_children():
    widget.grid_configure(padx=(18, 10), pady=5)
for widget in fourth_frame_info.winfo_children():
    widget.grid_configure(padx=1, pady=5)
 
win.mainloop()
