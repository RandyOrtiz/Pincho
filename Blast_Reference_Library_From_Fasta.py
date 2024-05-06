#!/usr/bin/env python3
import tkinter as tk
from tkinter import filedialog
import subprocess
import os
from os import path
import sys

def create_blast_db():
    # Determine the type of database based on the user's selection
    db_type = 'nucl' if db_var.get() == 'Nucleotide' else 'prot'

    # Get the file path from the label (updated when selecting file)
    file_path = file_label.cget("text")
    bin_dir = os.path.join(sys.path[0], 'bin')
    makeblastdb_dir = os.path.join(bin_dir, "ncbi-blast-2.3.0+-x64-linux/ncbi-blast-2.3.0+/bin/makeblastdb")

    # Execute the command
    try:
        subprocess.run(f"{makeblastdb_dir} -in {file_path} -dbtype {db_type} -parse_seqids", shell=True)
        result_label.config(text="Database created successfully!")
        root.after(2000, root.destroy)  # Wait for 2000 milliseconds before closing the window
    except subprocess.CalledProcessError as e:
        result_label.config(text=f"An error occurred: {e}")
        root.after(2000, root.destroy)  # Wait for 2000 milliseconds before closing the window

def select_file():
    filepath = filedialog.askopenfilename()
    if filepath:
        file_label.config(text=filepath)

# Create the main tk window
root = tk.Tk()
root.title("BLAST Database Creator")

frame = tk.Frame(root)
frame.pack(padx=10, pady=10)

db_var = tk.StringVar(value="Nucleotide")
tk.Radiobutton(frame, text="Nucleotide", variable=db_var, value="Nucleotide").grid(row=0, column=0)
tk.Radiobutton(frame, text="Protein", variable=db_var, value="Protein").grid(row=0, column=1)

tk.Button(frame, text="Select File", command=select_file).grid(row=1, column=0, columnspan=2)

file_label = tk.Label(frame, text="No file selected", relief="sunken", width=40, anchor="w")
file_label.grid(row=2, column=0, columnspan=2, pady=5)

tk.Button(frame, text="Create Database", command=create_blast_db).grid(row=3, column=0, columnspan=2)

result_label = tk.Label(frame, text="")
result_label.grid(row=4, column=0, columnspan=2, pady=5)

root.mainloop()