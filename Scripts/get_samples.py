from collections import defaultdict
from typing import NamedTuple, Tuple, List, Dict
from sys import argv
from os.path import exists, join, split
from os import remove, rename
import random
import subprocess

#Leaf51, Axenic

Entry = NamedTuple(
    "Entry",
    [
        ("study_accession", str),
        ("sample_accession", str),
        ("experiment_accession", str),
        ("run_accession", str),      
        ("read_count", int),
        ("submitted_ftp", str),      
        ("sample_title", str)    
    ]
)

def parse_fasta(seq_filename: str, strands: dict) -> Dict[str, Tuple[str, int]]:
    """
    Take a fasta file and write the data into the dictionary. Provides name, sequence and length

    :param seq_filename: Name of the fasta file
    :param strands:  dictionary to write the information into for appending

    :returns Dictionary containing {name: (sequence, length), ...}
    """
    with open(seq_filename, "r", encoding="utf-8") as file:
        entries = [entry.split("\n") for entry in file.read().split(">")[1:]]

        for entry in entries:
            label = entry[0].split()[0]  # Obtain just the name
            if "|" in label:  # Some sequences can have a pipe, so we want only the last entry
                label = label.split("|")[-1]
            data = "".join(entry[1:])

            strands[label] = (data, len(data))

    return strands


def parse_arguments() -> Tuple[str, str]:
    """
    Parse the arguments from commandline to check if they are valid files

    :returns the filenames in order provided
    """
    try:
        filename = argv[1]
        if not exists(filename):
            raise IndexError
    except IndexError as e:
        raise ValueError("Missing filename for samples") from e

    try:
        samples = list(argv[2].split(","))
    except IndexError as e:
        raise ValueError("Missing argument samples to download") from e


    try:
        location = argv[3]
    except IndexError as e:
        raise ValueError("Missing argument folder") from e

    
    try:
        if len(argv) > 4:
            remove = float(argv[4])
        else: 
            remove = 0
    except IndexError as e:
        raise ValueError("Error parsing argument") from e

    return filename, samples, location, remove


def get_sample_entries(filename: str) -> Dict[str, List[Entry]]:
    with open(filename, "r", encoding="utf-8") as samples_file:
        entries = defaultdict(list)
        for line in samples_file.readlines():
            curr_entry = Entry(*line.strip().split(","))

            curr_sample = curr_entry.sample_title.split("_")[0]
            
            entries[curr_sample].append(curr_entry)
        
        return entries

def download_samples(entries: Dict[str, List[Entry]], samples_to_get: List[str], download_location: str) -> Dict[str, List[str]]:
    filenames = defaultdict(list)

    for sample in samples_to_get:
        if sample not in entries.keys():
            raise Exception(f"Provided sample not present in list: {sample}")
        
        runs = entries[sample]
        for run in runs:
            link = run.submitted_ftp
            filename = link.split("/")[-1]
            file_location = join(download_location, filename)
            
            if not exists(file_location):
                print(f"Would have downloaded {link} at {file_location}")
                #result = subprocess.check_call(f"wget {link} -P {download_location}", shell=True)
                result = 0
                if result != 0:
                    raise Exception(f"wget failed with error code {result}")
            
            filenames[sample].append(file_location)
         
    return filenames        
        
def trim_files(filenames: Dict[str, List[str]], trim_amount: float, limit=3):
    files = defaultdict(list)
    for sample, files in filenames.items():
        for index, file in enumerate(files):
            if index >= limit:
                continue
            
            filename = split(file)[1][:-3]
            target_location = join("..", "Intermediate", "to_trim_" + filename)
            trimmed_location = join("..", "Intermediate", "trimmed_" + filename)

            if not exists(trimmed_location):

                # Unzip file to temp dir
                if not exists(target_location) and exists(file):
                    print(f"Would have uncompressed {filename} at {file} to {target_location}")
                    #result = subprocess.check_call(f"gzip -dkc {file_location} > {target_location}", shell=True)
                    result = 0
                    if result != 0:
                        raise Exception(f"gzip failed with error code {result}")
                else:
                    raise Exception(f"Could not find file: {file}")


                # Remove percentage of entries
                # write file to proper location with _trimmed suffix
                fastas = parse_fasta(target_location, {})
                with open(trimmed_location, "w") as trimmed_file:
                    samples = random.sample(fastas.items(), round(trim_amount*len(fastas)))
                    for title, seq in samples:
                        trimmed_file.write(f">{title}\n{seq}\n")

                remove(target_location)
            else:
                print(f"Found trimmed file at {trimmed_location}")

            files[sample].append(trimmed_location)

def rename_files(filenames: Dict[str, List[str]]):
    for sample, files in filenames.items():
        for index, file in enumerate(files):
            new_name = join(split(file)[0], f"{sample}_{index}.fastq")
            print(f"renaming {file} to {new_name}")
            #rename(file, new_name)
    

def main():
    sample_filename, samples_to_get, location, remove_percentage = parse_arguments()
    entries = get_sample_entries(sample_filename)
    file_locations = download_samples(entries, samples_to_get, location)
    if remove_percentage > 0:
        file_locations = trim_files(file_locations, remove_percentage)
    rename_files(file_locations)
    
if __name__ == "__main__":
    main()
    