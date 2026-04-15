import concurrent.futures
import copy
import os
import re
import shutil
import subprocess
import tempfile
import time
from collections import Counter
from pathlib import Path

import numpy as np
from Bio import Align, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from dna_utils import rnatoDNA, revcomp
from ensemble_helper import find_fully_enclosed_segments
from idt_api import findPlaceholderSequences

aligner1 = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-1)
aligner2 = Align.PairwiseAligner(
    mode='global',
    match_score=2,
    mismatch_score=-1,
    open_gap_score=-2.1,
    extend_gap_score=-0.5
)

# Runtime working directory for UNAFold temp files.
# Override via the UNAFOLD_WORK_DIR environment variable if needed.
UNAFOLD_WORK_DIR = Path(
    os.environ.get("UNAFOLD_WORK_DIR", str(Path(tempfile.gettempdir()) / "unafold_runs"))
)


def editableProbes(miRNA, TARGET_DNA_LENGTH):
    miDNA = rnatoDNA(miRNA)
    placeholder = revcomp(miDNA) + "TTTT"
    probe = revcomp(placeholder)[:-9]
    probe_end = probe[-7:]
    probe = "AAAA" + revcomp(probe_end) + "TATT" + "A"*(TARGET_DNA_LENGTH - 34) + probe
    return probe[:(15+(TARGET_DNA_LENGTH - 34))]

def uneditableProbes(miRNA):
    miDNA = rnatoDNA(miRNA)
    placeholder = revcomp(miDNA) + "TTTT"
    probe = revcomp(placeholder)[:-9]
    probe_end = probe[-7:]
    probe = "AAAA" + revcomp(probe_end) + "TATTAA" + probe
    return probe[17:]

#UNAFold File Processing Functions

def make_fasta(actual_population, fasta_file):
    """
    This function takes in a list of sequences (as strings) and creates a large 
    fasta file of all sequences and saves it in the specified location.
    """
    with open(fasta_file, 'w') as f:
        for i, seq in enumerate(actual_population, start=1):
            f.write(f'>sequence_{i}\n')
            f.write(f'{seq}\n')


def idt_hairpin_batch(SALT_CONC, MAG_CONC, folder_name,
                      fasta_file="sequences.fasta",
                      use_tm_constraint=True,
                      penalize_off_targets=True):
    """
    This function runs the UNAFold perl script. The hybrid-ss-min script allows
    for the use of the mfold option which can provide additional, suboptimal 
    structures if the suboptimality percentage is high enough (in this case 50).
    Additionally, the code is run with 150mM [Na+] at 25C. We set a tmax=26 because
    although we don't use the file 26C generates, if you have tmax=25, it only
    provides the most optimal structure in the ct file for some reason.
    """
    
    # Make sure we're in the correct folder so output files land there
    os.chdir(folder_name)

    # Temperature range depends ONLY on whether we're using a Tm constraint
    if use_tm_constraint:
        t_min = 23
        t_max = 27
    else:
        t_min = 25
        t_max = 26

    # Build the command
    base_cmd = (
        rf"hybrid-ss-min -n DNA "
        rf"-N {SALT_CONC} -M {MAG_CONC} "
        rf"-t {t_min} --tmax={t_max} "
    )

    # Use mfold only when we care about off-target structures
    if penalize_off_targets:
        base_cmd += r"--mfold=50 "

    command = base_cmd + fasta_file

    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print("UNAFold error in idt_hairpin_batch:")
        print(result.stderr)


def ct_to_dot_bracket(ct_lines):
    """
    This function converts ct files (uniquely with 6 columns instead of 4) 
    to dot bracket notation.
    """
    sequence = ""
    dot_bracket = ""

    paired_indices = [0] * len(ct_lines)

    for line in ct_lines[1:]:
        parts = line.split()
        if len(parts) >= 6:
            index = int(parts[0])
            paired_index = int(parts[4])

            paired_indices[index] = paired_index

    for i, paired_index in enumerate(paired_indices[1:], 1):
        if paired_index == 0:
            dot_bracket += "."
        elif paired_index > i:
            dot_bracket += "("
        else:
            dot_bracket += ")"

    return dot_bracket

def parse_ct_file(ct_file = r'sequences.fasta.25.ct'):
    """
    When using a collective fasta file of multiple sequences, the UNAFold software 
    outputs a conjoined ct file with multiple structures per sequence. This function
    parses through the conjoined ct file, and using ct_to_dot_bracket, it converts
    each ct file to dot bracket notation. It creates a dictionary where each key is
    the sequence index number (already provided by the fasta file) and each value
    is a list of dot bracket notation structures for the given sequence.
    """
    sequences = []
    with open(ct_file, 'r') as f:
        lines = f.readlines()
        lines = lines[::-1]

    all_cts = {}
    current_sequence_ct = []
    beginning, end = True, False
    for line in lines:
        if line.strip().split("\t")[1].split(" ")[0] != "dG":
            current_sequence_ct = [line] + current_sequence_ct
            last_num = line.strip().split("\t")[0]
            beginning = False
        elif end == True and sequence_name != last_sequence:
            break
        else:
            current_sequence_ct = [line] + current_sequence_ct
            sequence_name = line.strip().split("\t")[2]
            dot_bracket = ct_to_dot_bracket(current_sequence_ct)
            try:
                if dot_bracket not in all_cts[sequence_name]:
                    all_cts[sequence_name].append(dot_bracket)
            except:
                all_cts[sequence_name] = [dot_bracket]
            
            current_sequence_ct = []
            
            if sequence_name == "sequence_1":
                end = True
                last_sequence = sequence_name

    all_cts = dict(reversed(list(all_cts.items())))
    return all_cts

def parse_dG_file_and_calc_MT(ct_file = "sequences.fasta.dG"):
    """
    When using a collective fasta file of multiple sequences, the UNAFold software 
    outputs a conjoined dG file with multiple dG calculations per sequence and temperature. 
    This function parses through the conjoined dG file, and calculates the melting
    temperature at 25C using a series of entropy and enthalpy calculations. It creates and
    returns a list of melting temperatures whose order corresponds to the order of the
    population.
    """
    
    def find_melting_temp(current_dG, starting_temp=24):
        entropies, enthalpies, temp = [], [], starting_temp
        for temp_index in range(len(current_dG))[1:-1]:
            entropies.append(-(current_dG[temp_index+1] - current_dG[temp_index-1])/2)
            enthalpies.append(current_dG[temp_index] - ((273.15 + temp)*(current_dG[temp_index+1] - current_dG[temp_index-1])/2))
            temp += 1
        melting_temps = []
        for i in range(len(entropies)):
            melting_temps.append((enthalpies[i]/entropies[i]) - 273.15)
        melting_temp = np.mean(melting_temps)
        return round(melting_temp, 1)
    
    sequences = []
    with open(ct_file, 'r') as f:
        lines = f.readlines()

    melting_temps = []
    current_dG = []
    beginning = True
    for line in lines[1:]:
        if line.strip().split("\t")[0] == "23" and not beginning:
            melting_temp = find_melting_temp(current_dG)
            melting_temps.append(melting_temp)
            current_dG = []
            try:
                current_dG.append(float(line.strip().split("\t")[1]))
            except:
                print("error1")
        else:
            try:
                current_dG.append(float(line.strip().split("\t")[1]))
            except:
                print("error2")
            beginning = False
    
    melting_temp = find_melting_temp(current_dG)
    melting_temps.append(melting_temp)
    
    return melting_temps


def initialize_dna_strand(STARTER, POPULATION_SIZE):
    """
    This function initializes a DNA pool by mutating the STARTER.
    Ensures that lists are not shared across population members.
    """
    starting_pool = []
    for i in range(POPULATION_SIZE):
        member = [mutate(copy.deepcopy(j), 0.8) for j in STARTER]
        starting_pool.append(copy.deepcopy(member))  # Ensures full independence
    return starting_pool




def mutate(member_portion, mutation_rate):
    """
    This function mutates the editable region of a given strand by either mutating 
    bases in the middle or adding or subtracting bases in the middle and end. It also 
    has the option of removing or adding back the last base of the ineditable region.
    """
    
    if isinstance(member_portion, list):
        member_portion = copy.deepcopy(member_portion)
        # Either subtract a base or add a base back to the end of the uneditable region
        if np.random.rand()*2 < mutation_rate:
            if member_portion[-1] == len(member_portion)-2:
                member_portion[-1] = 0
            else:
                member_portion[-1] += 1 
    
    else:
        # Mutate in the middle
        for i in range(len(member_portion)):
            if np.random.rand() < mutation_rate:
                member_portion = member_portion[:i] + np.random.choice(['A', 'T', 'G', 'C'], p=[0.4, 0.4, 0.1, 0.1]) + member_portion[i+1:]

        # Add one base to the end
        if np.random.rand()*1.2 < mutation_rate:
            member_portion += np.random.choice(['A', 'T', 'G', 'C'], p=[0.4, 0.4, 0.1, 0.1])

        # Add two bases to the end
        elif np.random.rand()*2 < mutation_rate:
            member_portion += (np.random.choice(['A', 'T', 'G', 'C'], p=[0.4, 0.4, 0.1, 0.1]) + np.random.choice(['A', 'T', 'G', 'C'], p=[0.7, 0.1, 0.1, 0.1]))

        # Add one base in the middle
        if np.random.rand()*2 < mutation_rate:
            rand_index = np.random.choice(range(len(member_portion)-1))
            member_portion = member_portion[:rand_index] + np.random.choice(['A', 'T', 'G', 'C'], p=[0.25, 0.25, 0.25, 0.25]) + member_portion[rand_index:]

        # Subtract one base in the middle
        if np.random.rand()*2 < mutation_rate:
            rand_index = np.random.choice(range(len(member_portion)-1))
            member_portion = member_portion[:rand_index] + member_portion[rand_index+1:]

        # Subtract one base from the end of the editable region
        if np.random.rand() < mutation_rate:
            member_portion = member_portion[:-1]

        # Add one base to the end of the editable region
        elif np.random.rand()*2 < mutation_rate:
            member_portion = member_portion[:-2]

    
    return member_portion


def create_ideal_structure(TARGET_DNA_LENGTH):
    """
    This function creates the theoretical ideal structure of the probe in dot
    bracket notation. It relies on the target DNA length to calculate how many
    bounded pairs (first and second given) and unbounded pairs the sequence
    should have.
    """
    first_given = "...." + "("*round((TARGET_DNA_LENGTH-36)*(1/3)+6)
    second_given = ")"*round((TARGET_DNA_LENGTH-36)*(1/3)+6)
    desired_structure = first_given + "."*(TARGET_DNA_LENGTH-(len(first_given)+len(second_given))) + second_given
    return desired_structure


def align_mismatch_penalties(dynamic_penalty, desired_structure, current_seq):
    desired_structure, current_seq = desired_structure + ".", current_seq + "."
    penalize_gaps, start_of_bracket, dp_pointer, dscs_pointer, tuple_list = False, False, 0, 0, []
    
    while dp_pointer < len(dynamic_penalty):
        if dynamic_penalty[dp_pointer] == "[":
            penalize_gaps = True
            start_of_bracket = False
            penalty = re.search(r"\(\d\)", dynamic_penalty[dp_pointer:]).group()[1:-1]
            dp_pointer += 1
        elif dynamic_penalty[dp_pointer] == "]":
            penalize_gaps = False
            dp_pointer += len(re.search(r"\(\d\)", dynamic_penalty[dp_pointer:]).group())+1
        elif dynamic_penalty[dp_pointer] == "+":
            start_of_bracket = True
            dp_pointer += 1
        elif desired_structure[dscs_pointer] == "-" and penalize_gaps and not start_of_bracket:
            tuple_list.append((desired_structure[dscs_pointer], current_seq[dscs_pointer], 0))
            dscs_pointer += 1
        elif desired_structure[dscs_pointer] != "-" and penalize_gaps and not start_of_bracket:
            start_of_bracket = True
        elif desired_structure[dscs_pointer] == "-" and penalize_gaps and start_of_bracket:
            tuple_list.append((desired_structure[dscs_pointer], current_seq[dscs_pointer], penalty))
            dscs_pointer += 1
        elif desired_structure[dscs_pointer] == "-" and not penalize_gaps:
            tuple_list.append((desired_structure[dscs_pointer], current_seq[dscs_pointer], 0))
            dscs_pointer += 1
        elif desired_structure[dscs_pointer] != "-" and desired_structure[dscs_pointer] != current_seq[dscs_pointer]:
            tuple_list.append((desired_structure[dscs_pointer], current_seq[dscs_pointer], dynamic_penalty[dp_pointer]))
            dscs_pointer += 1
            dp_pointer += 1
        elif desired_structure[dscs_pointer] == current_seq[dscs_pointer]:
            tuple_list.append((desired_structure[dscs_pointer], current_seq[dscs_pointer], 0))
            dscs_pointer += 1
            dp_pointer += 1
        elif dynamic_penalty[dp_pointer] == "-" and dscs_pointer - len(desired_structure)-1 != 0:
            tuple_list.append((desired_structure[dscs_pointer], current_seq[dscs_pointer], penalty))
            dscs_pointer += 1
            
    total_penalty = 0
    for i in tuple_list:
        total_penalty += int(i[-1])
    return total_penalty

def compute_segment_alignment_penalty(desired_structure_str,
                                      desired_align_str,
                                      current_align_str,
                                      segment_ranges):
    """
    Treat each fine-tuned structural segment as a 'meta-base' and penalize
    mismatches between desired and current structures inside those segments.

    - desired_structure_str: original desired dot-bracket (no gaps)
    - desired_align_str: aligned desired structure from aligner2 (with '-')
    - current_align_str: aligned current structure from aligner2 (with '-')
    - segment_ranges: list of (start, end) indices in ORIGINAL coordinates
                      (0-based, end exclusive), e.g. from find_fully_enclosed_segments.

    Returns:
        segment_penalty (non-negative int): number of mismatched positions inside
        the structural segments.
    """
    if not segment_ranges:
        return 0

    n = len(desired_structure_str)
    # Map each original desired index -> list of aligned indices
    orig_to_aligned = [[] for _ in range(n)]
    orig_idx = 0
    for ai, ch in enumerate(desired_align_str):
        if ch == "-":
            continue
        if orig_idx < n:
            orig_to_aligned[orig_idx].append(ai)
            orig_idx += 1
        else:
            break

    mismatches = 0
    for (start, end) in segment_ranges:
        # Safety clipping
        start = max(0, start)
        end   = min(n, end)
        for orig in range(start, end):
            for ai in orig_to_aligned[orig]:
                if ai >= len(current_align_str):
                    continue
                if desired_align_str[ai] != current_align_str[ai]:
                    mismatches += 1

    return mismatches



def score_hairpin_batch(hairpin_dict, DESIRED_STRUCTURE, DYNAMIC_PENALTY, generation,
                        penalize_off_targets=True):
    """
    Score how well UNAFold-predicted structures match the desired structure.

    Now includes an extra penalty that enforces structural subsection fidelity:
    fine-tuned structural segments (from find_fully_enclosed_segments) are
    treated as 'meta-bases', and mismatches inside these segments are added
    to the alignment score.
    """

    scores = []
    desired_structure = Seq(DESIRED_STRUCTURE)

    # Guard against empty DESIRED_STRUCTURE
    if len(desired_structure) == 0:
        raise ValueError("DESIRED_STRUCTURE is empty in score_hairpin_batch")

    desired_structure_str = str(desired_structure)

    # --- NEW: compute fine structural segments once (as in quick_preopt_sequence) ---
    fake_seq = "N" * len(desired_structure_str)
    fine_segments = find_fully_enclosed_segments(
        fake_seq,
        desired_structure_str,
        max_len=10,
        min_len=3
    )
    # fine_segments is list of (start, end) in original coordinates

    max_alignment_score = aligner1.align(desired_structure, desired_structure)[0].score

    for hairpin_list in hairpin_dict.values():
        # If hairpin_list itself is empty, treat as very bad
        if not hairpin_list:
            scores.append([1e9, 1])
            continue

        # Last structure is treated as the primary (MFE) structure
        current_seq = hairpin_list[-1]

        # If UNAFold produced no structure / dot-bracket, treat as very bad
        if not current_seq or len(current_seq) == 0:
            scores.append([1e9, 1])
            continue

        # === Global structural mismatch part ===
        score = max_alignment_score - aligner1.align(desired_structure, current_seq)[0].score

        alignment = aligner2.align(desired_structure, current_seq)[0]
        desired_structure_align = str(alignment[0])
        current_structure_align = str(alignment[1])

        # Existing dynamic penalties
        penalties = align_mismatch_penalties(DYNAMIC_PENALTY,
                                             desired_structure_align,
                                             current_structure_align)
        score += penalties * (generation / 2)

        # --- NEW: segment-level structural fidelity penalty ---
        # Treat each fine structural segment as a meta-base and penalize mismatches
        segment_penalty = compute_segment_alignment_penalty(
            desired_structure_str,
            desired_structure_align,
            current_structure_align,
            fine_segments,
        )

        # You can scale this if it feels too strong/weak:
        # e.g., score += segment_penalty * 0.5
        score += segment_penalty

        # Off-target count only if we want to penalize it
        if penalize_off_targets:
            other_structs = len(hairpin_list)
        else:
            other_structs = 1

        scores.append([score, other_structs])

    return scores



def parse_seq_constraint(seq_constraint_string):
    defined_bases = ["A", "C", "T", "G"]
    in_bracket = False
    current = ""
    seq_list = []
    seq_list_nested = []
    for i in seq_constraint_string:
        if i == "[":
            in_bracket = True
            if current != "":
                seq_list.append(current)
                current = ""
        elif i != "]" and in_bracket:
            if i == "|":
                seq_list_nested.append(current)
                current = ""
            else:
                current += i
        elif i == "]":
            seq_list_nested.append(current)
            seq_list_nested.append(1)
            seq_list.append(seq_list_nested)
            seq_list_nested, in_bracket, current = [], False, ""
        else:
            current += i
    if current != "":
        seq_list.append(current)

    return seq_list



def fitness_calc(sequence_length, melting_temp, alignment_score, other_structs, generation,
                 TARGET_DNA_LENGTH, TARGET_MELTING_TEMPERATURE,
                 use_tm_constraint=True,
                 penalize_off_targets=True):
    """
    Compute fitness for a candidate.

    - If use_tm_constraint is False, melting_temp may be None and the Tm term is skipped.
    - If penalize_off_targets is False, the off-target structures term is skipped entirely.
    """

    # Optional Tm penalty
    tm_term = 0
    if use_tm_constraint and TARGET_MELTING_TEMPERATURE is not None and melting_temp is not None:
        tm_term = abs(melting_temp - TARGET_MELTING_TEMPERATURE) * 2

    # Optional off-target structures penalty
    off_term = 0
    if penalize_off_targets and (other_structs is not None):
        off_term = abs(other_structs - 1) * (generation ** 1.5)

    if not penalize_off_targets and not use_tm_constraint:
        #alignment_score_multiplier = generation + 1
        alignment_score_multiplier = 1
    else:
        alignment_score_multiplier = 1

    if TARGET_DNA_LENGTH > 50:
        length_nonimportance = 4
    else:
        length_nonimportance = 1

    score = (
        abs(sequence_length - TARGET_DNA_LENGTH) / length_nonimportance
        + tm_term
        + alignment_score * alignment_score_multiplier
        + off_term
    )
    return score


def encode_sequences(population):
    # Step 1: Align sequences using Clustal Omega
    with tempfile.NamedTemporaryFile(delete=False, mode="w") as input_file, tempfile.NamedTemporaryFile(delete=False, mode="w") as output_file:
        # Write sequences to the input FASTA file
        for i, (seq, _) in enumerate(population):
            input_file.write(f">seq{i}\n{seq}\n")
        input_file.flush()
        
        # Call Clustal Omega
        clustalomega_cline = ClustalOmegaCommandline(
            cmd="clustalo",
            infile=input_file.name,
            outfile=output_file.name,
            force=True,
            auto=True,
            verbose=False
        )
        clustalomega_cline()
    
    # Read the aligned sequences
    with open(output_file.name) as f:
        aligned_records = list(SeqIO.parse(f, "fasta"))
    
    # Clean up temporary files
    os.unlink(input_file.name)
    os.unlink(output_file.name)
    
    # Convert aligned sequences to MultipleSeqAlignment object
    alignment = MultipleSeqAlignment(
        [SeqRecord(Seq(str(record.seq)), id=record.id) for record in aligned_records]
    )

    # Step 2: Calculate conservation scores for each aligned position
    alignment_length = alignment.get_alignment_length()
    conservation_scores = []

    for i in range(alignment_length):
        column = [str(record.seq[i]) for record in alignment]
        # Exclude gaps ("-") from the column
        filtered_column = [base for base in column if base != "-"]
        if filtered_column:  # Only calculate if there are non-gap bases
            most_common = Counter(filtered_column).most_common(1)[0]
            conservation_scores.append((i, most_common[1]))
        else:  # If all bases are gaps, conservation score is 0
            conservation_scores.append((i, 0))

    
    # Sort positions by their conservation score
    sorted_positions = sorted(conservation_scores, key=lambda x: x[1], reverse=True)
    sorted_indices = [pos[0] for pos in sorted_positions]
    
    # Step 3: Rearrange bases in each sequence according to the most conserved positions
    rearranged_sequences = []
    for record in alignment:
        rearranged_seq = ''.join([record.seq[i] for i in sorted_indices])
        rearranged_sequences.append(rearranged_seq)
    
    # Step 4: Encoding and Decoding
    # Encode by storing the original indices and the rearranged sequences
    ends = [seq[1] for seq in population]
    encoded_sequences = [''.join([seq[i] for i in sorted_indices]) for seq in rearranged_sequences]
    
    return list(zip(encoded_sequences, ends)), sorted_indices, alignment_length


def decode_sequences(child, sorted_indices, max_len):
    decoded_child = ''.join(child[sorted_indices.index(i)] for i in range(max_len))
    decoded_child = decoded_child.replace("-", "")
    return decoded_child

def mut_rate(generation, NUM_GENERATIONS):
    return 0.15

def is_folder_safe_to_delete(folder_name):
    return not any(file.endswith(('.fasta', '.dG', '.ct')) for file in os.listdir(folder_name))

def clean_folder(folder_name, retries=5, delay=1):
    for file_name in os.listdir(folder_name):
        if file_name.endswith(('.ann', '.ct', '.plot', '.fasta.dG', '.fasta', '.fasta.run')):
            file_path = os.path.join(folder_name, file_name)
            for attempt in range(retries):
                try:
                    os.remove(file_path)
                    break  # Exit loop if successful
                except PermissionError:
                    if attempt < retries - 1:
                        time.sleep(delay)  # Wait before retrying
                    else:
                        print(f"Warning: Could not delete {file_path} after {retries} attempts.")

def delete_old_unafold_folders():
    """
    Deletes all run folders under UNAFOLD_WORK_DIR before the GA starts.
    This prevents conflicts from leftover data in previous runs.
    """
    if not UNAFOLD_WORK_DIR.exists():
        return
    for folder in UNAFOLD_WORK_DIR.glob("run_*"):
        if folder.is_dir():
            try:
                shutil.rmtree(folder)
                print(f"Deleted old folder: {folder}")
            except Exception as e:
                print(f"Warning: Could not delete {folder} - {e}")


def run_single_GA(i, miRNA, DESIRED_STRUCTURE, SEQ_CONSTRAINT, DYNAMIC_PENALTY, SALT_CONC, MAG_CONC,
                  TARGET_MELTING_TEMPERATURE, NUM_GENERATIONS, POPULATION_SIZE,
                  use_tm_constraint=True,
                  penalize_off_targets=True,
                  cancel_cb=None,
                  starter_override=None):

    
    try:
        folder_number = i + 1
        folder_name = str(UNAFOLD_WORK_DIR / f"run_{folder_number}")
        
        print(f"Running GA in folder {folder_name} (Run {i+1})...")

        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        
        seq_constraint_list = parse_seq_constraint(SEQ_CONSTRAINT)
        TARGET_DNA_LENGTH = len(DESIRED_STRUCTURE)

        # Build STARTER, optionally overriding with a pre-optimized token list
        if starter_override is not None:
            STARTER = starter_override
        else:
            STARTER = ["A" * len(tok) if isinstance(tok, str) else tok
                       for tok in seq_constraint_list]

        # Initialize the population with the requested size
        population = initialize_dna_strand(STARTER, 2000)

        
        for generation in range(NUM_GENERATIONS):
            # Cooperative cancellation
            if cancel_cb is not None and cancel_cb():
                print(f"Cancellation requested in run_single_GA (folder {folder_name}). Exiting early.")
                return None  # sentinel for "cancelled"

            actual_population = []
            for member in population:
                member_complete = ""
                for element in member:
                    if isinstance(element, list):
                        element_index = element[-1]
                        member_complete += element[element_index]
                    else:
                        member_complete += element
                actual_population.append(member_complete)

            fasta_file = os.path.join(folder_name, f'sequences_{folder_number}.fasta')
            make_fasta(actual_population, fasta_file)

            # UNAFold run
            idt_hairpin_batch(
                SALT_CONC,
                MAG_CONC,
                folder_name,
                fasta_file=fasta_file,
                use_tm_constraint=use_tm_constraint,
                penalize_off_targets=penalize_off_targets
            )

            hairpin_dict = parse_ct_file(os.path.join(folder_name, f'sequences_{folder_number}.fasta.25.ct'))

            # Structural scores (with or without off-target counting)
            hairpin_scores = score_hairpin_batch(
                hairpin_dict,
                DESIRED_STRUCTURE,
                DYNAMIC_PENALTY,
                generation,
                penalize_off_targets=penalize_off_targets
            )

            # Only calculate melting temperatures if we are actually using them
            if use_tm_constraint:
                temps = parse_dG_file_and_calc_MT(
                    os.path.join(folder_name, f'sequences_{folder_number}.fasta.dG')
                )
            else:
                temps = [None] * len(actual_population)

            pop_temps_and_scores = zip(actual_population, temps, hairpin_scores)

            fitness_scores = [
                fitness_calc(
                    len(dna[0]),        # sequence length
                    dna[1],             # melting_temp (may be None)
                    dna[2][0],          # alignment_score
                    dna[2][1],          # other_structs (1 if penalize_off_targets=False)
                    generation,
                    TARGET_DNA_LENGTH,
                    TARGET_MELTING_TEMPERATURE,
                    use_tm_constraint=use_tm_constraint,
                    penalize_off_targets=penalize_off_targets,
                )
                for dna in pop_temps_and_scores
            ]

            best_index = np.argmin(fitness_scores)
            best_temp = temps[best_index]  # may be None
            best_fitness = fitness_scores[best_index]
            best_dna_full = actual_population[best_index]

            # ===== NEW: early stop if we hit perfect (zero) fitness =====
            # Using <= 0 for robustness in case of tiny numerical noise.
            if best_fitness <= 0:
                hairpin_dict_best_seq = list(hairpin_dict.values())[best_index]
                structure = hairpin_dict_best_seq[-1]

                if penalize_off_targets:
                    other_structs = len(hairpin_dict_best_seq)
                else:
                    other_structs = 1

                if best_temp is not None:
                    temp_str = f"{best_temp:.2f}"
                else:
                    temp_str = "N/A"

                print(
                    f"(UNAFold: run {folder_number}) Generation {generation + 1}: "
                    f"Best DNA = {''.join(best_dna_full)}, "
                    f"Melting Temperature = {temp_str}, "
                    f"Fitness = {best_fitness:.2f}, Length = {len(best_dna_full)}"
                )
                print("Minimum Free Energy (MFE) Structure:", "\n", structure)
                print("Number of other structures: ", other_structs, "\n")

                # Immediately return this perfect solution for this run
                return best_fitness, best_dna_full, structure, best_temp, True, other_structs
            # ===== END NEW BLOCK =====

            if (generation + 1) % 10 == 0 or generation + 1 == NUM_GENERATIONS:
                hairpin_dict_best_seq = list(hairpin_dict.values())[best_index]
                structure = hairpin_dict_best_seq[-1]

                if penalize_off_targets:
                    other_structs = len(hairpin_dict_best_seq)
                else:
                    other_structs = 1

                if best_temp is not None:
                    temp_str = f"{best_temp:.2f}"
                else:
                    temp_str = "N/A"

                print(
                    f"(UNAFold: run {folder_number}) Generation {generation + 1}: "
                    f"Best DNA = {''.join(best_dna_full)}, "
                    f"Melting Temperature = {temp_str}, "
                    f"Fitness = {best_fitness:.2f}, Length = {len(best_dna_full)}"
                )
                print("Minimum Free Energy (MFE) Structure:", "\n", structure)
                print("Number of other structures: ", other_structs, "\n")

                if generation + 1 == NUM_GENERATIONS:
                    return best_fitness, best_dna_full, structure, best_temp, True, other_structs


            # parent selection, crossover, mutation, etc. unchanged...
            parents_indices = np.argsort(fitness_scores)[:int(POPULATION_SIZE * 0.5)]
            parents = [copy.deepcopy(population[i]) for i in parents_indices]

            mutating = POPULATION_SIZE - int(((3/2)*len(parents)) - (len(parents)/2 * (generation/NUM_GENERATIONS)))
            crossing_over_number = POPULATION_SIZE - len(parents) - mutating

            crossing_over_children = []
            parents_rearranged = parents

            for _ in range(crossing_over_number):
                parent1_index, parent2_index = np.random.choice(len(parents_rearranged), size=2, replace=False)
                parent1, parent2 = parents_rearranged[parent1_index], parents_rearranged[parent2_index]

                child = []
                for gene_idx in range(len(parent1)):
                    # Alt-blocks (lists) are left unchanged
                    if isinstance(parent1[gene_idx], list):
                        child.append(parent1[gene_idx])
                        continue

                    seq1 = parent1[gene_idx]
                    seq2 = parent2[gene_idx]

                    # If either sequence is too short, skip crossover and just copy seq1
                    if seq1 is None:
                        seq1 = ""
                    if seq2 is None:
                        seq2 = ""

                    min_len = min(len(seq1), len(seq2))

                    # Need at least 3 bases to choose an internal crossover point (1 .. min_len-2)
                    if min_len <= 2:
                        child_seq = seq1  # fallback: no crossover
                    else:
                        high = min_len - 1  # exclusive upper bound
                        # Now low=1, high>=2 → safe: low < high
                        crossover_point = np.random.randint(1, high)
                        child_seq = seq1[:crossover_point] + seq2[crossover_point:]

                    child.append(child_seq)

                # Mutate the resulting child
                child_mut = []
                for element in child:
                    child_mut.append(mutate(element, mut_rate(generation, NUM_GENERATIONS)))
                crossing_over_children.append(child_mut)


            mutating_children = []
            for parent in parents[:mutating]:
                child_mut = []
                for element in parent:
                    child_mut.append(mutate(element, mut_rate(generation, NUM_GENERATIONS)))
                mutating_children.append(child_mut)

            population = copy.deepcopy(parents) + copy.deepcopy(crossing_over_children) + copy.deepcopy(mutating_children)

            # Safely clean up ct files if they exist
            for j in range(24, 28):
                ct_path = os.path.join(folder_name, f'sequences_{folder_number}.fasta.{j}.ct')
                if os.path.exists(ct_path):
                    try:
                        os.remove(ct_path)
                    except PermissionError:
                        print(f"Warning: could not delete {ct_path} (permission error)")


    except Exception as e:
        print(f"Error in run_single_GA {folder_number}: {e}.")

        # If this is the Biopython "sequence has zero length" error (or similar),
        # it's not transient – don't recurse forever.
        msg = str(e)
        if "sequence has zero length" in msg or "zero length" in msg:
            # Return None so callers (quick_preopt_sequence / GA) can skip this run
            return None

        # For other errors (e.g. file-lock / transient IO), keep the old cleanup+retry
        if is_folder_safe_to_delete(folder_name):
            clean_folder(folder_name)

        return run_single_GA(
            i, miRNA, DESIRED_STRUCTURE, SEQ_CONSTRAINT, DYNAMIC_PENALTY,
            SALT_CONC, MAG_CONC, TARGET_MELTING_TEMPERATURE,
            NUM_GENERATIONS, POPULATION_SIZE,
            use_tm_constraint=use_tm_constraint,
            penalize_off_targets=penalize_off_targets,
            cancel_cb=cancel_cb,
            starter_override=starter_override
        )


def quick_preopt_sequence(miRNA,
                          desired_structure,
                          SALT_CONC,
                          MAG_CONC,
                          generations,
                          population_size,
                          cancel_cb=None,
                          # ---------- USER-TUNABLE PARAMETERS ----------
                          FINE_GENERATIONS=None,
                          FINE_RUNS_PER_SEGMENT=1,
                          COARSE_GENERATIONS=None,
                          COARSE_RUNS_PER_SEGMENT=2):
    """
    Two-level warm-start using ensemble_helper.find_fully_enclosed_segments.

    Fine-level segmentation:
        max_len = 10, min_len = 3

    Coarse-level segmentation:
        max_len = 50, min_len = 30

    This function NEVER calls GA(), only run_single_GA().
    """

    # ------------------------------------------------------------
    # Default generations if not overridden
    # ------------------------------------------------------------
    if FINE_GENERATIONS is None:
        FINE_GENERATIONS = 30

    if COARSE_GENERATIONS is None:
        COARSE_GENERATIONS = 50

    # Keep preopt population smaller than main GA
    quick_pop = min(100, population_size)

    n = len(desired_structure)
    if n == 0:
        return None

    fake_seq = "N" * n

    # ---------- Raw segments ----------
    raw_fine_segments = find_fully_enclosed_segments(
        fake_seq,
        desired_structure,
        max_len=10,
        min_len=6
    )

    raw_coarse_segments = find_fully_enclosed_segments(
        fake_seq,
        desired_structure,
        max_len=50,
        min_len=10
    )

    # ---------- Sanitize segments: enforce 0 <= start < end <= n ----------
    def _sanitize_segments(seg_list, n):
        cleaned = []
        for (s, e) in seg_list:
            if s is None or e is None:
                continue
            if not isinstance(s, int) or not isinstance(e, int):
                continue
            if s < 0 or e < 0:
                continue
            if s >= n or e > n:
                continue
            if e <= s:
                continue
            cleaned.append((s, e))
        return cleaned

    fine_segments = _sanitize_segments(raw_fine_segments, n)
    coarse_segments = _sanitize_segments(raw_coarse_segments, n)

    if not fine_segments and not coarse_segments:
        # Nothing reasonable to pre-optimize
        return None

    # Start with neutral sequence
    pre_seq = list("A" * n)

    # ==========================================================
    # 1) FINE SEGMENT GA PASS
    # ==========================================================
    for seg_idx, (start, end) in enumerate(fine_segments):
        if cancel_cb and cancel_cb():
            return None

        seg_len = end - start
        if seg_len <= 0:
            continue  # extra safety

        sub_struct = desired_structure[start:end]
        if len(sub_struct) == 0:
            continue  # guard against weird slicing

        sub_seq_constraint  = "N" * seg_len
        sub_dynamic_penalty = "0" * seg_len

        # Avoid folder conflicts
        base_run_index = 1000 + seg_idx * 10
        indices = [base_run_index + r for r in range(FINE_RUNS_PER_SEGMENT)]

        with concurrent.futures.ThreadPoolExecutor(
                max_workers=FINE_RUNS_PER_SEGMENT) as executor:

            futures = [
                executor.submit(
                    run_single_GA,
                    idx,
                    miRNA,
                    sub_struct,
                    sub_seq_constraint,
                    sub_dynamic_penalty,
                    SALT_CONC,
                    MAG_CONC,
                    0.0,                # Tm ignored
                    FINE_GENERATIONS,   # fine-level generations
                    quick_pop,          # warm-start pop size
                    False,              # use_tm_constraint
                    False,              # penalize_off_targets
                    cancel_cb,
                    None                # starter_override
                )
                for idx in indices
            ]

            results = [f.result() for f in futures]

        # Select best result
        best_result = None
        for r in results:
            if r is None:
                continue
            if best_result is None or r[0] < best_result[0]:
                best_result = r

        if best_result:
            best_fitness, best_seq, _, _, _, _ = best_result
            best_seq = best_seq[:seg_len].ljust(seg_len, "A")
            pre_seq[start:end] = list(best_seq)

    # ==========================================================
    # 2) COARSE SEGMENT GA PASS (seeded with fine-pass output)
    # ==========================================================
    for seg_idx, (start, end) in enumerate(coarse_segments):
        if cancel_cb and cancel_cb():
            return None

        seg_len = end - start
        if seg_len <= 0:
            continue

        sub_struct = desired_structure[start:end]
        if len(sub_struct) == 0:
            continue

        sub_seq_constraint  = "N" * seg_len
        sub_dynamic_penalty = "0" * seg_len

        # Seed this region with result of fine-pass
        starter_override = ["".join(pre_seq[start:end])]

        base_run_index = 2000 + seg_idx * 10
        indices = [base_run_index + r for r in range(COARSE_RUNS_PER_SEGMENT)]

        with concurrent.futures.ThreadPoolExecutor(
                max_workers=COARSE_RUNS_PER_SEGMENT) as executor:

            futures = [
                executor.submit(
                    run_single_GA,
                    idx,
                    miRNA,
                    sub_struct,
                    sub_seq_constraint,
                    sub_dynamic_penalty,
                    SALT_CONC,
                    MAG_CONC,
                    0.0,
                    COARSE_GENERATIONS,   # coarse-level generations
                    quick_pop,
                    False,
                    False,
                    cancel_cb,
                    starter_override
                )
                for idx in indices
            ]

            results = [f.result() for f in futures]

        best_result = None
        for r in results:
            if r is None:
                continue
            if best_result is None or r[0] < best_result[0]:
                best_result = r

        if best_result:
            best_fitness, best_seq, _, _, _, _ = best_result
            best_seq = best_seq[:seg_len].ljust(seg_len, "A")
            pre_seq[start:end] = list(best_seq)

    return "".join(pre_seq)



# --- GA signature: add generate_placeholder flag (default True to preserve old behavior)
def GA(miRNA,
       target_melting_temp=49,
       desired_structure="....((((((....................))))))",
       seq_constraint=None,
       struct_constraints=None,
       population_size=250,
       generations=100,
       runs=10,
       salt_conc=0.30,
       mag_conc=0.0,
       automatic_ph_mt=False,
       desired_ph_mt=55,
       generate_placeholder=True,
       batch_size=1,
       use_tm_constraint=True,      # existing
       penalize_off_targets=True,   # NEW
       cancel_cb=None):
    """
    Top-level GA driver.
    - penalize_off_targets: if False, GA does NOT penalize extra suboptimal structures.
    """

    delete_old_unafold_folders()

    TARGET_MELTING_TEMPERATURE = target_melting_temp
    DESIRED_STRUCTURE = desired_structure
    TARGET_DNA_LENGTH = len(desired_structure)
    if struct_constraints is None:
        struct_constraints = "0" * TARGET_DNA_LENGTH
    if seq_constraint is None:
        seq_constraint = "N" * TARGET_DNA_LENGTH
    POPULATION_SIZE = population_size
    NUM_GENERATIONS = generations

    SALT_CONC = salt_conc
    MAG_CONC  = mag_conc

    # -----------------------------------------
    # Editable-region length and warm-start GA
    # -----------------------------------------
    seq_constraint_list = parse_seq_constraint(seq_constraint)

    # Treat string tokens as editable blocks; list tokens are section-alts we leave alone.
    editable_len = sum(len(tok) for tok in seq_constraint_list if isinstance(tok, str))

    starter_override = None

    if editable_len > 70:
        # Use ensemble_helper-based segment GA to get a warm-start sequence
        pre_seq = quick_preopt_sequence(
            miRNA,
            DESIRED_STRUCTURE,
            SALT_CONC,
            MAG_CONC,
            generations=NUM_GENERATIONS,
            population_size=POPULATION_SIZE,
            cancel_cb=cancel_cb,
        )

        if pre_seq is not None and len(pre_seq) == TARGET_DNA_LENGTH:
            # Convert the pre-optimized flat sequence into a token list
            # aligned with parse_seq_constraint output.
            starter_override = []
            idx = 0
            for tok in seq_constraint_list:
                if isinstance(tok, str):
                    L = len(tok)
                    starter_override.append(pre_seq[idx:idx + L])
                    idx += L
                else:
                    # Section-alts: keep them as-is, but advance index by their length.
                    if tok and isinstance(tok[0], str):
                        L = len(tok[0])
                    else:
                        L = 0
                    starter_override.append(tok)
                    idx += L


    def run_in_batches(runs_to_do, cancel_cb=None):
        fitnesses = []
        failed_runs = []
        cancelled = False

        bs = max(1, int(batch_size))

        for batch_start in range(0, len(runs_to_do), bs):
            if cancel_cb is not None and cancel_cb():
                cancelled = True
                break

            batch = runs_to_do[batch_start:batch_start + bs]
            with concurrent.futures.ThreadPoolExecutor(max_workers=bs) as executor:
                results = list(executor.map(
                    run_single_GA,
                    batch,
                    [miRNA] * len(batch),
                    [desired_structure] * len(batch),
                    [seq_constraint] * len(batch),
                    [struct_constraints] * len(batch),
                    [SALT_CONC] * len(batch),
                    [MAG_CONC] * len(batch),
                    [TARGET_MELTING_TEMPERATURE] * len(batch),
                    [NUM_GENERATIONS] * len(batch),
                    [POPULATION_SIZE] * len(batch),
                    [use_tm_constraint] * len(batch),
                    [penalize_off_targets] * len(batch),
                    [cancel_cb] * len(batch),
                    [starter_override] * len(batch),
                ))


            for i, result in enumerate(results):
                if cancel_cb is not None and cancel_cb():
                    cancelled = True
                    break

                if result is None:
                    continue

                fitness, sequence, mfe_structure, melt_temp, save, alt_structs = result
                if save:
                    fitnesses.append([fitness, sequence, mfe_structure, melt_temp, alt_structs])
                else:
                    failed_runs.append(batch[i])

            if cancelled:
                break

        return fitnesses, failed_runs, cancelled

    all_runs = list(range(runs))
    fitnesses, failed_runs, cancelled = run_in_batches(all_runs, cancel_cb=cancel_cb)

    while failed_runs and not (cancel_cb is not None and cancel_cb()):
        retry_fitnesses, failed_runs, cancelled2 = run_in_batches(failed_runs, cancel_cb=cancel_cb)
        fitnesses.extend(retry_fitnesses)
        if cancelled2:
            break

    fitnesses = sorted(fitnesses, key=lambda x: x[0])

    reports = []
    for rank, fit in enumerate(fitnesses):
        seq   = fit[1]
        tm    = fit[3]          # may be None
        fitv  = float(fit[0])
        struct = fit[2]
        alt_structs = int(fit[4])

        report = {
            "rank": rank,
            "sequence": seq,
            "structure": struct,
            "melting_temp": float(tm) if tm is not None else None,
            "fitness": fitv,
            "length": len(seq),
            "alt_structures": alt_structs,
        }

        if generate_placeholder and miRNA:
            probe_tm_for_placeholder = tm if tm is not None else 0.0
            placeholder_seq, placeholder_melt, mirna_placeholder_melt = findPlaceholderSequences(
                miRNA, seq, probe_tm_for_placeholder, automatic_ph_mt, desired_ph_mt
            )
            report.update({
                "placeholder": placeholder_seq,
                "placeholder_tm": float(placeholder_melt),
                "mirna_tm": float(mirna_placeholder_melt),
            })

        reports.append(report)

        tm_display = f"{tm:.2f}" if tm is not None else "N/A"
        print(f"#{rank}: DNA: {seq}, Melting Temperature: {tm_display}, Fitness: {round(fitv, 2)}")
        print(f"MFE Structure:\n{struct}")
        if "placeholder" in report:
            print(
                f"placeholder strand: {report['placeholder']}\n"
                f"placeholder Melting Temp: {report['placeholder_tm']}\n"
                f"miRNA Melting Temp: {report['mirna_tm']}\n"
            )

    return reports





