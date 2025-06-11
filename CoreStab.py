import subprocess
import sys
import re
import os
import shutil
import argparse
from Bio.PDB import PDBParser, NeighborSearch

three_to_one = {
    'ALA': 'A', 'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'PHE': 'F',
    'TRP': 'W', 'PRO': 'P', 'GLY': 'G', 'SER': 'S', 'THR': 'T', 'CYS': 'C',
    'TYR': 'Y', 'ASN': 'N', 'GLN': 'Q', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K',
    'ARG': 'R', 'HIS': 'H'
}

hydrophobic_residues = ['A', 'V', 'L', 'I', 'F', 'W', 'M']

size_rank_volumeC = {
    'A': 1, 'V': 2, 'I': 3, 'L': 4, 'M': 5, 'F': 6, 'W': 7
}


def clean_pdb(pdb_file):
    cleaned_lines = []
    current_chain = 'A'

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith(('ATOM', 'HETATM')):
                res_name = line[17:20].strip()
                if line.startswith('HETATM') and res_name == 'HOH':
                    continue
                line = line[:21] + current_chain + line[22:]
                cleaned_lines.append(line)

    cleaned_pdb_file = pdb_file.replace('.pdb', '_cleaned.pdb')
    with open(cleaned_pdb_file, 'w') as file:
        file.writelines(cleaned_lines)

    return cleaned_pdb_file


def run_freesasa(pdb_file):
    result = subprocess.run(['freesasa', pdb_file, '--format=seq'], capture_output=True, text=True)
    return result.stdout


def parse_freesasa_output(output, threshold=1.0):
    core_residues = []
    for line in output.split('\n'):
        match = re.match(r'SEQ\s*(\w?)\s*(\d+)\s*(\w+)\s*:\s*([\d\.]+)', line)
        if match:
            chain, res_num, res_name, sasa = match.groups()
            chain = 'A'
            if float(sasa) < threshold:
                core_residues.append((chain, res_num, res_name))
    return core_residues


def get_excluded_residues(pdb_file, functional_res_list, distance_threshold=4.5, shell_level=1):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    model = next(structure.get_models())
    chain = model.child_list[0]  # Only chain in monomeric PDB

    atoms = list(chain.get_atoms())
    ns = NeighborSearch(atoms)

    excluded = set()
    res_id_map = {res.id[1]: res for res in chain if res.id[0] == ' '}

    first_shell = set()
    for res_num_str in functional_res_list:
        res_num = int(res_num_str)
        if res_num not in res_id_map:
            continue

        res = res_id_map[res_num]
        excluded.add(('A', str(res_num)))

        if 'CA' not in res:
            continue

        neighbors = ns.search(res['CA'].coord, distance_threshold)
        for neighbor in neighbors:
            n_res = neighbor.get_parent()
            if n_res.id[0] == ' ':
                first_shell.add(n_res.id[1])
                excluded.add(('A', str(n_res.id[1])))

    if shell_level == 2:
        second_shell = set()
        for res_num in first_shell:
            if res_num not in res_id_map:
                continue
            res = res_id_map[res_num]
            if 'CA' not in res:
                continue
            neighbors = ns.search(res['CA'].coord, distance_threshold)
            for neighbor in neighbors:
                n_res = neighbor.get_parent()
                if n_res.id[0] == ' ':
                    excluded.add(('A', str(n_res.id[1])))

    return excluded


def generate_mutation_list(core_residues, if_hydrophobic=False, excluded_residues=set()):
    mutations = []
    mutation_dict = {}

    for chain, res_num, res_name in core_residues:
        if (chain, res_num) in excluded_residues:
            continue

        original_res = three_to_one.get(res_name, 'X')
        if original_res not in hydrophobic_residues and if_hydrophobic:
            continue

        for target_res in hydrophobic_residues:
            if target_res == original_res:
                continue
            if if_hydrophobic:
                if size_rank_volumeC.get(target_res, 0) <= size_rank_volumeC.get(original_res, 0):
                    continue
            mutation = f'{original_res}{chain}{res_num}{target_res};'
            mutations.append(mutation)
            mutation_dict.setdefault(f'{chain}-{res_name}{res_num}', []).append(target_res)

    with open('individual_list.txt', 'w') as f:
        f.write('\n'.join(mutations))

    return mutations, mutation_dict


def run_foldx(pdb_file):
    with open('foldx_log.txt', 'w') as log_file:
        result = subprocess.run(
            [
                'foldx',
                '--command=BuildModel',
                '--pdb=' + pdb_file,
                '--mutant-file=individual_list.txt',
                '--numberOfRuns=5'
            ],
            stdout=log_file,
            stderr=log_file
        )

    if result.returncode != 0:
        raise RuntimeError(
            f"\n❌ FoldX failed with return code {result.returncode}. "
            f"Check '{os.path.abspath('foldx_log.txt')}' for errors.\n"
            f"Current working directory: {os.getcwd()}"
        )


def parse_foldx_output(pdb_file, mutations):
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    output_file = f'Average_{base_name}.fxout'

    if not os.path.exists(output_file):
        cleaned_base_name = base_name.replace('_cleaned', '')
        output_file = f'Average_{cleaned_base_name}.fxout'

    if not os.path.exists(output_file):
        raise FileNotFoundError(
            f"\n FoldX output file not found. Tried: '{output_file}'.\n"
            f"Make sure FoldX completed successfully and generated results.\n"
            f"Expected path: {os.path.abspath(output_file)}"
        )

    negative_ddg_mutations = []
    with open(output_file, 'r', encoding='utf-8', errors='ignore') as file:
        lines = file.readlines()

    mutation_map = {}
    for i, line in enumerate(lines):
        if line.startswith(f'{base_name}_'):
            parts = line.strip().split()
            if len(parts) > 2:
                try:
                    ddg = float(parts[2])
                    mutation_map[i] = ddg
                except ValueError:
                    continue

    mutation_index = 0
    for i in sorted(mutation_map.keys()):
        if mutation_index < len(mutations):
            ddg = mutation_map[i]
            if ddg < 0:
                negative_ddg_mutations.append((mutations[mutation_index], ddg))
            mutation_index += 1

    return negative_ddg_mutations


def create_run_directory(pdb_file):
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    dir_name = base_name
    counter = 1

    while os.path.exists(dir_name):
        dir_name = f"{base_name}_{counter}"
        counter += 1

    try:
        os.makedirs(dir_name, exist_ok=True)
    except Exception as e:
        raise RuntimeError(f"Failed to create directory '{dir_name}': {e}")

    return dir_name


def main():
    parser = argparse.ArgumentParser(description="Mutational scanning pipeline using FreeSASA and FoldX.")
    parser.add_argument('pdb_file', help='Path to the PDB file')
    parser.add_argument('--ifHydrophobic', action='store_true', help='Restrict mutations to more hydrophobic residues only')
    parser.add_argument('--sasaThreshold', type=float, default=1.0, help='SASA threshold to define core residues (default=1.0)')
    parser.add_argument('--functionalRes', nargs='+', help='Functionally important residue numbers (e.g., 45 128 193)')
    parser.add_argument('--contactShell', type=int, choices=[1, 2], default=1, help='Contact shell to exclude: 1 (default) or 2')
    parser.add_argument('--contactCutoff', type=float, default=4.5, help='Distance cutoff in Å for defining contacts (default=4.5)')
    parser.add_argument('--verbose', action='store_true', help='Print verbose logs')

    args = parser.parse_args()

    pdb_file = args.pdb_file
    if_hydrophobic = args.ifHydrophobic
    sasa_threshold = args.sasaThreshold
    verbose = args.verbose

    run_dir = create_run_directory(pdb_file)
    original_dir = os.getcwd()
    os.chdir(run_dir)

    shutil.copy(os.path.join(original_dir, pdb_file), '.')

    if verbose:
        print(f"📁 Created run directory: {run_dir}")

    cleaned_pdb_file = clean_pdb(pdb_file)
    if verbose:
        print(f"🧽 Cleaned PDB file: {cleaned_pdb_file}")

    sasa_output = run_freesasa(cleaned_pdb_file)
    if verbose:
        print("📊 FreeSASA output:\n", sasa_output)

    core_residues = parse_freesasa_output(sasa_output, threshold=sasa_threshold)

    excluded_residues = set()
    if args.functionalRes:
        excluded_residues = get_excluded_residues(
            cleaned_pdb_file,
            args.functionalRes,
            distance_threshold=args.contactCutoff,
            shell_level=args.contactShell
        )
        if verbose:
            print("\n⛔ Excluding functional residues and their contacts from mutation:")
            for chain, res_num in sorted(excluded_residues, key=lambda x: int(x[1])):
                print(f"{chain}-{res_num}")

    if verbose:
        print("\n📊 Core residues identified (buried):")
        for chain, res_num, res_name in core_residues:
            print(f"{chain}-{res_name}{res_num}")

    mutations, mutation_dict = generate_mutation_list(core_residues, if_hydrophobic=if_hydrophobic, excluded_residues=excluded_residues)

    if verbose:
        print("\n📊 Mutants created:")
        for key, value in mutation_dict.items():
            print(f"{key} --> {''.join(value)}")

    run_foldx(cleaned_pdb_file)
    negative_ddg_mutations = parse_foldx_output(cleaned_pdb_file, mutations)

    print("\n📊 Residues selected for mutation from PDB:")
    for chain, res_num, res_name in core_residues:
        if (chain, res_num) not in excluded_residues:
            print(f"{chain}-{res_name}{res_num}")

    print("\n📊 Mutations with negative ΔΔG values:")
    for mutation, ddg in sorted(negative_ddg_mutations, key=lambda x: x[1]):
        print(f'{mutation}: {ddg:.3f}')

    os.chdir(original_dir)


if __name__ == '__main__':
    main()

