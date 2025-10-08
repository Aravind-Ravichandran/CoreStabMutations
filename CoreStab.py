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

ROSETTA_EXE = "~/Downloads/softwares/rosetta.source.release-371/main/source/bin/cartesian_ddg.linuxgccrelease"


def clean_pdb(pdb_file, keep_ions=False):
    """
    Clean PDB file:
    - Keeps only ATOM records (and optionally ions if keep_ions=True).
    - Removes water (HOH) and other hetero residues.
    - Forces all residues into chain A.
    - Renumbers residues sequentially starting from 1.
    Returns path to cleaned PDB.
    """
    cleaned_lines = []
    current_chain = 'A'
    residue_map = {}  # old_resnum -> new_resnum
    new_resnum = 1
    prev_old_res = None

    with open(pdb_file, 'r') as file:
        for line in file:
            if not line.startswith("ATOM") and not line.startswith("HETATM"):
                continue

            res_name = line[17:20].strip()

            # skip water
            if res_name == "HOH":
                continue

            # skip all HETATM unless ion retention requested
            if line.startswith("HETATM") and not keep_ions:
                continue

            # extract old residue number + insertion code
            old_resnum = line[22:27].strip()  # residue number + insertion

            if old_resnum != prev_old_res:
                residue_map[old_resnum] = new_resnum
                prev_old_res = old_resnum
                new_resnum += 1

            # replace chain with A
            line = line[:21] + current_chain + line[22:]

            # replace residue number with renumbered value
            renum_str = str(residue_map[old_resnum]).rjust(4)
            line = line[:22] + renum_str + line[26:]

            cleaned_lines.append(line)

    cleaned_pdb_file = pdb_file.replace(".pdb", "_cleaned.pdb")
    with open(cleaned_pdb_file, "w") as file:
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
    chain = model.child_list[0]
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
        raise RuntimeError("❌ FoldX failed. Check foldx_log.txt.")


def parse_foldx_output(pdb_file, mutations):
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    output_file = f'Average_{base_name}.fxout'
    if not os.path.exists(output_file):
        cleaned_base_name = base_name.replace('_cleaned', '')
        output_file = f'Average_{cleaned_base_name}.fxout'
    if not os.path.exists(output_file):
        raise FileNotFoundError(f"\n FoldX output file not found. Expected: '{output_file}'\n")
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


# ---------------- Rosetta Integration ---------------- #

def generate_rosetta_mutfile(mutations, pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    model = next(structure.get_models())
    chain = model.child_list[0]
    res_id_map = {res.id[1]: res for res in chain if res.id[0] == ' '}
    with open("mut_file.txt", "w") as f:
        f.write(f"total {len(mutations)}\n")
        for mut in mutations:
            mut = mut.strip(";")
            wt = mut[0]
            chain_id = mut[1]
            resnum = mut[2:-1]
            mutaa = mut[-1]
            resnum_int = int(resnum)
            if resnum_int not in res_id_map:
                raise ValueError(f"Residue {resnum} not found in PDB for mutation {mut}")
            pdb_wt = res_id_map[resnum_int].get_resname()
            pdb_wt_one = three_to_one[pdb_wt]
            if pdb_wt_one != wt:
                print(f"⚠️ Warning: PDB WT {pdb_wt_one} != mutation WT {wt} at {resnum}")
            f.write("1\n")
            f.write(f"{pdb_wt_one} {resnum} {mutaa}\n")


def run_rosetta_ddg(pdb_file):
    os.makedirs("_rosetta", exist_ok=True)
    cmd = [
        os.path.expanduser(ROSETTA_EXE),
        "-in:file:s", pdb_file,
        "-ddg:mut_file", "mut_file.txt",
        "-ddg:iterations", "1",
        "-score:weights", "ref2015_cart",
        "-fa_max_dis", "9.0",
        "-ex1", "-ex2",
        "-use_input_sc",
        "-flip_HNQ",
        "-no_optH", "false",
        "-ignore_unrecognized_res", "true",
        "-ddg:dump_pdbs", "false",
        "-ddg:mean", "true",
        "-ddg:output_silent", "true",
        "-out:file:scorefile", "_rosetta/ddg_results.sc",
        "-out:suffix", "_ddg",
        "-renumber_pdb", "false"
    ]
    with open("_rosetta/rosetta_log.txt", "w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=log)
    if result.returncode != 0:
        raise RuntimeError("❌ Rosetta ddg run failed. Check _rosetta/rosetta_log.txt.")


def parse_rosetta_ddg(ddg_file, pdb_file):
    if not os.path.exists(ddg_file):
        raise FileNotFoundError("Rosetta ddg file not found.")
    wt_energies = []
    mut_energies = {}
    with open(ddg_file, "r") as f:
        for line in f:
            if line.startswith("COMPLEX:"):
                parts = line.split()
                tag = parts[2].strip(":")
                score = float(parts[3])
                if tag == "WT":
                    wt_energies.append(score)
                elif tag.startswith("MUT_"):
                    mut = tag.replace("MUT_", "")
                    mut_energies.setdefault(mut, []).append(score)
    if not wt_energies:
        raise RuntimeError("No WT energies found in Rosetta ddg file.")
    avg_wt = sum(wt_energies) / len(wt_energies)
    negative_ddg = []
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    model = next(structure.get_models())
    chain = model.child_list[0]
    res_id_map = {res.id[1]: res for res in chain if res.id[0] == ' '}
    for mut, vals in mut_energies.items():
        avg_mut = sum(vals) / len(vals)
        ddg_reu = avg_mut - avg_wt
        ddg_kcal = ddg_reu / 2.94  
        if ddg_kcal < 0:
            resnum = int(re.findall(r'\d+', mut)[0])
            mut3 = mut.replace(str(resnum), "").upper()
            mutaa = three_to_one.get(mut3[:3], "?")
            wt_one = three_to_one[res_id_map[resnum].get_resname()]
            negative_ddg.append((f"{wt_one}{resnum}{mutaa}", ddg_kcal))
    return negative_ddg


# ----------------------------------------------------- #

def create_run_directory(pdb_file):
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    dir_name = base_name
    counter = 1
    while os.path.exists(dir_name):
        dir_name = f"{base_name}_{counter}"
        counter += 1
    os.makedirs(dir_name, exist_ok=True)
    return dir_name


def main():
    parser = argparse.ArgumentParser(description="Mutational scanning pipeline using FreeSASA, FoldX, and optionally Rosetta.")
    parser.add_argument('pdb_file', help='Path to the PDB file')
    parser.add_argument('--ifHydrophobic', action='store_true')
    parser.add_argument('--sasaThreshold', type=float, default=1.0)
    parser.add_argument('--functionalRes', nargs='+')
    parser.add_argument('--contactShell', type=int, choices=[1, 2], default=1)
    parser.add_argument('--contactCutoff', type=float, default=4.5)
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--scorer', choices=['FoldX', 'AddRosetta'], default='FoldX',
                        help="Choose scoring method: FoldX (default) or AddRosetta (FoldX + Rosetta)")
    args = parser.parse_args()

    pdb_file = args.pdb_file
    run_dir = create_run_directory(pdb_file)
    original_dir = os.getcwd()
    os.chdir(run_dir)
    shutil.copy(os.path.join(original_dir, pdb_file), '.')
    cleaned_pdb_file = clean_pdb(pdb_file)
    sasa_output = run_freesasa(cleaned_pdb_file)
    core_residues = parse_freesasa_output(sasa_output, threshold=args.sasaThreshold)
    excluded_residues = set()
    if args.functionalRes:
        excluded_residues = get_excluded_residues(
            cleaned_pdb_file,
            args.functionalRes,
            distance_threshold=args.contactCutoff,
            shell_level=args.contactShell
        )
    mutations, mutation_dict = generate_mutation_list(core_residues, if_hydrophobic=args.ifHydrophobic,
                                                     excluded_residues=excluded_residues)

    # Run FoldX
    run_foldx(cleaned_pdb_file)
    foldx_negative = parse_foldx_output(cleaned_pdb_file, mutations)
    print("\n📊 Mutations with negative ddG values (FoldX):")
    for mutation, ddg in sorted(foldx_negative, key=lambda x: x[1]):
        mut = mutation.strip(";")
        wt = mut[0]
        resnum = mut[2:-1]
        mutaa = mut[-1]
        print(f"{wt}{resnum}{mutaa}: {ddg:.3f} kcal/mol")

    # Run Rosetta if requested
    if args.scorer == "AddRosetta":
        generate_rosetta_mutfile(mutations, cleaned_pdb_file)
        run_rosetta_ddg(cleaned_pdb_file)
        rosetta_negative = parse_rosetta_ddg("mut_file.ddg", cleaned_pdb_file)
        print("\n📊 Mutations with negative ddG values (Rosetta):")
        for mutation, ddg in sorted(rosetta_negative, key=lambda x: x[1]):
            print(f"{mutation}: {ddg:.3f} kcal/mol")

    os.chdir(original_dir)


if __name__ == '__main__':
    main()

