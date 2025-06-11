# CoreStab

Mutation pipeline to stabilize buried core.

Requires FreeSASA and Foldx installed and added to bashrc

usage: CoreStab.py [-h] [--ifHydrophobic] [--sasaThreshold SASATHRESHOLD] [--functionalRes FUNCTIONALRES [FUNCTIONALRES ...]]
                   [--contactShell {1,2}] [--contactCutoff CONTACTCUTOFF] [--verbose]
                   pdb_file

Mutational scanning pipeline using FreeSASA and FoldX.

positional arguments:
  pdb_file              Path to the PDB file

optional arguments:
  -h, --help            show this help message and exit
  --ifHydrophobic       Restrict mutations to more hydrophobic residues only
  --sasaThreshold SASATHRESHOLD
                        SASA threshold to define core residues (default=1.0)
  --functionalRes FUNCTIONALRES [FUNCTIONALRES ...]
                        Functionally important residue numbers (e.g., 45 128 193)
  --contactShell {1,2}  Contact shell to exclude: 1 (default) or 2
  --contactCutoff CONTACTCUTOFF
                        Distance cutoff in Å for defining contacts (default=4.5)
  --verbose             Print verbose logs
