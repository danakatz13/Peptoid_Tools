# How to clean up / standardize all hits that match our macrocycle definition (macrocycle with peptide like backbone between 7-16 residues)

## Script to find initial hits
macrocycle_finder.py : to be used against CSD entries
macrocycle_finder-noCSD.py : uses SMILES from different .csv files to do same substructure match

## Script to process all hits
peptoid_playbook.py : cleans up and standardizes all smiles, retains one ID column and the source it came from and includes different molecular descriptors
Uses data from raw_hits folder

## Script to analyze dihedral angles of hits from CSD
csd_dihedrals.py : must be run with access to all entries from CSD, because uses ID to retrieve entry and then analyze torsions

## Dihedral angles from CSD hits
csd_backbone_torsions.csv
