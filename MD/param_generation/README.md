# How to parametize non-standard amino acids in AMBER

## Generating non-standard amino acid residue
Need to setup a custom enviorment (**add requirments.txt**)
```
module purge
module load anaconda3/2020.07
module load amber/openmpi/intel/22.03
conda activate /scratch/dk4128/penv
```
### Use SMILES and generate_amber_mol2.py
To generate these consistent structure files, we use their SMILES, written as either a primary amine (peptoid) or a peptide residue (NC(R)C(=O)) for a non-canonical amino acid.

generate_amber_mol2 will generate a .mol of the capped residue. 
```
python generate_amber_mol2.py smiles.csv -o output_dir
```
smiles.csv has 2 columns : code + SMILES with comma delimeter. Code should be 3 characters and NEVER the same as a canonical amino acid.
In output_dir you will now have a .mol file for each residue. Originally, it generated the files as PDBs, but it seemed to cause some downstream consequence due to some of the double bonds not being properly labeled (PDBs have no bond order). MOL2 should correctly store bond orders

### Generate GAFF atom types for your capped residue
This step is necessary to view the GAFF atom names that AMBER will give your atoms. Visualizing this in PyMol will be helpful for later steps to match atom names. 

For one single residue at a time : 
```
antechamber -i 003.mol2 -fi mol2 -o 003.gaff2.mol2 -fo mol2 -rn 003 -at gaff2 -an yes -dr no -pf yes -c bcc -nc 0
```
For all residues in a directory : 
```
for f in *.mol2; do      base="${f%.mol2}";     antechamber -i "$f" -fi mol2 -o "${base}.gaff2.mol2" -fo mol2 -rn "$base" -at gaff2 -an yes -dr no -pf yes -c bcc -nc 0; done
```

## Parameterization
These following steps generate partial charges/atom names/bond angles/bond dihedrals. All of these files are necessary for non-canonical amino acids, because we need to add these parameters in to the exsisiting force field so it knows how to handle a non-canonical. 

For one single residue at a time : 
```
antechamber -i 003.mol2 -fi mol2 -bk 003 -fo ac -o 003.ac -c bcc -at amber -nc 0
```
For all residues in a directory : 
```
for f in *.mol2; do base="${f%.mol2}"; antechamber -i "$f" -fi mol2 -bk "$base" -fo ac -o "$base.ac" -c bcc -at amber -nc 0; done
```
After these .ac files are generated, you need to check that the backbone nitrogen was classified as 'N' type. These different atom types will affect the bond lengths and angles. I am not sure why some residues label it as 'NT' or 'N2' and some don't, so I don't know how to automate this step. Just check each .ac file and make sure the first nitrogen listed has an 'N' in the last column. If not change it, and make sure the indentation lines up with the rest of the rows. (007)

Next, we need to generate a mainchain file. This file will tell AMBER, which residues are actually present in the sequence and which are considered 'DUMMY' atoms and will be replaced with the preceeding/proceeding residue. 

This is an example of what atoms are considered what along with the mainchain file **insert peptoid mainchain**

<img width="155" height="212" alt="image" src="https://github.com/user-attachments/assets/c3ec77e1-652e-464d-bae6-68455fac48f4" />
<img width="91" height="215" alt="image" src="https://github.com/user-attachments/assets/d0481570-489a-4c25-a2f3-2371eff3987d" />

**The atom names HAVE to match the atom names in the .gaff2.mol2 structure you generated. View that structure in PyMol to decide which atoms need to be defined in the mainchain file. If you use the decompose2smiles.py script, the names should be consistent and you can use a single mainchain file**

Once you have the mainchain file, run this: 
```
prepgen -i 003.ac -o 003.prepin -m peptoid_mainchain.mc -rn 003
```
```
for f in *.ac; do
    base=$(basename "$f" .ac)
    prepgen -i "$f" -o "${base}.prepin" -m peptoid_mainchain.mc -rn "$base"
done
```
Next to generate .frcmod file:
```
parmchk2 -i 003.prepin -f prepi -o 003.frcmod -a Y
```
```
for f in *.prepin; do
    base=$(basename "$f" .prepin)
    parmchk2 -i "$f" -f prepi -o "${base}.frcmod" -a Y
done
```
These 2 files (.prepin and .frcmod) describe all the parameters needed for AMBER to build a non-standard residue. Most errors will involve missing parameters for these 2 files.

## Preparing Structure

To clean up the .pdb:
```
pdb4amber -i {name}_truncated.pdb -o {name}_clean.pdb --dry --reduce
```
To change the PDB atom types in the clean structure to match the GAFF atom types AMBER uses for parameterization, the python script rename_GAFF_types.py uses the atom names in the .gaff2.mol2 file generated and matches them to atoms in PDB structure based on their residue code and creates an atom map. This is a quick (mostly) automated way to rename the atoms to match AMBER GAFF types.
In a directory with your prepared structure, add the .gaff2.mol2 files and run this script:
```
rename_GAFF_types.py
```
It matches the PDB file in the directory and the .gaff2.mol2 files and prepares a new structure named file_renamed.pdb. This is now the structure you should use in tleap.

**Note: For some reason, the GAFF atom type for chlorine atoms is 'CL' but in the prepin file, its generated as 'Cl'. In the PDB structure, just change the 'CL' in the atom name to 'Cl'. The atom number is unchanged (ex. CL16 --> Cl16)
**Note: For some reason, for residues with an aromatic ring directly attached to Nitrogen, it incorrectly reads the 'CA' and has some parameters missing. Add this to .frcmod for that residue, and that should eliminate any problems.

Under torsions:
```
CA-N -C     70.0   120.0
C -N -CA    70.0   120.0
```

Under dihedrals: 
```
CA-CA-N -C    4   1.800   180.000   2.000
C -N -CA-CA   4   1.800   180.000   2.000
```
### A third-level heading

Style	Syntax	Keyboard shortcut	Example	Output
Bold	** ** or __ __	Command+B (Mac) or Ctrl+B (Windows/Linux)	**This is bold text**	This is bold text
Italic	* * or _ _     	Command+I (Mac) or Ctrl+I (Windows/Linux)	_This text is italicized_	This text is italicized
Strikethrough	~~ ~~ or ~ ~	None	~~This was mistaken text~~	This was mistaken text
Bold and nested italic	** ** and _ _	None	**This text is _extremely_ important**	This text is extremely important
All bold and italic	*** ***	None	***All this text is important***	All this text is important
Subscript	<sub> </sub>	None	This is a <sub>subscript</sub> text	This is a subscript text
Superscript	<sup> </sup>	None	This is a <sup>superscript</sup> text	This is a superscript text
Underline	<ins> </ins>	None	This is an <ins>underlined</ins> text	This is an underlined text
Some basic Git commands are:
```
git status
git add
git commit
```


