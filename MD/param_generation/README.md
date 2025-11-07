# How to parametize non-standard amino acids in AMBER

## Generating non-standard amino acid residue
To start, we need to generate a .mol2 file of the capped residue by itself. If you have this available, skip this step.

There are 2 ways you can generate these pdb/.mol2 capped residue files. If you use #2, they should have a consistent atom naming scheme and generating a mainchain file in a later step should be the same for all residues. If you have many different residues, I suggest you use this approach. 

### Option #1 : Use ChemDraw
### Option #2 : Use SMILES and decompose2smiles.py
To generate these consistent structure files, we use their SMILES, written as either a primary amine (peptoid) or a peptide residue (NC(R)C(=O)) for a non-canonical amino acid.

decompse2smiles.py will generate a .pdb of the capped residue. 

### Generate GAFF atom types for your capped residue
This step is necessary to view the GAFF atom names that AMBER will give your atoms. Visualizing this in PyMol will be helpful for later steps to match atom names. 

For one single residue at a time : 
```
antechamber -i ACG.pdb -fi pdb -o ACG.gaff2.mol2 -fo mol2 -rn ACG -at gaff2 -an yes -dr no -pf yes -c bcc -nc 0
```
For all residues in a directory : 
```
for f in *.pdb; do base="${f%.pdb}"; antechamber -i "$f" -fi pdb -o "${base}.gaff2.mol2" -fo mol2 -rn "$base" -at gaff2 -an yes -dr no -pf yes -c bcc -nc 0; done
```

## Parameterization
These following steps generate partial charges/atom names/bond angles/bond dihedrals. All of these files are necessary for non-canonical amino acids, because we need to add these parameters in to the exsisiting force field so it knows how to handle a non-canonical. 

For one single residue at a time : 
```
antechamber -i ACG.pdb -fi pdb -bk ACG -fo ac -o ACG.ac -c bcc -at amber -nc 0
```
For all residues in a directory : 
```
for f in *.pdb; do base="${f%.pdb}"; antechamber -i "$f" -fi pdb -bk "$base" -fo ac -o "$base.ac" -c bcc -at amber -nc 0; done
```
After these .ac files are generated, you need to check that the backbone nitrogen was classified as 'N' type. These different atom types will affect the bond lengths and angles. I am not sure why some residues label it as 'NT' or 'N2' and some don't, so I don't know how to automate this step. Just check each .ac file and make sure the first nitrogen listed has an 'N' in the last column. If not change it, and make sure the indentation lines up with the rest of the rows.

Next, we need to generate a mainchain file. This file will tell AMBER, which residues are actually present in the sequence and which are considered 'DUMMY' atoms and will be replaced with the preceeding/proceeding residue. 

This is an example of what atoms are considered what along with the mainchain file 

<img width="155" height="212" alt="image" src="https://github.com/user-attachments/assets/c3ec77e1-652e-464d-bae6-68455fac48f4" />
<img width="91" height="215" alt="image" src="https://github.com/user-attachments/assets/d0481570-489a-4c25-a2f3-2371eff3987d" />

**The atom names HAVE to match the atom names in the .gaff2.mol2 structure you generated. View that structure in PyMol to decide which atoms need to be defined in the mainchain file. If you use the decompose2smiles.py script, the names should be consistent and you can use a single mainchain file**

Once you have the mainchain file, run this: 
```
prepgen -i ACG.ac -o ACG.prepin -m ACG_mainchain.mc -rn ACG
```
```
for f in *.ac; do
    base=$(basename "$f" .ac)
    prepgen -i "$f" -o "${base}.prepin" -m "${base}_mainchain.mc" -rn "$base"
done
```
Next to generate .frcmod file:
```
parmchk2 -i ACG.prepin -f prepi -o ACG.frcmod -a Y
```
```
for f in *.prepin; do
    base=$(basename "$f" .prepin)
    parmchk2 -i "$f" -f prepi -o "${base}.frcmod" -a Y
done
```
These 2 files (.prepin and .frcmod) describe all the parameters needed for AMBER to build a non-standard residue. Most errors will involve missing parameters for these 2 files.

## Preparing Structure
Okay, this is the most annoying part. We have to manually change the .pdb atom names to match with the GAFF atom names we generated so AMBER can understand.

To clean up the .pdb:
```
pdb4amber -i {name}_truncated.pdb -o {name}_clean.pdb --dry --reduce
```
Now, we need to edit the atom names for any non-canonical residues we built parameters for. I find the easiest way to do this is have the clean structure open in one PyMol window, and the .gaff2.mol2 structures for each residue open in another one. Going through the PDB, you have to change the atom name to match the .gaff2.mol2 structure. 

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


