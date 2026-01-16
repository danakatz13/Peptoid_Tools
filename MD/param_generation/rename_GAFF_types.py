import networkx as nx
from networkx.algorithms import isomorphism
import numpy as np
import sys
import os
import glob

OUTPUT_SUFFIX = "_renamed.pdb"  # Result will be inputname_renamed.pdb

def get_element_from_type(atom_type):
    """Infers element from GAFF atom types."""
    atom_type = atom_type.lower()
    if atom_type.startswith('cl'): return 'Cl'
    if atom_type.startswith('br'): return 'Br'
    if atom_type.startswith('c'): return 'C'
    if atom_type.startswith('n'): return 'N'
    if atom_type.startswith('o'): return 'O'
    if atom_type.startswith('h'): return 'H'
    if atom_type.startswith('s'): return 'S'
    if atom_type.startswith('p'): return 'P'
    if atom_type.startswith('f'): return 'F'
    if atom_type.startswith('i'): return 'I'
    return 'X'

def load_mol2_graph(mol2_file):
    """Parses a .mol2 file into a NetworkX graph."""
    G = nx.Graph()
    with open(mol2_file, 'r') as f:
        lines = f.readlines()

    section = None
    for line in lines:
        line = line.strip()
        if line.startswith("@<TRIPOS>ATOM"):
            section = "ATOM"
            continue
        elif line.startswith("@<TRIPOS>BOND"):
            section = "BOND"
            continue
        elif line.startswith("@<TRIPOS>"):
            section = None
            continue

        if section == "ATOM" and line:
            parts = line.split()
            atom_id = int(parts[0])
            name = parts[1]
            atom_type = parts[5]
            element = get_element_from_type(atom_type)
            G.add_node(atom_id, element=element, name=name)
            
        elif section == "BOND" and line:
            parts = line.split()
            a1 = int(parts[1])
            a2 = int(parts[2])
            G.add_edge(a1, a2)
    return G

def parse_pdb_structure(pdb_file):
    """Parses PDB into residues and connectivity."""
    residues = {} 
    connectivity = []
    
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            serial = int(line[6:11])
            name = line[12:16].strip()
            resName = line[17:20].strip()
            chain = line[21:22].strip()
            resSeq = int(line[22:26])
            
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            
            element = line[76:78].strip()
            if not element:
                element = ''.join([c for c in name if c.isalpha()])[:1]

            key = (resName, chain, resSeq)
            if key not in residues:
                residues[key] = {'atoms': [], 'coords': {}, 'elements': {}, 'serials': []}
            
            residues[key]['atoms'].append(line)
            residues[key]['coords'][serial] = np.array([x, y, z])
            residues[key]['elements'][serial] = element
            residues[key]['serials'].append(serial)

        elif line.startswith("CONECT"):
            parts = line.split()
            source = int(parts[1])
            for target_str in parts[2:]:
                target = int(target_str)
                connectivity.append((source, target))

    return residues, connectivity, lines

def build_pdb_residue_graph(residue_data, full_connectivity):
    """Builds a subgraph for a single PDB residue."""
    G = nx.Graph()
    serials = residue_data['serials']
    coords = residue_data['coords']
    elements = residue_data['elements']
    
    for s in serials:
        G.add_node(s, element=elements[s])

    # Filter connectivity to only include bonds strictly within this residue
    relevant_bonds = 0
    serial_set = set(serials)
    for (s1, s2) in full_connectivity:
        if s1 in serial_set and s2 in serial_set:
            G.add_edge(s1, s2)
            relevant_bonds += 1


    if relevant_bonds == 0 and len(serials) > 1:
        s_list = list(serials)
        for i in range(len(s_list)):
            for j in range(i + 1, len(s_list)):
                s1, s2 = s_list[i], s_list[j]
                dist = np.linalg.norm(coords[s1] - coords[s2])
                is_h = 'H' in [elements[s1], elements[s2]]
                threshold = 1.6 if is_h else 1.95
                if dist < threshold:
                    G.add_edge(s1, s2)
    return G

def find_matching_mol2(resName, available_mol2s):
    # 1. Exact Name match (e.g. "006.mol2")
    for m in available_mol2s:
        if os.path.basename(m).split('.')[0] == resName:
            return m
            
    # 2. Contains Name (e.g. "006.gaff2.mol2")
    for m in available_mol2s:
        if resName in os.path.basename(m):
            return m
            
    return None

def process_single_pdb(pdb_file, available_mol2s):
    print(f"\nProcessing {pdb_file}...")
    residues, connectivity, all_lines = parse_pdb_structure(pdb_file)
    

    serial_to_line_idx = {}
    for idx, line in enumerate(all_lines):
        if line.startswith(("ATOM", "HETATM")):
            serial = int(line[6:11])
            serial_to_line_idx[serial] = idx

    print(f"  Found {len(residues)} residues.")
    

    loaded_mol2_graphs = {}

    for res_key, res_data in residues.items():
        resName, chain, resSeq = res_key
        
        mol2_path = find_matching_mol2(resName, available_mol2s)
        
        if not mol2_path:
            print(f"  ⚠ Skipping {resName}-{resSeq}: No matching .mol2 file found in directory.")
            continue

        if mol2_path not in loaded_mol2_graphs:
            loaded_mol2_graphs[mol2_path] = load_mol2_graph(mol2_path)
        
        G_mol2 = loaded_mol2_graphs[mol2_path]
        G_pdb = build_pdb_residue_graph(res_data, connectivity)

        nm = isomorphism.categorical_node_match("element", "X")
        GM = isomorphism.GraphMatcher(G_mol2, G_pdb, node_match=nm)

        if GM.subgraph_is_isomorphic():
            mapping = GM.mapping # {Mol2_ID: PDB_Serial}
            
            for mol2_id, pdb_serial in mapping.items():
                new_name = G_mol2.nodes[mol2_id]['name']
                
                line_idx = serial_to_line_idx[pdb_serial]
                old_line = all_lines[line_idx]
                
                # Format Name (Center aligned 4 chars)
                if len(new_name) == 4: formatted = new_name
                else: formatted = f" {new_name:<3}"[:4]

                new_line = old_line[:12] + f"{formatted:^4}" + old_line[16:]
                all_lines[line_idx] = new_line
            
            # print(f"    ✔ Renamed {resName}-{resSeq} using {os.path.basename(mol2_path)}")
        else:
            print(f"    ✘ FAILED {resName}-{resSeq}. Structure does not match {os.path.basename(mol2_path)}")

    base_name = os.path.splitext(pdb_file)[0]
    output_name = f"{base_name}{OUTPUT_SUFFIX}"
    
    with open(output_name, 'w') as f:
        f.writelines(all_lines)
    print(f"  Done. Saved to {output_name}")

def main():
    # 1. Find all PDB files (excluding already renamed ones)
    all_pdbs = glob.glob("*.pdb")
    pdbs_to_process = [p for p in all_pdbs if OUTPUT_SUFFIX not in p]
    
    if not pdbs_to_process:
        print("No .pdb files found in this directory!")
        return

    # 2. Find all Mol2 files (library of fragments)
    available_mol2s = glob.glob("*.mol2")
    if not available_mol2s:
        print("No .mol2 files found! Cannot rename atoms without templates.")
        return
        
    print(f"Found {len(available_mol2s)} template (.mol2) files.")

    # 3. Process each PDB
    for pdb_file in pdbs_to_process:
        process_single_pdb(pdb_file, available_mol2s)

if __name__ == "__main__":
    main()
