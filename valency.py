import sympy as sp

def read_connectivity(filename):
    connectivity = {}
    max_atom_index = 0

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = list(map(int, line.split()))
            atom = parts[0]
            neighbors = parts[1:]

            max_atom_index = max(max_atom_index, atom, *neighbors)


            if atom not in connectivity:
                connectivity[atom] = set()
            connectivity[atom].update(neighbors)

            for neighbor in neighbors:
                if neighbor not in connectivity:
                    connectivity[neighbor] = set()
                connectivity[neighbor].add(atom)

    return max_atom_index, {k: sorted(v) for k, v in connectivity.items()}

def build_symbolic_huckel_matrix(n_atoms, connectivity):
    x = sp.Symbol('x')
    H = sp.zeros(n_atoms)

    for i in range(n_atoms):
        H[i, i] = x  

    for i in range(1, n_atoms+1):
        for j in connectivity.get(i, []):
            H[i-1, j-1] = 1
            H[j-1, i-1] = 1  

    return H

def infer_hybridization(valency):
    if valency == 1:
        return "primary"       
    elif valency == 2:
        return "secondary"      
    elif valency == 3:
        return "tertiary"    
    elif valency == 4:
        return "quarternary"    
    else:
        return "unknown" 
    
def renumber_by_valency(connectivity):
    
    valencies = {atom: len(neighbors) for atom, neighbors in connectivity.items()}
    sorted_atoms = sorted(valencies, key=lambda atom: (-valencies[atom], atom))
    new_numbers = {atom: i + 1 for i, atom in enumerate(sorted_atoms)}
    new_connectivity = {
        new_numbers[atom]: sorted([new_numbers[n] for n in neighbors])
        for atom, neighbors in connectivity.items()
    }
    print(valencies)
    print(sorted_atoms)
    return dict(sorted(new_connectivity.items()))
    
if __name__ == "__main__":
    filename = "coronene.txt" 
    n_atoms, connectivity = read_connectivity(filename)

    print(f"Number of atoms: {n_atoms}")
    print("Connectivity:", connectivity)

    print("Atom Valency and Hybridization:\n")
    for atom in range(1, n_atoms + 1):
        neighbors = connectivity.get(atom, [])
        valency= infer_hybridization(len(neighbors))
        print(f"Atom {atom}: no_of_atom = {len(neighbors)}, valency = {valency}")
     
    new_connectivity = renumber_by_valency(connectivity)

    print("Renumbered Connectivity (sorted by valency):")
    for atom, neighbors in new_connectivity.items():
        print(f"{atom} : {neighbors}")




H = build_symbolic_huckel_matrix(n_atoms, connectivity)

print("\n Huckel matrix:")
sp.pprint(H)

