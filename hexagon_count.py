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

import networkx as nx


def count_hexagons(connectivity):
    G = nx.Graph()
    for atom, neighbors in connectivity.items():
        for neighbor in neighbors:
            if atom < neighbor:
                G.add_edge(atom, neighbor)

    hexagons = set()

    def dfs(path, start):
        if len(path) == 6:
            if start in G[path[-1]]:
                hexagons.add(tuple(sorted(path)))
            return
        for neighbor in G[path[-1]]:
            if neighbor not in path:
                dfs(path + [neighbor], start)

    for node in G.nodes:
        dfs([node], node)

    return len(hexagons), list(hexagons)


if __name__ == "__main__":
    filename = "coronene.txt" 
    n_atoms, connectivity = read_connectivity(filename)
    print(count_hexagons(connectivity))
    count, hexagons = count_hexagons(connectivity)
    print("Number of hexagons:", count)
