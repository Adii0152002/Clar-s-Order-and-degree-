import networkx as nx
from collections import defaultdict

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
            connectivity.setdefault(atom, set()).update(neighbors)
            for n in neighbors:
                connectivity.setdefault(n, set()).add(atom)
    return max_atom_index, {k: sorted(v) for k, v in connectivity.items()}

# #  Matching-Kekule 
def has_perfect_matching(G: nx.Graph) -> bool:
    if len(G) % 2 != 0:
        return False
    M = nx.algorithms.matching.max_weight_matching(G, maxcardinality=True)
    return len(M) * 2 == len(G)

#  Hexagon (face) detection 
def chordless_hexagons(G: nx.Graph):
    hexes = set()
    adj = {u: set(G.neighbors(u)) for u in G.nodes()}

    def is_chordless(cycle):
        
        n = 6
        for i in range(n):
            for j in range(i+1, n):
                if (j == (i+1) % n) or (i == (j+1) % n):
                    continue
                if cycle[j] in adj[cycle[i]]:
                    return False
        return True

    for start in G.nodes():
        stack = [(start, [start])]
        while stack:
            node, path = stack.pop()
            if len(path) < 6:
                for nbr in G[node]:
                    if nbr == start:
                        continue
                    if nbr in path:
                        continue
                    stack.append((nbr, path + [nbr]))
            elif len(path) == 6:
                if start in G[node]:
                    cycle = path[:]
                    if len(set(cycle)) != 6:
                        continue
                    if not is_chordless(cycle):
                        continue
                    def rotations(t):
                        return [tuple(t[i:] + t[:i]) for i in range(6)]
                    rots = rotations(cycle)
                    rrev = rotations(list(reversed(cycle)))
                    canon = min(rots + rrev)
                    hexes.add(canon)
    return sorted(hexes)

#  Sextet Polynomial 
class SextetPolynomial:
    def __init__(self, graph: nx.Graph):
        self.graph = graph.copy()
        self.hexagons = chordless_hexagons(self.graph)
        self.result = defaultdict(int)
        self.clar_order = None
        self.clar_degree = None

    def remove_sextet(self, G: nx.Graph, cycle):
        H = G.copy()
        H.remove_nodes_from(cycle)
        return H

    def traverse(self, G: nx.Graph, idx_start=0, used_vertices=None, chosen=None):
        if used_vertices is None:
            used_vertices = set()
        if chosen is None:
            chosen = []

        # IMPORTANT: count the current chosen set if the leftover graph is Kekulé.
        # This counts all Clar covers (including the empty set).
        if has_perfect_matching(G):
            self.result[len(chosen)] += 1

        # Try to extend by adding further vertex-disjoint hexagons
        for i in range(idx_start, len(self.hexagons)):
            h = self.hexagons[i]
            Hset = set(h)
            if used_vertices & Hset:
                continue
            G2 = self.remove_sextet(G, h)
            used2 = used_vertices | Hset
            chosen2 = chosen + [tuple(sorted(h))]
            self.traverse(G2, i + 1, used2, chosen2)

    def compute(self):
        self.traverse(self.graph, 0, set(), [])
        if self.result:
            self.clar_order = max(self.result.keys())
            self.clar_degree = self.result[self.clar_order]
        else:
            self.clar_order = 0
            self.clar_degree = 0
        return dict(self.result)

#  Main 
if __name__ == "__main__":
    filename = "coronene.txt"
    _, connectivity = read_connectivity(filename)

    G = nx.Graph()
    for atom, nbrs in connectivity.items():
        for n in nbrs:
            if atom < n:
                G.add_edge(atom, n)

    sp = SextetPolynomial(G)
    poly = sp.compute()

    print("Detected (chordless) hexagons:", sp.hexagons)
    print("Sextet polynomial coefficients:", poly)
    print("Clar order (max sextets):", sp.clar_order)
    print("Clar degree (number of Clar structures):", sp.clar_degree)

    if poly:
        terms = []
        for k in sorted(poly.keys()):
            v = poly[k]
            if k == 0:
                terms.append(str(v))
            elif k == 1:
                terms.append("x" if v == 1 else f"{v}*x")
            else:
                terms.append(f"{v}*x^{k}")
        print("S(H; x) =", " + ".join(terms))
    else:
        print("S(H; x) = 0")
        
        



