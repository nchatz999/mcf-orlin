#!/usr/bin/env python3
"""
Δημιουργία δικτύων δοκιμών για αλγορίθμους minimum cost flow.

Χρήση:
    python generator.py <nodes> <sources> <sinks> <supply> <output> [--seed N]

Παράδειγμα:
    python generator.py 10 2 2 100 network.txt
"""

import networkx as nx
import random
import argparse
import sys


def make_balanced_flow_network(n_nodes, n_sources, n_sinks, total_supply):
    nodes = list(range(n_nodes))
    all_possible = [(u, v) for u in nodes for v in nodes if u != v]

    while True:
        G = nx.DiGraph()
        G.add_nodes_from(nodes)
        edges = random.sample(all_possible, len(all_possible) // 2)
        G.add_edges_from(edges)
        if nx.is_strongly_connected(G):
            break

    sources = nodes[:n_sources]
    sinks = nodes[-n_sinks:]

    remaining = total_supply
    for i, s in enumerate(sources[:-1]):
        amt = random.randint(1, remaining - (len(sources) - i - 1))
        G.nodes[s]['demand'] = -amt
        remaining -= amt
    G.nodes[sources[-1]]['demand'] = -remaining

    remaining = total_supply
    for i, t in enumerate(sinks[:-1]):
        amt = random.randint(1, remaining - (len(sinks) - i - 1))
        G.nodes[t]['demand'] = amt
        remaining -= amt
    G.nodes[sinks[-1]]['demand'] = remaining

    for node in nodes[n_sources:-n_sinks]:
        G.nodes[node]['demand'] = 0

    for u, v in G.edges():
        G[u][v]['weight'] = random.randint(1, 100)

    return G, sources, sinks


def export_to_file(G, filename):
    n = G.number_of_nodes()
    m = G.number_of_edges()

    with open(filename, 'w') as f:
        f.write(f"{n} {m}\n")
        supplies = [str(-G.nodes[node].get('demand', 0)) for node in range(n)]
        f.write(" ".join(supplies) + "\n")
        for u, v in G.edges():
            f.write(f"{u} {v} {G[u][v]['weight']}\n")


def main():
    parser = argparse.ArgumentParser(description='Δημιουργία ισορροπημένου δικτύου ροής')
    parser.add_argument('n_nodes', type=int, help='Αριθμός κόμβων')
    parser.add_argument('n_sources', type=int, help='Αριθμός κόμβων πηγής')
    parser.add_argument('n_sinks', type=int, help='Αριθμός κόμβων κατανάλωσης')
    parser.add_argument('total_supply', type=int, help='Συνολικό supply (= συνολικό demand)')
    parser.add_argument('output_file', type=str, help='Αρχείο εξόδου')
    parser.add_argument('--seed', type=int, default=None, help='Seed τυχαιότητας')
    args = parser.parse_args()

    if args.n_sources + args.n_sinks > args.n_nodes:
        print("Σφάλμα: n_sources + n_sinks πρέπει να είναι <= n_nodes", file=sys.stderr)
        sys.exit(1)

    if args.seed is not None:
        random.seed(args.seed)

    G, sources, sinks = make_balanced_flow_network(
        args.n_nodes,
        args.n_sources,
        args.n_sinks,
        args.total_supply
    )

    optimal_cost = nx.min_cost_flow_cost(G)
    export_to_file(G, args.output_file)

    print(f"Δημιουργήθηκε: {args.output_file}")
    print(f"  Κόμβοι: {args.n_nodes}, Ακμές: {G.number_of_edges()}")
    print(f"  Πηγές: {sources}, Καταναλωτές: {sinks}")
    print(f"  Συνολικό supply: {args.total_supply}")
    print(f"  Βέλτιστο κόστος: {optimal_cost}")


if __name__ == "__main__":
    main()
