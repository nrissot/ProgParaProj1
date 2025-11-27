from random import randint, shuffle

EXPORT_PATH             : str = "exported_graph.dot"
NUMBER_OF_COMMUNITY     : int = 4
NUMBER_OF_NODES         : int = 100
INTER_COMMUNITY_WEIGHT  : int = 10
INTRA_COMMUNITY_WEIGHT  : tuple[int, int] = (1, 3)

type node  = str
type edge  = tuple[node, node, int]
type graph = tuple[list[node], list[edge]]

def main() :
    # generate a loose graph that will be the basis for the inter-community exchanges
    communities_as_node : list[node] = [ f"C{x}" for x in range(NUMBER_OF_COMMUNITY) ]

    # expand the graph to convert the community nodes into tight knit community of nodes
    final_graph : graph = ([], [])
    ionodes : list[tuple[node,node,node]] = []

    for community in communities_as_node :
        community_nodes : list[node] = []
        for node_number in range(NUMBER_OF_NODES // NUMBER_OF_COMMUNITY) :
            community_nodes.append(f"{community}N{node_number}") # format C2N8
        
        new_edges_wo_weight : list[tuple[node, node]] = []
        for n in community_nodes :
            for i in range(randint((NUMBER_OF_NODES // NUMBER_OF_COMMUNITY)//4, (NUMBER_OF_NODES // NUMBER_OF_COMMUNITY)//2)):
                new_edge_without_weight : tuple[node, node] = (n, community_nodes[randint(0, NUMBER_OF_NODES // NUMBER_OF_COMMUNITY-1)])
                if (new_edge_without_weight[0] != new_edge_without_weight[1] and not new_edge_without_weight in new_edges_wo_weight) :
                    new_edges_wo_weight.append(new_edge_without_weight)
        
        community_edges : list[edge] = []
        for e in new_edges_wo_weight :
            community_edges.append((e[0], e[1], randint(INTRA_COMMUNITY_WEIGHT[0],INTRA_COMMUNITY_WEIGHT[1])))

        final_graph[0].extend(community_nodes)
        final_graph[1].extend(community_edges)
        ionodes.append((community_nodes[0], community_nodes[1]))
    
    for i in range(NUMBER_OF_COMMUNITY) :
        final_graph[1].append((ionodes[i][1], ionodes[(i+1) % NUMBER_OF_COMMUNITY][0], INTER_COMMUNITY_WEIGHT))
    
    shuffle(final_graph[0])
    shuffle(final_graph[1])

    prelude = "graph graphe_genere {\n\tnode [shape=circle, style=filled, color=lightyellow, fontcolor=black];\n\tedge [color=black, fontcolor=blue];\n"
    postlude = "\n}\n"

    # print the graph to a file
    with open(EXPORT_PATH, "w") as out_file :
        print(prelude, file=out_file)
        for n in final_graph[0] :
            print(f"\t{n} [label=\"{n}\"];", file=out_file)
        print("\n", file=out_file)
        for e in final_graph[1] :
            print(f"\t{e[0]} -- {e[1]} [label=\"{e[2]}\", weight={e[2]}];", file=out_file)
        print(postlude, file=out_file)



if __name__ == "__main__" :
    print("Generating a graph")
    main()