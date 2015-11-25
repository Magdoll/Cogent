__author__ = 'etseng@pacb.com'
import pdb

def detect_and_replace_cycle(G, path_d, weight_d, mermap, max_node_index, kmer_size):
    """
    For each path in path_d, if the same node occurs multiple times,
    replace the whole occurrence of <n> -> ... -> <n> with a new node

    UGLY but temporary fix :(

    FOr each cycle replacement:
    (a) add prev -> newnode -> succ in G
    (b) reduce one weight in the path and if falls to zero, remove nodes
    (c) update path
    """
    for seqid, path in path_d.iteritems():
        has_cycle = True
        while has_cycle:
            has_cycle = False
            for n in path:
                if path.count(n) > 1:
                    has_cycle = True
                    # find first and last occurrence of n
                    i = path.index(n)
                    j = len(path) - path[::-1].index(n)  # is actually one beyond the last occurrence
                    newnode = max_node_index + 1
                    newmer = mermap[path[i]]
                    for p in path[i+1:j]: newmer += mermap[p][kmer_size-1:]
                    mermap[newnode] = newmer
                    # update G before update path
                    G.add_edge(path[i-1], newnode, weight=weight_d[seqid])
                    G.add_edge(newnode, path[j], weight=weight_d[seqid])
                    for k in xrange(i-1, j):
                        s, t = path[k], path[k+1]
                        G[s][t]['weight'] -= weight_d[seqid]
                        if G[s][t]['weight'] == 0:
                            G.remove_edge(s, t)
                    # now we can update the path
                    path_d[seqid] = path[:i] + [newnode] + path[j:]
                    path = path_d[seqid]
                    break



