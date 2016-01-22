__author__ = 'etseng@pacb.com'
import pdb
import splice_graph as sp

def detect_and_replace_cycle(G, path_d, weight_d, mermap, max_node_index, kmer_size, debug=False):
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
                    # add in sanity check
                    oldseq = sp.stitch_string_from_path(path, mermap)
                    has_cycle = True
                    if debug:
                        pdb.set_trace()
                    # find first and last occurrence of n
                    # NOTE: need to handle special case when i=0 and j=last
                    i = path.index(n)
                    j = len(path) - path[::-1].index(n)  # is actually one beyond the last occurrence
                    newnode = max_node_index + 1
                    max_node_index += 1
                    newmer = mermap[path[i]]
                    for p in path[i+1:j]: newmer += mermap[p][kmer_size-1:]
                    mermap[newnode] = newmer
                    # update G before update path
                    if i > 0: # if i==0, then nothing to connect with prev
                        G.add_edge(path[i-1], newnode, weight=weight_d[seqid])
                    if j < len(path): # if j=last, then nothing to connect with after
                        G.add_edge(newnode, path[j], weight=weight_d[seqid])
                    for k in xrange(max(0,i-1), min(j, len(path)-1)):
                        s, t = path[k], path[k+1]
                        G[s][t]['weight'] -= weight_d[seqid]
                    # now we can update the path
                    path_d[seqid] = path[:i] + [newnode] + path[j:]
                    path = path_d[seqid]


                    assert sp.stitch_string_from_path(path, mermap) == oldseq
                    break

    # clean up zero-weight edges and zero degree nodes
    # TODO here
    for s,t,d in G.edges(data=True):
        if d['weight']==0:
            G.remove_edge(s, t)
    bad = []
    for n in G.nodes_iter():
        if G.degree(n) == 0: bad.append(n)
    G.remove_nodes_from(bad)

