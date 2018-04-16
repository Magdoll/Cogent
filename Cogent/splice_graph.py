import pdb
import logging
from Cogent import settings as cc_settings
from Cogent.settings import DEBUG_FLAG
from Cogent import splice_align
from Cogent.process_path import make_in_same_path, stitch_string_from_path, path_finder
from Cogent.sanity_checks import sanity_check_is_chain


log = logging.getLogger('Cogent.splice_graph')

def add_seq_to_graph(G, node_d, path_d, seq, seqid, seq_weight):
    """
    Add a new sequence (seqid + seq) to graph G and update node_d & path_d

    node_d --- dict of sequence --> node index
    path_d --- dict of seqid --> paths through G to reconstruct sequence
    """
    #log.info("Sanity check: in add_seq_to_graph, k-mer size is: {0}".format(cc_settings.KMER_SIZE))
    seq = seq.upper()

    max_node_index = max(node_d.itervalues()) + 1
    
    prev = seq[:cc_settings.KMER_SIZE]
    if prev in node_d:
        s = node_d[prev]
    else:
        s = max_node_index
        max_node_index += 1
        node_d[prev] = s
    path = [s]

    for i in xrange(1, len(seq)-cc_settings.KMER_SIZE+1):
        next = seq[i:i+cc_settings.KMER_SIZE]
        if next in node_d:
            t = node_d[next]
        else:
            t = max_node_index
            max_node_index += 1
            node_d[next] = t
        if G.has_edge(s, t):
            G[s][t]['weight'] += seq_weight
        else:
            G.add_edge(s, t, weight=seq_weight)
        path.append(t)
        s = t
    path_d[seqid] = path


def untangle_homopolymer_helper(G, path_d, mermap, seqweights, node):
    """
    Homopolymers longer than kmer size will be self-edges with no indication how long it is.

    Instead must look at the paths to figure out how long they are.
    In the paths they will be:

    [ ....!node, node, node, node, !node, ...., node, node ... ]

    For each occurrence, replace it with the full homopolymer:

    [ ....!node, new_node1, !node .... new_node2 ... ]

    Homopolymers of the same size can use the same node ID (use mermap to find them).
    G must be updated to reflect the path changes.

    The original node can be deleted if there are no single instances of it in path_d.
    """
    max_node_index = max(G.nodes_iter()) + 1
    homo_rev_mermap = {}
    for k,path in path_d.iteritems():
        weight = seqweights[k]
        if node in path and path.count(node) >= 2:
            # find the longest stretch of repeated node in path
            i = path.index(node)
            while i < len(path) - 1:
                for j in xrange(i+1, len(path)):
                    if path[j]!=node:
                        if j - i == 1: # just a single occurrence, skip over it
                            break

                        # replace w/ path[i-1] -> new_homopolymer_node -> path[j]
                        # also update G by adding path[i-1] -> newnode -> path[j]
                        # we delete node from G & mermap at the very end
                        newmer = mermap[node] + mermap[node][-1]*(j-i-1)
                        if newmer in homo_rev_mermap:
                            newnode = homo_rev_mermap[newmer]
                        else:
                            newnode = max_node_index
                            max_node_index += 1
                            homo_rev_mermap[newmer] = newnode
                            mermap[newnode] = newmer
                        if i > 0:
                            s, t = path[i-1], newnode
                            if G.has_edge(s, t):
                                G[s][t]['weight'] += weight
                            else:
                                G.add_edge(path[i-1], newnode, weight=weight)
                        if j < len(path):
                            s, t = newnode, path[j]
                            if G.has_edge(s, t):
                                G[s][t]['weight'] += weight
                            else:
                                G.add_edge(newnode, path[j], weight=weight)

                        #### DEBUG ####
                        if DEBUG_FLAG:
                            pdb.set_trace()

                        path_d[k] = path_d[k][:i] + [newnode] + path_d[k][j:]
                        path = path_d[k]
                        break
                    elif j == len(path)-1:  # at the end of the path
                        assert j - i > 0
                        # replace w/ path[i-1] -> new_homopolymer_node
                        newmer = mermap[node] + mermap[node][-1]*(j-i)
                        if newmer in homo_rev_mermap:
                            newnode = homo_rev_mermap[newmer]
                        else:
                            newnode = max_node_index + 1
                            max_node_index += 1
                            homo_rev_mermap[newmer] = newnode
                            mermap[newnode] = newmer
                        if i > 0:
                            G.add_edge(path[i-1], newnode, weight=weight)

                        #### DEBUG ####
                        if DEBUG_FLAG:
                            pdb.set_trace()

                        path_d[k] = path_d[k][:i] + [newnode]
                        path = path_d[k]
                        break
                try:
                    i = path.index(node, i+1)
                except ValueError:  # we have searched through the whole path and there are no more of this node
                    break

    # can safely remove the self-edge
    G.remove_edge(node, node)
    can_del_node = True
    # can only remove the node if it is not used anymore in G
    for path in path_d.itervalues():
        if node in path:
            can_del_node = False

    if can_del_node:
        G.remove_node(node)
        del mermap[node]


def contract_sinks(G, path_d, mermap):
    """
    contract all cases where ... -> pred -> sink
    (where n' has only one outgoing to sink, and sink has only one incoming that is n')

    update G to: ... -> new_node where new_node = pred -> sink

    update all paths, which, if contains sink, must be either:

    path = [ ..., pred, sink ] update to [ ..., pred ]
    or
    path = [ sink ] update to [ pred ]
    """
    sinks = filter(lambda n: G.out_degree(n)==0 and G.in_degree(n)==1, G.nodes_iter())
    for sink in sinks:
        pred = G.predecessors(sink)[0]
        # confirm that pred must have no other outgoing edges
        if G.out_degree(pred) > 1: continue # skip it, has other outgoing edges
        # contract it by simply updating pred and removing sink from G
        log.debug("contract sink {0},{1}".format(pred, sink))
        mermap[pred] = mermap[pred] + mermap[sink][(cc_settings.KMER_SIZE-1):]
        # delete sink from G and path_d and mermap
        for k in path_d:
            if sink in path_d[k]:
                if len(path_d[k]) == 1:
                    path_d[k] = [pred]
                else:
                    assert path_d[k][-1] == sink and path_d[k][-2] == pred
                    path_d[k] = path_d[k][:-1]
        # now can safely remove sink from G and mermap
        G.remove_node(sink)
        del mermap[sink]


def find_dangling_sinks(G, path_d, mermap):
    """
    For isoforms w/ 3' alt ends that have a longer last exon, it shows up as a *branch* from the path

    ex:
    ... -> pred -> sink_node (sink_node has only one incoming edge)
    ... -> pred -> n' -> ...other path...
    (pred only has two outgoing edges)

    we can "tuck" sink into n' if and only if sink is a substring of n'

    pred = [prefix] + [suffix]
    sink_node = [suffix] + [extra]
    n' = [suffix] + [...]

    update by deleting sink_node, and updating:
    pred = [prefix] + [suffix] + [extra]
    n' = [just use last k-mer of extra] + [...]
    """
    cand_sinks = filter(lambda n: G.out_degree(n)==0 and G.in_degree(n)==1, G.nodes_iter())
    for sink in cand_sinks:
        pred = G.predecessors(sink)[0]
        for n in G.successors(pred):
            if n == sink or n not in G: continue
            if splice_align.node_is_similar(mermap[sink], mermap[n][:len(mermap[sink])]):
                log.debug("tugging dangling sink: {0}->{1}(sink), {0}->{2}".format(pred, sink, n))
                # sink is just a shortened version of <n>
                # just update all paths with presence of <sink> to <n>
                # and safely remove <sink> from G
                for k in path_d:
                    if sink in path_d[k]:
                        assert path_d[k][-1] == sink
                        path_d[k] = path_d[k][:-1] + [n]
                del mermap[sink]
                G.remove_node(sink)
                break


def collapse_chain(G, chain, mermap, path_d):
    """
    chain -- list of nodes that consist of the unipath

    contract u -> n1 -> n2 ... -> nk -> v down to  u -> v
    n1, n2, .. nk must have exactly one incoming and one outgoing edge

    (1) create a new mer that reprents n1 -> n2 -> ... -> nk
    (2) update G: replace all instances of u -> .. -> v with u -> v and remove the collapsed nodes
                 (weight is updated as the max weight in the path)
    (3) update path_d: replace all path instances down to u -> v
    """
    if len(chain) <= 2: return   # nothing to do

    sanity_check_is_chain(G, chain)

    weights_in_path = [G.get_edge_data(chain[0], chain[1])['weight']]
    n0 = chain[0]
    newmer = mermap[chain[0]]
    for i in xrange(1, len(chain)-1):
        n = chain[i]
        newmer += mermap[n][(cc_settings.KMER_SIZE-1):]
        weights_in_path.append(G.get_edge_data(chain[i], chain[i+1])['weight'])
        del mermap[n]
        G.remove_node(n)
    mermap[n0] = newmer
    nlast = chain[-1]
    avg_weights_in_path = sum(weights_in_path)/len(weights_in_path)
    G.add_edge(n0, nlast, weight=avg_weights_in_path)
    log.debug("shortening {0} -> {1} weight avg: {2}".format(n0, nlast, avg_weights_in_path))

    # updating the paths
    # note that the paths may include only part of the collapsed path
    # ex:  .... -> some_node -> n2 -> n3 -> some_node -> ...
    # which should be updated to ... -> some_node -> u -> v -> some_node -> ...
    for k in path_d:
        oldp = path_d[k]
        if len(oldp) == 1 and oldp[0] in chain[1:-1]: # special case
            path_d[k] = [n0, nlast]
            continue

        newp = []   # use this to hold the first part of the path
        newp2 = []  # use this to hold the second part of the path
        for i,x in enumerate(oldp):
            if x not in chain: newp.append(x)
            else: # first encounter of something in chain
                break
        for j in xrange(len(oldp)-1, i-1, -1):
            x = oldp[j]
            if x not in chain: newp2.append(x)
            else: 
                break
        newp2 = newp2[::-1]  # make it back into the right order
        if i > j or i == len(oldp)-1:
            continue #nothing to do
        elif i == j:  # ex: ... -> some_node -> {must be u or v} -> some_node -> ...
            assert oldp[j] in (n0, nlast)
            newp = oldp  # no changes required
        else: # i < j
            newp += [n0, nlast] + newp2
        path_d[k] = newp



def reachability(G, mermap, visited,  path_d):
    """
    Find all unipaths and contract them!
    G must NOT have cycles so we can successfully traverse through all of them!
    """
    sources = filter(lambda n: G.in_degree(n)==0, G.nodes_iter())
    for source in sources:
        reachability_helper(G, source, [], visited, mermap, path_d)


def reachability_helper(G, cur, chain, visited, mermap, path_d):
    """
    a chain started must have outdeg = 1  (indeg does not matter)
    a chain ender must have indeg = 1  (outdeg does not matter)
    a chain member (neither start/end) must have indeg = 1 and outdeg = 1
    """
    if cur in visited:
        if len(chain) >= 2:
            log.debug("chain found! {0}".format(chain + [cur]))
            collapse_chain(G, chain+[cur], mermap, path_d)
        return
    visited[cur] = 1
    indeg = G.in_degree(cur)
    outdeg = G.out_degree(cur)
    
    if indeg == 1: 
        if outdeg == 0 or outdeg > 1:
            # chain ender
            if len(chain) >= 2:
                log.debug("chain found! {0}".format(chain + [cur]))
                collapse_chain(G, chain+[cur], mermap, path_d)
            # start another possible chain, does not include itself because outdeg is not 1
            for n in G.successors_iter(cur):
                reachability_helper(G, n, [], visited, mermap, path_d)
        else: # outdeg == 1, continue the chain
            reachability_helper(G, G.successors_iter(cur).next(), chain + [cur], visited, mermap, path_d)
    elif outdeg == 1: # indeg is 0 or > 1, outdeg == 1
        if len(chain) >= 2:
            log.debug("chain found! {0}".format(chain + [cur]))
            collapse_chain(G, chain+[cur], mermap, path_d)
        # possible chain starter, includes itself becuz outdeg is 1
        reachability_helper(G, G.successors_iter(cur).next(), [cur], visited, mermap, path_d)
    else: # outdeg!=1, indeg!=1
        for n in G.successors_iter(cur):
            reachability_helper(G, n, [], visited, mermap, path_d)


def find_source_bubbles(G, path_d, mermap):
    """
    Find all cases where
       src1 --> n3 
       ...> path1 --> n3
    and that
    <i> src1 and path1 each has only one outgoing edge to n3
    <ii> src1 and path1 are similar

    path1: can also be a source
    """
    def traceback(cur):
        """
        Retrace path of n1 -> n2 ... -> cur
        where n1, n2....all have exactly one outgoing edge
        """
        acc = []
        while True:
            acc.append(cur)
            preds = G.predecessors(cur)
            if len(preds) == 0 or len(preds) > 1 or G.out_degree(preds[0]) > 1:
                break
            cur = preds[0]
        return acc[::-1]

    def replace_node(n_to_del, path_to_replace_with):
        for k in path_d:
            if n_to_del in path_d[k]:
                i = path_d[k].index(n_to_del)
                path_d[k] = path_d[k][:i] + path_to_replace_with + path_d[k][i+1:]
        G.remove_node(n_to_del)
        del mermap[n_to_del]

    in_same_path = make_in_same_path(path_d)
    sources = filter(lambda n: G.in_degree(n) == 0, G.nodes_iter())
    for src1 in sources:
        if src1 not in G: continue # deleted in the loop below
        succ = G.successors(src1)
        if len(succ) == 1:
            n3 = succ[0]
            cands = G.predecessors(n3)
            for n in cands:
                if src1 not in G: break  # deleted, jump out of this
                if n not in G: continue # deleted in the loop below
                if n!=src1 and n not in in_same_path[src1] and (n in sources or G.out_degree(n)==1):
                    t = traceback(n)
                    seq1 = mermap[src1]
                    seq2 = stitch_string_from_path(t, mermap)
                    minlen = min(len(seq1), len(seq2))
                    # we know that seq1 and seq2 both have the same successor so they must share the same last (KMER_SIZE-1) suffix
                    if DEBUG_FLAG:
                        pdb.set_trace()
                    if splice_align.node_is_similar(seq1[::-1][:minlen], seq2[::-1][:minlen]):
                        # should collapse src1 into n
                        # to do so: 
                        # (1) if both are sources, replace the shorter src with the longer src
                        # (2) if one is src other is path, replace the src with path
                        # -- delete the replaced node from G
                        # -- for all path in path_d that uses the deleted node, update with replacement
                        log.debug("should collapse {0},{1}".format(src1, t))
                        if len(t) == 1 and t[0] in sources: # both are sources
                            if len(mermap[t[0]]) > len(mermap[src1]):
                                replace_node(n_to_del=src1, path_to_replace_with=t)
                            else: # src1 is longer, use src1
                                replace_node(t[0], [src1])
                        else: # src1 is a source, <t> is not a source node but a path
                            if len(seq2) > len(seq1):  # let's just "tuck" src1 into <t>
                                replace_node(src1, t)
                            else: # we don't know if nodes in <t> branch out so we can't collapse them
                                log.debug("should NOT collapse {0},{1}".format(src1, t))
                    else: # src1 and <t> are not similar enough, DO NOT collapse
                        log.debug("should NOT collapse {0},{1}".format(src1, t))

    

def find_bubbles(G, path_d, mermap):
    """
    We find all cases where n' -> n1 -> n3
                            n' -> n2 -> n3
    (that is, n3 has > 1 incoming) and n1, n2 each have only one incoming and one outgoing
    <i> make sure that n1 and n2 is not used in the same path  (which indicates in-gene repeat?)
    <ii> retrace n1, n2 to make sure that they are largely similar
    """
    def has_common_unique_pred(n1, n2):
        """
        Case:
         pred -> n1 -> common succ
         pred -> n2 -> common succ
        """
        preds1 = G.predecessors(n1)
        preds2 = G.predecessors(n2)
        return len(preds1) == 1  and len(preds2) == 1 and preds1[0] == preds2[0]

    def traceback_path(n1, n2):
        """
        Find a common pred where
         pred -> n1
         pred -> some_node -> n2
        """
        assert G.in_degree(n1) == 1
        pred = G.predecessors(n1)[0]
        return path_finder(G, n2, pred, [n2], 2)

    def replace_node(n_to_del, n_to_replace_with):
        #pdb.set_trace()
        G.remove_node(n_to_del)
        del mermap[n_to_del]
        for k in path_d:
            if n_to_del in path_d[k]:
                i = path_d[k].index(n_to_del)
                path_d[k] = path_d[k][:i] + [n_to_replace_with] + path_d[k][i+1:]

    def replace_path_w_node(path_to_del, n_to_replace_with, common_succ):
        """
        
        """
        # first, it's possible that the last node in <path_to_del> has other successors
        # ex: path_to_del = x1 -> x2 -> common_succ
        #          also has       x2 -> another node x3
        # so must change to n_to_replace_with -> x3
        last_n_in_path = path_to_del[-1]
        for s,t,data in G.out_edges(last_n_in_path, data=True):
            if G.has_edge(n_to_replace_with, t):
                G[n_to_replace_with][t]['weight'] += data['weight']
            else:
                G.add_edge(n_to_replace_with, t, weight=data['weight'])

        # for every predecssor of path_to_del, replace with n_to_replace_with
        # ex:        pred -> x1 -> x2 -> ...
        # becomes    pred -> n_to_replace_with -> ...
        for pred in G.predecessors(path_to_del[0]):
            G.add_edge(pred, n_to_replace_with, weight=G.get_edge_data(pred, path_to_del[0])['weight'])

        path_len = len(path_to_del)
        for k in path_d:
            if path_to_del[0] in path_d[k]:
                i = path_d[k].index(path_to_del[0])
                m = min(i+path_len, len(path_d[k]))
                if path_d[k][i:m] == path_to_del[:(m-i)]:
                    path_d[k] = path_d[k][:i] + [n_to_replace_with] + path_d[k][i+path_len:]
        # now delete all non branching nodes in path_to_del
        # note: this filter must be done simultaneously because G.remove_node will dynamically change the degrees!

        nodes_in_path = set()
        for path in path_d.itervalues():
            nodes_in_path = nodes_in_path.union(path)
        safe_to_remove = filter(lambda x: G.out_degree(x)<=1 and x not in nodes_in_path, path_to_del)
        #safe_to_remove = filter(lambda x: G.in_degree(x)<=1 and G.out_degree(x)<=1 and x not in nodes_in_path, path_to_del)
        for node in safe_to_remove:
            log.debug("safe to delete from G: {0}".format(node))
            G.remove_node(node)
            del mermap[node]
    

    in_same_path = make_in_same_path(path_d)
    cands = filter(lambda n: G.in_degree(n)>=2, G.nodes_iter())
    for n in cands:
        if n not in G: continue # deleted in loop below
        _pred = G.predecessors(n)
        if len(_pred) >= 2:
            for i, n1 in enumerate(_pred):
                if n1 not in G: continue
                for n2 in _pred[i+1:]:
                    if n1 not in G or n2 not in G or n1 in in_same_path[n2]: continue

                    if has_common_unique_pred(n1, n2):
                        # what is known: common pred -> n1 -> common succ
                        #                common pred -> n2 -> common succ
                        # so they must share the same first (KMER_SIZE-1) and the last (KMER_SIZE-1)
                        if DEBUG_FLAG:
                            pdb.set_trace()
                        if splice_align.node_is_similar(mermap[n1], mermap[n2]):
                            mermap[n1] = splice_align.get_consensus_through_voting(mermap[n1],\
                                                                                     G.get_edge_data(n1, n)['weight'],\
                                                                                     mermap[n2],\
                                                                                     G.get_edge_data(n2, n)['weight'])
                            replace_node(n_to_del=n2, n_to_replace_with=n1)
                        else:
                            flag, is_skipped = splice_align.node_is_skipping(mermap[n1], mermap[n2], cc_settings.KMER_SIZE)
                            if is_skipped:
                                if flag == "SEQ1":  # seq1 is the one with retained exon
                                    replace_node(n_to_del=n2, n_to_replace_with=n1)
                                else:
                                    replace_node(n_to_del=n1, n_to_replace_with=n2)
                            else:
                                log.debug("should NOT collapse {0},{1}".format(n1, n2))
                    else:
                        if G.in_degree(n1) == 1:
                            p2 = traceback_path(n1, n2)
                            if p2 is not None:
                                # common pred -> n1 -> common succ
                                # common pred -> another node -> n2 -> common succ
                                s1 = mermap[n1]
                                s2 = stitch_string_from_path(p2, mermap)
                                if DEBUG_FLAG:
                                    pdb.set_trace()
                                if splice_align.node_is_similar(s1, s2):
                                    mermap[n1] = splice_align.get_consensus_through_voting(s1,\
                                                                        G.get_edge_data(n1, n)['weight'],\
                                                                        s2,
                                                                        G.get_edge_data(n2, n)['weight'])
                                    replace_path_w_node(p2, n1, common_succ=n)

                                else:
                                    flag, is_skipped = splice_align.node_is_skipping(s1, s2, cc_settings.KMER_SIZE)
                                    if is_skipped:
                                        log.debug("path collapse possible {0},{1}".format(n1, p2))
                                        mermap[n1] = s1 if flag == 'SEQ1' else s2
                                        replace_path_w_node(p2, n1, common_succ=n)
                                    else:
                                        log.debug("should NOT collapse {0},{1}".format(n1, n2))
                            else:
                                log.debug("should NOT collapse {0},{1}".format(n1, n2))
                        elif G.in_degree(n2) == 1:
                            p1 = traceback_path(n2, n1)
                            if p1 is not None:
                                s1 = stitch_string_from_path(p1, mermap)
                                s2 = mermap[n2]
                                if DEBUG_FLAG:
                                    pdb.set_trace()
                                if splice_align.node_is_similar(s1, s2):
                                    mermap[n2] = splice_align.get_consensus_through_voting(s1,\
                                                                                             G.get_edge_data(n1, n)['weight'],\
                                                                                             s2,
                                                                                             G.get_edge_data(n2, n)['weight'])
                                    log.debug("path collapse possible: {0},{1}".format(p1, n2))
                                    replace_path_w_node(p1, n2, common_succ=n)
                                else:
                                    flag, is_skipped = splice_align.node_is_skipping(s1, s2, cc_settings.KMER_SIZE)
                                    if is_skipped:
                                        mermap[n2] = s1 if flag=='SEQ1' else s2
                                        replace_path_w_node(p1, n2, common_succ=n)
                                    else:
                                        log.debug("should NOT collapse {0},{1}".format(n1, n2))
                            else:
                                log.debug("should NOT collapse {0},{1}".format(n1, n2))
                        else:
                            log.debug("should NOT collapse {0},{1}".format(n1, n2))








