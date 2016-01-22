import logging
from collections import defaultdict
from Cogent import all_simple_paths
from Cogent.settings import KMER_SIZE
from pulp import LpProblem, LpMinimize, LpVariable, LpInteger

log = logging.getLogger('Cogent.process_path')


def make_in_same_path(path_d):
    """
    go through path_d which is dict of seqid --> list of path nodes (in G)
    return a dict of node id --> set of nodes in the same path
    """
    in_same_path = defaultdict(lambda: set())
    for path in path_d.itervalues():
        for i,n1 in enumerate(path):
            for n2 in path[i+1:]:
                in_same_path[n1].add(n2)
                in_same_path[n2].add(n1)
    return in_same_path


def stitch_string_from_path(path, mermap):
    """
    mermap --- dict of <node> --> <sequence>
    path --- list of nodes ex: [0, 71, 100]
    """
    seq = mermap[path[0]]
    for i in xrange(1, len(path)):
        seq += mermap[path[i]][(KMER_SIZE-1):]
    return seq


# CURRENTLY UNUSED, for debugging only
# def count_subpath_freq(path_d, path):
#     count = 0
#     for k,p in path_d.iteritems():
#         try:
#             i = p.index(path[0])
#             if all(a==b for (a,b) in itertools.izip(p[i+1:], path[1:])):
#                 count += 1
#         except ValueError:
#             pass
#     return count

# CURRENTLY UNUSED, for debugging only
# def greedy_path_walk(G, sink, path=[]):
#     """
#     Start from <sink> node, continuously walk down until hits and end.
#     Whenever there are multiple successors always choose the one with higher weight.
#     Return: (nodes_in_path)
#     """
#     if G.out_degree(sink) == 0:
#         yield path
#     else:
#         it = G.out_edges_iter(sink, data=True)
#         for x in it:
#             for p in greedy_path_walk(G, x[1], path+[x[1]]): yield p



def path_finder(G, n, target, path, move_left):
    if G.in_degree(n) == 0 or move_left == 0: return None # failed, did not find it
    preds = G.predecessors(n)
    for pred in preds:
        if pred == target: # found it!
            return path
        else: # keep trying
            p = path_finder(G, pred, target, [pred]+path, move_left-1)
            if p is not None: return p
    return None


def path_match(target_path, match_path):
    if len(target_path) == 0: return True
    if target_path[0] in match_path:
        i = match_path.index(target_path[0])
        return path_match(target_path[1:], match_path[i+1:])
    else:
        return False


def make_into_lp_problem(good_for, N):
    """
    Helper function for solve_with_lp_and_reduce()

    N --- number of isoform sequences
    good_for --- dict of <isoform_index> --> list of matched paths index
    """
    prob = LpProblem("The Whiskas Problem",LpMinimize)

    # each good_for is (isoform_index, [list of matched paths index])
    # ex: (0, [1,2,4])
    # ex: (3, [2,5])
    variables = [LpVariable(str(i),0,1,LpInteger) for i in xrange(N)]

    # objective is to minimize sum_{Xi}
    prob += sum(v for v in variables)

    # constraints are for each isoform, expressed as c_i * x_i >= 1
    # where c_i = 1 if x_i is matched for the isoform
    # ex: (0, [1,2,4]) becomes t_0 = x_1 + x_2 + x_4 >= 1
    for t_i, p_i_s in good_for:
        #c_i_s = [1 if i in p_i_s else 0 for i in xrange(N)]
        prob += sum(variables[i]*(1 if i in p_i_s else 0) for i in xrange(N)) >= 1
    prob.writeLP('aloha.lp')
    return prob


def solve_with_lp_and_reduce(good_for, paths, mermap, outfile='aloha.fa'):
    prob = make_into_lp_problem(good_for, len(paths))
    prob.solve()
    for v in prob.variables(): log.debug("{0} = {1}".format(v.name, v.varValue))
    touse = []
    for v in prob.variables():
        if v.varValue == 1: touse.append(int(v.name))

    #HACK
    #touse = range(len(paths))

    with open(outfile, 'w') as f:
        for u in touse:
            seq = stitch_string_from_path(paths[u], mermap)
            f.write(">path{0}\n{1}\n".format(u, seq))


def find_minimal_path_needed_to_explain_pathd(G, path_d, keys, max_G_size=50):
    used_path = set()

    if G.number_of_nodes() <= max_G_size:
        sources = filter(lambda n: G.in_degree(n) == 0, G.nodes_iter())
        sinks = filter(lambda n: G.out_degree(n) == 0, G.nodes_iter())
        paths = []
        for src in sources:
            for sink in sinks:
                if src == sink: paths += [[src]]
                else:
                    paths += [p for p in all_simple_paths.all_simple_paths(G, src, sink)]
                log.info("number of paths now: {0}".format(len(paths)))
    else:
        log.info("Number of nodes exceeds {0}. Not dealing with pathological case now. Just use the paths!".format(max_G_size))
        paths = path_d.values()

    good_for = defaultdict(lambda: [])
    for j,k in enumerate(keys):
        kk = path_d[k]
        for i,p in enumerate(paths):
            if path_match(kk, p):
                good_for[j].append(i)
                used_path.add(i)

    for j,k in enumerate(keys):
        if j not in good_for:
            log.warning("[BUG] Missing good_for for {0}. Maybe known issue (graph cycles). Use own path.".format(k))
            i = len(paths)
            paths.append(path_d[k])
            good_for[j] = [i]
            used_path.add(i)

    # trim paths down to the used paths
    used_path = list(used_path)
    path_map = dict((_old,_new) for _new,_old in enumerate(used_path))# old path index -> new path index
    for k in good_for:
        good_for[k] = [path_map[_old] for _old in good_for[k]]
    new_paths = {}
    for _new,_old in enumerate(used_path):
        new_paths[_new] = paths[_old]
    return good_for.items(), new_paths