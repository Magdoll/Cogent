__author__ = 'etseng@pacb.com'

'''
Drawing functions using matplotlib for visualizing k-mer family finding alg.
Not core function. Not tested. Used by developer only.
'''

import matplotlib.pyplot as plt

def draw_reduced_graph(G):
    plt.figure(figsize=(14, 14))
    plt.clf()
    pos = nx.spring_layout(G)

    # generate N points along the same line

def draw_assignment(G, pos, labels, labels2, answer_d, output_prefix):
    plt.figure(figsize=(14, 14))
    plt.clf()
    groups = np.unique(labels2)
    m = len(groups)
    for i in range(m):
        nx.draw_networkx_nodes(G, pos, nodelist=list(labels[labels2==groups[i]]), \
                node_color=plt.cm.jet(i*1./m), node_size=300, alpha=.5)

    ans_labels = {}
    for k,v in answer_d.items():
        for x in v: ans_labels[x] = k
    nx.draw_networkx_labels(G, pos, ans_labels, font_size=12)

    plt.axis('off')
    plt.savefig(output_prefix +'.png')
    plt.show()

    # show only edges that are > 0.9
    pp = [a_b_d for a_b_d in G.edges(data=True) if a_b_d[2]['weight']>=.9]
    nx.draw_networkx_edges(G, pos, edgelist=pp, width=1, alpha=.5)
    plt.savefig(output_prefix+'.with_edges.png')
    plt.show()

def assign_pos_by_answer(G, pos, answer_d):
    m = math.sqrt(len(answer_d)) + 1
    keys = list(answer_d.keys())
    for i in range(len(keys)):
        _row, _col = i/m, i%m
        x_min = _row * 1. / m
        x_max = x_min + 1./m
        y_min = _col * 1. / m
        y_max = y_min + 1./m
        pts = generate_points_in_limit(len(answer_d[keys[i]]), x_min, x_max, y_min, y_max)
        for j,v in enumerate(answer_d[keys[i]]): pos[v] = pts[j,:]

def generate_points_in_limit(n, x_min, x_max, y_min, y_max):
    p = np.random.random((n*2, 2))
    func = np.bitwise_and
    good = p[func(func(x_min <= p[:,0], p[:,0] < x_max), func(y_min <= p[:,1], p[:,1] < y_max))]
    if len(good) >= n: return good
    while True:
        p = np.random.random((n, 2))
        tmp = p[func(func(x_min <= p[:,0], p[:,0] < x_max), func(y_min <= p[:,1], p[:,1] < y_max))]
        good = np.vstack([good, tmp])
        if len(good) >= n: break
    return good[:n, :]


def plot_proportion_related_vs_unrelated(seq_d, output_filename):
    f = open(output_filename, 'w')
    f.write("id1\tid2\ttag\tsim\n")

    n = len(seq_d)
    keys = list(seq_d.keys())
    for i in range(n-1):
        x1 = seq_d[keys[i]]
        for j in range(i+1, n):
            x2 = seq_d[keys[j]]
            sim = len(x1.intersection(x2))*1. / min(len(x1), len(x2))
            if keys[i].split('.')[1] == keys[j].split('.')[1]: # ex: PB.1.3 vs PB.1.4
                tag = 'related'
            else:
                tag = 'unrelated'
            f.write("{0}\t{1}\t{2}\t{3}\n".format(keys[i], keys[j], tag, sim))
    f.close()
