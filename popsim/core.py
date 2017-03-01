# -*- coding: utf-8 -*-

import click
import numpy as np
import os
import popsim.plink as plink_mod
from itertools import tee


class Individual:
    def __init__(self, tfam_idx, fam_str, id_str):
        self.tfam_idx = tfam_idx
        self.fam_str = fam_str
        self.id_str = id_str
        self.chrome_pair = None


def individual_from_tfam(tfam_lines, id_str):
    for idx, line in enumerate(tfam_lines):
        cols = line.rstrip().split()

        if cols[1] == id_str:
            return Individual(idx, cols[0], cols[1])


def centimorgans_from_bp(bp_dist):
    return bp_dist / 1000000


def swap_p(cmorgan_dist):
    return (1 - np.exp(-2 * cmorgan_dist / 100)) / 2


def extract_markers_info(tped_lines, nb_markers):
    markers_info = [None for _ in range(nb_markers)]

    for idx, line in enumerate(tped_lines):
        cols = line.rstrip().split()

        markers_info[idx] = (cols[1], int(cols[3]))

    return markers_info


def untracked_chrome_pairs_from_tped(tped_lines, tfam_indices, nb_markers):
    untracked_chrome_pairs = [([None for _ in range(nb_markers)], [None for _ in range(nb_markers)])
                              for _ in tfam_indices]

    for idx, line in enumerate(tped_lines):
        gtypes = line.rstrip().split()[4:]
        for ind_idx, tfam_idx in enumerate(tfam_indices):
            untracked_chrome_pairs[ind_idx][0][idx] = gtypes[2 * tfam_idx]
            untracked_chrome_pairs[ind_idx][1][idx] = gtypes[2 * tfam_idx + 1]

    return untracked_chrome_pairs


def reorder_tuples(order_iterable, tuple_iterable):
    return [[a_tuple[idx] for idx in order]
            for order, a_tuple in zip(order_iterable, tuple_iterable)]


def reorder_chrome_pair(order_iterable, chrome_pair):
    (a, b), (c, d) = map(tee, chrome_pair)
    return reorder_tuples(order_iterable, zip(a, b, c, d))


def pairwise(iterable):
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def gamete_set_from_chrome_pair(markers_info, chrome_pair):
    import itertools
    import popsim.utility as util
    pairs = [(0, 2), (0, 3), (1, 2), (1, 3)]

    order_iterable = [None for _ in markers_info]
    order_iterable[0] = list(range(4))
    order_changes = []

    for idx, (snp, next_snp) in enumerate(pairwise(markers_info)):
        cmorgan_dist = centimorgans_from_bp(next_snp[1] - snp[1])
        p_of_swap = swap_p(cmorgan_dist)

        if np.random.rand() <= p_of_swap:
            choice_idx = np.random.choice(4)
            order_iterable[idx + 1] = util.swap(order_iterable[idx], *pairs[choice_idx])
            order_changes.append(idx)
        else:
            order_iterable[idx + 1] = order_iterable[idx]

    reordered = reorder_chrome_pair(order_iterable, chrome_pair)

    print('Order changed {} times at indices: {}.'.format(len(order_changes), order_changes))
    print('Consecutive order values: {}'.format([k for k, g in itertools.groupby(order_iterable)]))
    # print('Pair of chromes', 'Order', 'Gametes')
    # for a, b, c in zip(zip(*chrome_pair), order_iterable, reordered):
    #     print(a, b, c)

    return list(zip(*reordered))


def chrome_pair_from_gamete_sets(inds):
    gamete_indices = [np.random.choice(4) for _ in inds]

    print('Picked gametes with indices {}'.format(gamete_indices))

    chrome_pair = [ind.gamete_set[gamete_idx] for ind, gamete_idx in zip(inds, gamete_indices)]

    distinct_families_pair = [set(list(zip(*chrome))[1]) for chrome in chrome_pair]

    if any([len(dfs) - 1 for dfs in distinct_families_pair]):
        print('-- Got yourself a MIX! --')

    print('The child has {} and {} mixed in.'.format(*distinct_families_pair))

    return chrome_pair


def child_from_parents(tpd_stem, child_id_str, parent_id_strings):
    import popsim.snp_trace as st

    wd_path = os.path.dirname(tpd_stem)
    tped_path = '{}.tped'.format(tpd_stem)
    tfam_path = '{}.tfam'.format(tpd_stem)

    with open(tped_path) as tped_f:
        nb_markers = sum(1 for _ in tped_f)

    with open(tped_path) as tped_f:
        markers_info = extract_markers_info(tped_f, nb_markers)

    ind_pair = [None for _ in parent_id_strings]

    for idx, id_str in enumerate(parent_id_strings):
        with open(tfam_path) as tfam_f:
            ind_pair[idx] = individual_from_tfam(tfam_f, id_str)

    for ind in ind_pair:
        if not st.check_snp_trace(wd_path, ind):
            st.init_snp_trace(wd_path, ind, nb_markers)

    with open(tped_path) as tped_f:
        untracked_chrome_pairs = untracked_chrome_pairs_from_tped(tped_f,
                                                                  [ind.tfam_idx for ind in ind_pair],
                                                                  nb_markers)

    for ind_idx, ind in enumerate(ind_pair):
        ind.chrome_pair = st.attach_snp_trace(untracked_chrome_pairs[ind_idx], wd_path, ind)
        ind.gamete_set = gamete_set_from_chrome_pair(markers_info, ind.chrome_pair)
        if child_id_str == 'SIM_chrome_1_gen_1_1':
            print(len(ind.gamete_set))
            for alleles in zip(*ind.gamete_set):
                print([t[0] for t in alleles])

    new_chrome_pair = chrome_pair_from_gamete_sets(ind_pair)

    child = Individual(0, child_id_str, child_id_str)
    child.chrome_pair = new_chrome_pair

    return child


def replace_random_parents(parent_id_strings, id_strings, mom, dad):
    for idx, parent_id_str in enumerate([mom, dad]):
        if parent_id_str and parent_id_str in id_strings:
            parent_id_strings[idx] = parent_id_str
            print('Overriding a previously picked parent with {}.'.format(parent_id_str))
        else:
            print('Proceeding with a random parent.')


@click.command()
@click.argument('stem', type=click.Path(readable=True, resolve_path=True))
@click.argument('chrome_idx', type=click.INT)
@click.argument('gen_idx', type=click.INT)
@click.argument('nb_children', type=click.INT)
@click.option('--seed', type=click.INT, help='Pass the seed to produce deterministic results.')
@click.option('--mom', type=click.STRING, help='Pass the first parent id string.')
@click.option('--dad', type=click.STRING, help='Pass the second parent id string.')
@click.option('--plink', type=click.STRING, help='Path to plink binary.')
def root(stem, chrome_idx, gen_idx, nb_children, seed, mom, dad, plink):
    from collections import OrderedDict
    from itertools import chain

    if seed is not None:
        np.random.seed(seed)

    plink_path = './plink'

    if plink:
        plink_path = plink

    plink_mod.ensure_tpd(plink_path, stem)
    tfam_path = '{}.tfam'.format(stem)

    merge_list_path = '{}/SIM_chrome_{}_gen_{}.merge_list'.format(os.path.dirname(stem), chrome_idx, gen_idx)

    with open(tfam_path) as tfam_f, open(merge_list_path, 'w') as merge_f:
        id_strings_by_fam_str = OrderedDict()

        for line in tfam_f:
            cols = line.rstrip().split()
            if cols[0] not in id_strings_by_fam_str:
                id_strings_by_fam_str[cols[0]] = []

            id_strings_by_fam_str[cols[0]].append(cols[1])

        fam_strings = list(id_strings_by_fam_str.keys())
        id_strings = list(chain(*id_strings_by_fam_str.values()))

        for idx in range(nb_children):
            child_id_str = 'SIM_chrome_{}_gen_{}_{}'.format(chrome_idx, gen_idx, idx + 1)

            print('---')
            print('| Working on creating {}.'.format(child_id_str))
            print('---')

            parent_fam_strings = [None, None]
            parent_fam_strings[0] = np.random.choice(fam_strings)
            parent_fam_strings[1] = parent_fam_strings[0]

            while parent_fam_strings[1] == parent_fam_strings[0]:
                parent_fam_strings[1] = np.random.choice(fam_strings)

            parent_id_strings = [np.random.choice(id_strings_by_fam_str[fam_str]) for fam_str in parent_fam_strings]

            for parent_id_str, parent_fam_str in zip(parent_id_strings, parent_fam_strings):
                print('Picked {} from {} as a parent.'.format(parent_id_str, parent_fam_str))

            replace_random_parents(parent_id_strings, id_strings, mom, dad)

            child = child_from_parents(stem, child_id_str, parent_id_strings)

            child_stem = plink_mod.save_ped(os.path.dirname(stem), child)
            plink_mod.copy_bim_as_map(stem, child_stem)
            plink_mod.pd_to_bpd(plink_path, child_stem, child_stem)

            print('{}.bed'.format(child_id_str), '{}.bim'.format(child_id_str), '{}.fam'.format(child_id_str),
                  file=merge_f)

    plink_mod.merge(plink_path, stem, merge_list_path, '{}_with_gen_{}'.format(stem, gen_idx))
