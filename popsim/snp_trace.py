# -*- coding: utf-8 -*-

import os


def snp_trace_path(dir_path, ind):
    return '{}.snp_trace'.format(os.path.join(dir_path, ind.id_str))


def check_snp_trace(dir_path, ind):
    return os.path.isfile(snp_trace_path(dir_path, ind))


def init_snp_trace(dir_path, ind, nb_snp):
    with open(snp_trace_path(dir_path, ind), 'w') as sf:
        for x in range(nb_snp):
            print(ind.fam_str, ind.fam_str, file=sf)


def get_snp_trace(dir_path, ind):
    with open(snp_trace_path(dir_path, ind), 'r') as sf:
        return [list(x) for x in zip(*[line.rstrip().split() for line in sf])]


def attach_snp_trace(chrome_pair, wd_path, ind):
    chrome_a, chrome_b = chrome_pair
    snp_trace_a, snp_trace_b = get_snp_trace(wd_path, ind)

    return [list(x) for x in [zip(chrome_a, snp_trace_a), zip(chrome_b, snp_trace_b)]]
