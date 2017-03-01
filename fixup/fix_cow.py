#!/usr/bin/env python

import sys


def pairs_from_rmap(rmap_f):
    for line in rmap_f:
        s = line.rstrip().split()

        snp_id, *chrom_and_pos = s[:3]

        yield snp_id, chrom_and_pos


if __name__ == '__main__':
    bim_path = '../data/cow.bim'
    snp_list_path = '../data/positioned_snp.list'
    nbim_path = '../data/fixed_cow.bim'
    rmap_path = '../data/cattle_rmap.txt'

    with open(rmap_path, 'r') as rmap_f:
        rmap_dict = {
            snp_id: chrom_and_pos
            for snp_id, chrom_and_pos in pairs_from_rmap(rmap_f)
        }

    with_info_snp_list = set()

    with open(bim_path, 'r') as bim_f, open(nbim_path, 'w') as nbim_f:
        for line in bim_f:
            s = line.rstrip().split()

            chrom, snp_id, cmorgans, bp_pos, *gtypes = s

            with_info_snp_list.add(snp_id)

            if snp_id in rmap_dict:
                nchrom, nbp_pos = rmap_dict[snp_id]
                del rmap_dict[snp_id]
            elif 'ARS-{}'.format(snp_id) in rmap_dict:
                nchrom, nbp_pos = rmap_dict['ARS-{}'.format(snp_id)]
                del rmap_dict['ARS-{}'.format(snp_id)]
            else:
                with_info_snp_list.remove(snp_id)
                continue

            print(
                nchrom,
                snp_id,
                cmorgans,
                nbp_pos,
                *gtypes,
                sep='\t',
                file=nbim_f)

    with open(snp_list_path, 'w') as snp_list_f:
        for snp_id in with_info_snp_list:
            print(snp_id, file=snp_list_f)
