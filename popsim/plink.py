# -*- coding: utf-8 -*-

import os
import subprocess as sp

tpd_exts = ['tped', 'tfam']


class NonZeroReturnCodeError(Exception):
    def __init__(self, info):
        super(NonZeroReturnCodeError, self).__init__(info['bin_name'], info['stdout'], info['stderr'])

        print('{} has finished with a non-zero error code: {}'.format(info['bin_name'], info['return_code']))
        print(info['stdout'])
        print(info['stderr'])


def files_exist(stem, exts):
    return all(os.path.isfile('{}.{}'.format(stem, ext)) for ext in exts)


def ensure_tpd(plink_path, bpd_stem):
    if files_exist(bpd_stem, tpd_exts):
        return

    process = sp.Popen([plink_path, '--bfile', bpd_stem, '--out', bpd_stem, '--recode', '--transpose', '--noweb'],
                       cwd=os.path.dirname(bpd_stem), stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise NonZeroReturnCodeError(
            {'bin_name': 'plink', 'stdout': stdout, 'stderr': stderr, 'return_code': process.returncode})


def copy_bim_as_map(src_stem, dst_stem):
    with open('{}.bim'.format(src_stem)) as bim_f, open('{}.map'.format(dst_stem), 'w') as map_f:
        for line in bim_f:
            print('\t'.join(line.rstrip().split('\t')[:-2]), file=map_f)


def save_ped(dir_path, ind):
    ped_stem = os.path.join(dir_path, ind.id_str)

    with open('{}.ped'.format(ped_stem), 'w') as ped_f, open('{}.snp_trace'.format(ped_stem), 'w') as st_f:
        print(ind.id_str, ind.id_str, 0, 0, 0, -9, file=ped_f, end=' ')

        for (allele_a, snp_trace_a), (allele_b, snp_trace_b) in zip(*ind.chrome_pair):
            print(allele_a, allele_b, file=ped_f, end=' ')
            print(snp_trace_a, snp_trace_b, file=st_f)

    return ped_stem


def pd_to_bpd(plink_path, src_stem, dst_stem):
    process = sp.Popen([plink_path, '--file', src_stem, '--out', dst_stem, '--make-bed', '--noweb'],
                       cwd=os.path.dirname(src_stem), stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise NonZeroReturnCodeError(
            {'bin_name': 'plink', 'stdout': stdout, 'stderr': stderr, 'return_code': process.returncode})


def merge(plink_path, bpd_stem, merge_list_path, out_stem):
    process = sp.Popen(
        [plink_path, '--bfile', bpd_stem, '--merge-list', merge_list_path, '--out', out_stem, '--make-bed', '--noweb'],
        cwd=os.path.dirname(bpd_stem), stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise NonZeroReturnCodeError(
            {'bin_name': 'plink', 'stdout': stdout, 'stderr': stderr, 'return_code': process.returncode})
