# -*- coding: utf-8 -*-

import popsim.core as ps

import pytest


@pytest.fixture
def tped_lines():
    return ['1 Hapmap43437-BTA-101873 0 135098 1 1 2 2',
            '1 ARS-BFGL-NGS-16466 0 267940 4 4 2 1',
            '1 Hapmap34944-BES1_Contig627_1906 0 393248 3 3 2 1',
            '1 ARS-BFGL-NGS-98142 0 471078 4 3 3 3']


def test_chrome_markers_gtr(tped_lines):
    assert ps.extract_markers_info(tped_lines, 4) == [('Hapmap43437-BTA-101873', 135098),
                                                      ('ARS-BFGL-NGS-16466', 267940),
                                                      ('Hapmap34944-BES1_Contig627_1906', 393248),
                                                      ('ARS-BFGL-NGS-98142', 471078)]


def test_chrome_markers_gtr_throws():
    wrong_tped_lines = ['1 Hapmap43437-BTA-101873 0 NOT_INT 1 1 2 2']

    with pytest.raises(ValueError):
        list(ps.extract_markers_info(wrong_tped_lines, 1))


def test_chrome_gtypes(tped_lines):
    assert ps.untracked_chrome_pairs_from_tped(tped_lines, [0, 1], 4) == [(['1', '4', '3', '4'], ['1', '4', '3', '3']),
                                                                          (['2', '2', '2', '3'], ['2', '1', '1', '3'])]


def test_reorder_tuples():
    from collections import deque

    x = iter(range(5))
    y = ['a', 'b', 'c', 'd', 'e']
    z = ['z', 'y', 'x', 'w', 'v']

    order_tuples = [deque(range(3)) for _ in range(5)]
    for i, r in enumerate(order_tuples):
        r.rotate(i)

    assert ps.reorder_tuples(order_tuples, zip(x, y, z)) == [[0, 'a', 'z'],
                                                             ['y', 1, 'b'],
                                                             ['c', 'x', 2],
                                                             [3, 'd', 'w'],
                                                             ['v', 4, 'e']]


def test_reorder_chrome_pair():
    from collections import deque

    chrome_pair = (['a', 'b', 'c', 'd', 'e'],
                   ['z', 'y', 'x', 'w', 'v'])

    order_tuples = [deque(range(4)) for _ in range(5)]
    for i, r in enumerate(order_tuples):
        r.rotate(i)

    assert ps.reorder_chrome_pair(order_tuples, chrome_pair) == [['a', 'a', 'z', 'z'],
                                                                 ['y', 'b', 'b', 'y'],
                                                                 ['x', 'x', 'c', 'c'],
                                                                 ['d', 'w', 'w', 'd'],
                                                                 ['e', 'e', 'v', 'v']]
