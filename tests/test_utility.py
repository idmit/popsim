# -*- coding: utf-8 -*-

import pytest

import popsim.utility as ps


@pytest.fixture
def ordered_list():
    return list(range(4))


def test_swap(ordered_list):
    assert ps.swap(ordered_list, 0, 2) == [2, 1, 0, 3]


def test_repeated_swap(ordered_list):
    assert ps.swap(ps.swap(ordered_list, 0, 2), 0, 3) == [3, 1, 2, 0]
