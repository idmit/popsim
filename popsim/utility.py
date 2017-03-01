# -*- coding: utf-8 -*-


def swap(lst, a, b):
    new_lst = list(lst)

    if new_lst[a] != b:
        if new_lst[a] != a:
            new_lst[new_lst[a]], new_lst[a] = new_lst[a], new_lst[new_lst[a]]

        if new_lst[b] != b:
            new_lst[new_lst[b]], new_lst[b] = new_lst[b], new_lst[new_lst[b]]

    new_lst[a], new_lst[b] = new_lst[b], new_lst[a]

    return new_lst
