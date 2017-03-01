# -*- coding: utf-8 -*-

import pytest

from popsim.core import takeoff


def test_takeoff():
    assert takeoff() == 'Takeoff complete!'
