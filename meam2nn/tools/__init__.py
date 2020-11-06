# coding: utf-8

from .atomic_info import atomic_number

def aslist(val):
    if isinstance(val, str):
        return [val]
    else:
        return list(val)