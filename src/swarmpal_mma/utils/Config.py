#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 10:48:21 2024

@author: ngp
"""

"""
Configuration dictionary
"""

import os
import json
from datetime import timedelta as timedelta
from datetime import datetime as dt

ROOT = os.path.abspath(os.path.dirname(__file__))
LIB = os.path.join(ROOT, 'kernels') 
DEFAULTS = {
    'params.max_gm_lat': 60.0,
    'params.min_gm_lat': 0,
    'params.LT_limit':   6.0,
    'params.n_max':      3,
    'm_max':     3,
    'ms':        [0, 1, -1, 0, 1, -1, 2, -2, 0, 1, -1, 2, -2, 3, -3],
    'ns':        [1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3],
    'params.dt':        8,
    'params.coordinates': 'GG' ,
    'file.q_file': os.path.join(LIB,'Q_08hours_kernels_1D_Grayver2017.h5'),
    'params.n_lag_days': 180, 
    'params.tfin':dt.now().date() - timedelta(days=11),
    'params.tini':dt.now().date() - timedelta(days=4)
     }


class Config(dict):
    """
    Store configuration data in dictionary and access via dot notation.

    >>> contact = Config({"address": {"street": "Baker St", "number": 221}})

    >>> print(contact['address']['number'])  # the traditional way
    221

    >>> print(contact.address.number)  # much nicer
    221

    >>> print(contact.address)  # pretty printing of content
    {
      "number": 221,
      "street": "Baker St"
    }

    >>> contact.address.number = 111  # update of values
    >>> print(contact.address)
    {
      "number": 111,
      "street": "Baker St"
    }

    >>> contact.address = Config({"street": "Downing", "number": 12})
    >>> print(contact.address.street)
    Downing
    """

    def __init__(self, *args, **kwargs):
        """ Construct the same way a plain Python dict is created. """
        wrap = lambda v: Config(v) if type(v) is dict else v
        kvdict = {k: wrap(v) for k, v in dict(*args, **kwargs).items()}
        super(Config, self).__init__(kvdict)
        self.__dict__ = self

    @staticmethod
    def load(filepath):
        """Load configuration from given JSON filepath"""
        with open(filepath) as f:
            return Config(json.load(f))

    def save(self, filepath):
        """Save configuration to given JSON filepath"""
        with open(filepath, 'w') as f:
            json.dump(self, f, indent=2)

    def __repr__(self):
        """Pretty string representation of configuration"""
        return json.dumps(self, sort_keys=True, indent=2)