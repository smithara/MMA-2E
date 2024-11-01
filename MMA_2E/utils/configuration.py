#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:50:19 2024

@author: ngp
"""
"""
**Parameters**
   'max_gm_lat': 65.0,
   'min_gm_lat': 0.0,
   'LT_limit':   6.0,
    'n_max':     3,
    'm_max':     3,
    'ms':        [0, 1, -1, 0, 1, -1, 2, -2, 0, 1, -1, 2, -2, 3, -3],
    'ns':        [1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3],
    'dt':        8,
    'coordinates': 'GG' ,
    'q_file': os.path.join(LIB,'Q_08hours_kernels_1D_Grayver2017.h5'),
    'n_lag_days': 180, 
    'tfin':dt.now().date() - timedelta(days=11),
    'tini':dt.now().date() - timedelta(days=4)
 ======================  =============  =======================================
 Value                   Type           Description
 ======================  =============  =======================================
 'params.max_gm_lat'    `float`         Set the data selection to GM latitudes 
                                        in deg lower than this value.
 'params.min_gm_lat'    `float`         Set the data selection to GM latitudes 
                                        in deg larger than this value.
 'params.LT_limit'      `float`         Set the data selection to values in the 
                                        local time (LT) between LT_l in hours and
                                        24h-LT_l.
 'n_max':               `int`           Maximum degree used for the MMA model
 'm_max':               `int`           Maximum order  used for the MMA model
 'dt':                  `int`           Integration time for the SHA model in hours
 'params.r_surf'         `float`        Reference radius in kilometers
                                        (defaults to Earth's surface radius of
                                        6371.2 km).
 'params.r_cmb'          `float`        Core-mantle boundary radius in
                                        kilometers (defaults to 3485.0 km).
 'params.dipole'         `list`,        Coefficients of the dipole (used for
                         `ndarray`,     GSM/SM coordinate transformations) in
                         `shape (3,)`   units of nanotesla.
 'params.ellipsoid'      `list`,        Equatorial (index 0) and polar radius
                         `ndarray`      (index 1) of the spheroid describing
                         `shape (2,)`   Earth (WGS84) in units of kilometers.
 'params.CHAOS_version'  `str`          Version of the CHAOS model that was
                                        constructed with the built-in RC-index.
 'params.cdf_to_mjd'     `int`          Number of days on Jan 01, 2000 since
                                        Jan 01, 0000 (CDF start epoch)
 ======================  =============  =======================================

**Files**

 ==========================  ============  ====================================
 Value                       Type          Description
 ==========================  ============  ====================================
 'file.Q_Matrix_kernel'      `HDF5-file`, kernel to calculate the Q matrices
 ==========================  ============  ====================================

**Plots**

 ==========================  ===========  =====================================
 Value                       Type         Description
 ==========================  ===========  =====================================
 'plots.figure_width'        `float`      Plot width in inches (defaults to 6.3
                                          inches or 16 cm)
 ==========================  ===========  =====================================

.. autosummary::
    :toctree: classes
    :template: myclass.rst

    BasicConfig

"""


import os
# import re
import json
import numpy as np
import warnings
from contextlib import contextmanager
from datetime import timedelta as timedelta
from datetime import datetime as dt

ROOT = os.path.abspath(os.path.dirname(__file__))
LIB = os.path.join(ROOT, 'lib')

def check_path_exists(s):
    """Check that path to file exists."""
    if s is None or s == 'None':
        return None
    if os.path.exists(s):
        return s
    else:
        raise FileNotFoundError(f'{s} does not exist.')


def check_float(s):
    """Convert to float."""
    try:
        return float(s)
    except ValueError:
        raise ValueError(f'Could not convert {s} to float.')

def check_date(s):
    """Convert to date."""
    if isinstance(s,dt.date()):
        return s
    else:
        try:
            return dt.strptime(s,'%Y-%m-%d')
        except ValueError:
            raise ValueError(f'Could not convert {s} to date, try %Y-%m-%d format.')

def check_int(s):
    """Convert to integer."""
    try:
        return int(s)
    except ValueError:
        raise ValueError(f'Could not convert {s} to integer.')


def check_string(s):
    """Convert to string."""
    try:
        return str(s)
    except ValueError:
        raise ValueError(f'Could not convert {s} to string.')


def check_vector(s, len=None):
    """Check that input is vector of required length."""
    try:
        s = np.array(s)
        assert s.ndim == 1
        if len is not None:
            if s.size != len:
                raise ValueError(f'Wrong length: {s.size} != {len}.')
        return s
    except Exception as err:
        raise ValueError(f'Not a valid vector. {err}')


# def check_version_string(s):
#     """Check correct format of version string."""

#     s = check_string(s)

#     match = re.search(r'\d+\.\d+', s)
#     if match:
#         return s
#     else:
#         raise ValueError(f'Not supported version format "{s}".'
#                          'Must be of the form "x.x" with x an integer.')

def get_degOr(n_max,m_max):
    ms=[]
    ns=[]
    count=0
    for n in range(1,n_max+1):
        for m in range(min(m_max,n)+1):
            ms.append(m)
            ns.append(n)
            count=count+1
            if m!=0:
                ms.append(-m)
                ns.append(n)
                count=count+1
    return ns, ms

defparams = {
    'max_gm_lat': [65.0, check_float],
    'min_gm_lat': [0.0, check_float],
    'LT_limit':   [6.0, check_float],
    'n_max':      [3, check_int],
    'm_max':      [3, check_int],
    'ms':         [np.array([0, 1, -1, 0, 1, -1, 2, -2, 0, 1, -1, 2, -2, 3, -3]),
                   lambda x: check_vector(x, len=15)],
    'ns':        [np.array([1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3]),
                   lambda x: check_vector(x, len=15)],
                  'dt':        [8, check_int],
    'params.coordinates': ['GG', check_string],
    'q_file': [os.path.join(LIB,'Q_08hours_kernels_1D_Grayver2017.h5'),
                                check_path_exists],
    'n_lag_days': [180, check_int], 
    'tfin': [dt.now().date() - timedelta(days=11), check_date],
    'tini': [dt.now().date() - timedelta(days=4),check_date]
     }

class BasicConfig(dict):
    """Class for creating MMA 2E configuration dictionary."""

    defaults = defparams

    def __init__(self, *args, **kwargs):
        super().update(*args, **kwargs)
        wrap = lambda v: BasicConfig(v) if type(v) is dict else v
        kvdict = {k: wrap(v) for k, v in dict(*args, **kwargs).items()}
        super(BasicConfig, self).__init__(kvdict)
        self.__dict__ = self
        
    def __setitem__(self, key, value):
        """Set and check value before updating dictionary."""

        try:
            try:
                cval = self.defaults[key][1](value)
            except ValueError as err:
                raise ValueError(f'Key "{key}": {err}')
            super().__setitem__(key, cval)
        except KeyError:
            raise KeyError(f'"{key}" is not a valid parameter.')

    def __str__(self):
        return '\n'.join(map('{0[0]}: {0[1]}'.format, sorted(self.items())))

    def reset(self, key):
        """
        Load default values.

        Parameters
        ----------
        key : str
            Single keyword that is reset to the default.

        """
        self.__setitem__(key, self.defaults[key][0])

    def fullreset(self):
        """
        Load all default values.

        """
        super().update({key: val for key, (val, _) in self.defaults.items()})

    def load(self, filepath):
        """
        Load configuration dictionary from file.

        Parameters
        ----------
        filepath : str
            Filepath and name to json-formatted configuration txt-file.

        """

        with open(filepath, 'r') as f:
            kwargs = json.load(f)

        if len(kwargs) == 0:
            warnings.warn(
                'Configuration dictionary loaded from file is empty.')

        for key, value in kwargs.items():
            # check format and set key value pairs
            self.__setitem__(key, value)

    def save(self, filepath):
        """
        Save configuration dictionary to a file.

        Parameters
        ----------
        filepath : str
            Filepath and name of the textfile that will be saved with the
            configuration values.

        """

        def default(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()

        with open(filepath, 'w') as f:
            json.dump(self, f, default=default, indent=4, sort_keys=True)

        print(f'Saved configuration textfile to {filepath}.')

    @contextmanager
    def context(self, key, value):
        """
        Use context manager to temporarily change setting.

        Parameters
        ----------
        key : str
            BasicConfig configuration key.
        value
            Value compatible with ``key``.

        Examples
        --------
        Temporarily change the radius of Earth's surface for a computation
        and then change it back to the original value.

        .. code-block:: python

          from chaosmagpy import basicConfig

          print('Before: ', basicConfig['params.r_surf'])

          # change Earth's radius to 10 km
          with basicConfig.context('params.r_surf', 10):
              # do something at r_surf = 10 km ...
              print('Inside: ', basicConfig['params.r_surf'])

          print('After: ', basicConfig['params.r_surf'])

        """
        old_value = self.__getitem__(key)
        self.__setitem__(key, value)
        yield
        self.__setitem__(key, old_value)

