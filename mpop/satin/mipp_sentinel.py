#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c14526.ad.smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Sentinel-1 reader
"""

import ConfigParser
import os

from mipp.read_geotiff import read_geotiff
from datetime import datetime
import numpy as np

from glob import glob
import mpop
from mpop import CONFIG_PATH
import logging

LOG = logging.getLogger(__name__)

from mpop.plugin_base import Reader

class SentinelGRDChannel(mpop.channel.GenericChannel):

    def __init__(self, name='unknown'):
        mpop.channel.GenericChannel.__init__(self)
        self._is_loaded = False
        self.name = name
        self.data = None
        self.shape = None
        self._projectables = []
        self._projectables.append(name)

    def is_loaded(self):
        return self._is_loaded

    def set_loaded(self):
        self._is_loaded = not self._is_loaded

    def project(self, coverage):
        """Project what can be projected in the product.
        """
        
        import copy
        res = copy.copy(self)
        res.data = coverage.project_array(self.data)
        return res

class GeoTiffReader(Reader):

    pformat = "mipp_sentinel"

    def load(self, satscene, **kwargs):

        LOG.debug('channels to load: ' + str(satscene.channels_to_load))
        conf = ConfigParser.ConfigParser()
        conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))
        options = {}

        for option, value in conf.items(satscene.instrument_name + "-level2",
                                        raw=True):
            options[option] = value

        options["resolution"] = kwargs.get("resolution", 'unknown')
        values = {"orbit": satscene.orbit}

        path_template = datetime.strftime(satscene.time_slot, options["dir"]) % values

        dirlist = glob(path_template)
        
        if len(dirlist) != 1:
            raise IOError("Couldn't identify a unique measurments directory!")

        dirname = dirlist[0]
        if not os.path.exists(dirname):
            raise IOError('Directory ' + str(dirname) + ' does not exist')

        filelist = glob(os.path.join(dirname, options["filename"]))
        if len(filelist) == 0:
            LOG.warning('No files found!')
            
        lons = np.array([])
        lats = np.array([])
        data = np.array([])
        # Concatenation of granules
        for granule_file in filelist:
            LOG.debug('Load granule: ' + str(granule_file))
            local_lons, local_lats, local_data = self.load_granule(granule_file)
            lons = np.concatenate([lons, local_lons.ravel()])
            lats = np.concatenate([lats, local_lats.ravel()])
            data = np.concatenate([data, local_data.ravel()])
    
        # fill the scene object
        from pyresample import geometry
        satscene.area = geometry.SwathDefinition(lons=lons, lats=lats)

        # Take only the fist - assume only one channel!
        # FIXME!
        for chname in satscene.channels_to_load:
            chn = SentinelGRDChannel(chname)
            chn.set_loaded()
            chn.data = np.ma.masked_array(data)
            chn.shape = data.shape
            satscene[chname] = chn
            break


    def load_granule(self, filename):
        """Load one sentinel granule file"""
        from geotiepoints import SatelliteInterpolator

        params, data = read_geotiff(filename)
        tie_lons = params['tiepoints']['lons']
        tie_lats = params['tiepoints']['lats']
        tie_cols = params['tiepoints']['cols']
        tie_rows = params['tiepoints']['rows']

        fine_cols = np.arange(0, data.shape[1])
        fine_rows = np.arange(0, data.shape[0])
        interpolator = SatelliteInterpolator((tie_lons, tie_lats),
                                             (tie_rows, tie_cols),
                                             (fine_rows, fine_cols),
                                             1, 3)
        lons, lats = interpolator.interpolate()
        
        return lons, lats, data
