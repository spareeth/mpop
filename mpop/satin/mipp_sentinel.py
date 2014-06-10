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

from mipp import read_geotiff
from datetime import datetime

from glob import glob
from mpop import CONFIG_PATH
import logging

LOG = logging.getLogger(__name__)

from mpop.plugin_base import Reader

class GeoTiffReader(Reader):

    pformat = "mipp_sentinel"
    
    def load(self, satscene, **kwargs):

        conf = ConfigParser.ConfigParser()
        conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))
        options = {}

        for option, value in conf.items(satscene.instrument_name + "-level2",
                                        raw=True):
            options[option] = value

        options["resolution"] = kwargs.get("resolution", 'unknown')
        values = {"orbit": satscene.orbit,
                  "mode": '*',
                  "product_type": '*',
                  "resolution": '*'
                  }

        path_template = datetime.strftime(satscene.time_slot, options["dir"]) % values
        print path_template

        #options["filename"] = kwargs.get("filename")


