# (C) British Crown Copyright 2011 - 2012, Met Office
#
# This file is part of cartopy.
#
# cartopy is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cartopy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with cartopy.  If not, see <http://www.gnu.org/licenses/>.

import warnings

import matplotlib.pyplot as plt

import cartopy.crs as ccrs

from cartopy.tests.mpl import ImageTesting


@ImageTesting(['gridliner1'])
def test_gridliner():
    desired_gridline_prj = [ccrs.PlateCarree(), ccrs.OSGB()]
    projections = [ccrs.PlateCarree(), ccrs.OSGB(), ccrs.RotatedPole(37, 50)]
    ny, nx = 2, 4

    plt.figure(figsize=(10, 10))

    ax = plt.subplot(nx, ny, 1, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()
    ax.gridlines()

    ax = plt.subplot(nx, ny, 2, projection=ccrs.OSGB())
    ax.set_global()
    ax.coastlines()
    ax.gridlines()

    ax = plt.subplot(nx, ny, 3, projection=ccrs.OSGB())
    ax.set_global()
    ax.coastlines()
    ax.gridlines(ccrs.PlateCarree(), color='blue', linestyle='-')
    ax.gridlines(ccrs.OSGB())

    ax = plt.subplot(nx, ny, 4, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()
    ax.gridlines(ccrs.NorthPolarStereo(), alpha=0.5,
                 linewidth=1.5, linestyle='-')

    ax = plt.subplot(nx, ny, 5, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()
    osgb = ccrs.OSGB()
    ax.set_extent(tuple(osgb.x_limits) + tuple(osgb.y_limits), crs=osgb)
    ax.gridlines(osgb)

    ax = plt.subplot(nx, ny, 6, projection=ccrs.NorthPolarStereo())
    ax.set_global()
    ax.coastlines()
    ax.gridlines(alpha=0.5, linewidth=1.5, linestyle='-')

    ax = plt.subplot(nx, ny, 7, projection=ccrs.NorthPolarStereo())
    ax.set_global()
    ax.coastlines()
    osgb = ccrs.OSGB()
    ax.set_extent(tuple(osgb.x_limits) + tuple(osgb.y_limits), crs=osgb)
    ax.gridlines(osgb)

    ax = plt.subplot(nx, ny, 8,
                     projection=ccrs.Robinson(central_longitude=135))
    ax.set_global()
    ax.coastlines()
    ax.gridlines(ccrs.PlateCarree(), alpha=0.5, linewidth=1.5, linestyle='-')

    delta = 1.5e-2
    plt.subplots_adjust(left=0 + delta, right=1 - delta,
                        top=1 - delta, bottom=0 + delta)


@ImageTesting(['gridliner_labels'])
def test_grid_labels():
    plt.figure(figsize=(8, 10))

    crs_pc = ccrs.PlateCarree()
    crs_merc = ccrs.Mercator()
    crs_osgb = ccrs.OSGB()

    ax = plt.subplot(3, 2, 1, projection=crs_pc)
    ax.coastlines()
    ax.gridlines(add_labels=True)

    ax = plt.subplot(3, 2, 2, projection=crs_pc)
    ax.coastlines()
    with warnings.catch_warnings(record=True) as w:
        ax.gridlines(crs=crs_merc, add_labels=True)
        assert len(w) == 1
        message = w[0].message
        assert isinstance(message, UserWarning)
        assert str(message).find('labelling cancelled') >= 0

    ax = plt.subplot(3, 2, 3, projection=crs_merc)
    ax.coastlines()
    ax.gridlines(add_labels=True)

    ax = plt.subplot(3, 2, 4, projection=crs_osgb)
    ax.coastlines()
    with warnings.catch_warnings(record=True) as w:
        ax.gridlines(add_labels=True)
        assert len(w) == 1
        message = w[0].message
        assert isinstance(message, UserWarning)
        assert str(message).find('labelling cancelled') >= 0

    ax = plt.subplot(3, 2, 5, projection=crs_pc)
    ax.set_extent([-20, 10.0, 45.0, 70.0])
    ax.coastlines()
    ax.gridlines(add_labels=True)

    ax = plt.subplot(3, 2, 6, projection=crs_merc)
    ax.set_extent([-20, 10.0, 45.0, 70.0], crs=crs_pc)
    ax.coastlines()
    ax.gridlines(add_labels=True)

    # Increase margins between plots to stop them bumping into one another.
    plt.subplots_adjust(wspace=0.25, hspace=0.25)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
