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

import numpy as np
from matplotlib.testing.decorators import image_comparison as mpl_image_comparison
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.path import Path
import shapely.geometry

import cartopy.crs as ccrs
import cartopy.mpl_integration.patch as cpatch

from cartopy.tests.mpl import image_comparison


@image_comparison(baseline_images=['poly_interiors'])
def test_polygon_interiors():
    
    ax = plt.subplot(211, projection=ccrs.PlateCarree())
    ax.coastlines()    
    ax.set_global() # XXX could be the default???
   
    pth = Path([[0, 45], [60, 45], [60, -45], [0, -45], [0, -45], 
                [10, 20], [10, -20], [40, -20], [40, 20], [10, -20]], 
               [1, 2, 2, 2, 79, 1, 2, 2 , 2, 79])
    
    patches_native = []
    patches = []
    for geos in cpatch.path_to_geos(pth):
        for pth in cpatch.geos_to_path(geos):
            patches.append(mpatches.PathPatch(pth))
        
        # buffer by 10 degrees (leaves a small hole in the middle)
        geos_buffered = geos.buffer(10)
        for pth in cpatch.geos_to_path(geos_buffered):
            patches_native.append(mpatches.PathPatch(pth))
    
    collection = PatchCollection(patches_native, facecolor='red', alpha=0.4,
                                 transform=ax.projection
                                 )
    ax.add_collection(collection)

    collection = PatchCollection(patches, facecolor='yellow', alpha=0.4,
                                 transform=ccrs.Geodetic()
                                 )
    
    ax.add_collection(collection)
    
    
    # test multiple interior polygons
    ax = plt.subplot(212, projection=ccrs.PlateCarree(), xlim=[-5, 15], ylim=[-5, 15])
    ax.coastlines()
    
    exterior = np.array(shapely.geometry.box(0, 0, 12, 12).exterior.coords)
    interiors = [
                 np.array(shapely.geometry.box(1, 1, 2, 2, ccw=False).exterior.coords),
                 np.array(shapely.geometry.box(1, 8, 2, 9, ccw=False).exterior.coords),
                 ]

    poly = shapely.geometry.Polygon(exterior, interiors)
    
    patches = []
    for pth in cpatch.geos_to_path(poly):
        patches.append(mpatches.PathPatch(pth))
    
    collection = PatchCollection(patches, facecolor='yellow', alpha=0.4,
                                 transform=ccrs.Geodetic()
                                 )
    ax.add_collection(collection)
    

@image_comparison(baseline_images=['contour_with_interiors'])
def test_contour_interiors():
#    ############## produces a polygon with multiple holes:
    nx, ny = 10, 10
    numlev = 2
    lons, lats = np.meshgrid(np.linspace(-50, 50, nx), np.linspace(-45, 45, ny))
    data = np.sin(np.sqrt(lons**2 + lats**2))
    
    ax = plt.subplot(221, projection=ccrs.PlateCarree())
    ax.set_global()
    plt.title("Native projection")
    plt.contourf(lons, lats, data, numlev, transform=ccrs.PlateCarree())
    ax.coastlines()
    
    plt.subplot(222, projection=ccrs.PlateCarree())
    plt.title("Non-native projection")
    ax = plt.gca()
    ax.set_global()
    plt.contourf(lons, lats, data, numlev, transform=ccrs.Geodetic())
    ax.coastlines()    
    
    ############## produces singular polygons (zero area polygons)
    
    numlev = 2
    x, y = np.meshgrid(np.arange(-5.5, 5.5, 0.25), np.arange(-5.5, 5.5, 0.25))
    dim = x.shape[0]
    data = Z = np.sin(np.sqrt(x**2 + y**2))
    lats = np.arange(dim) + 30
    lons = np.arange(dim) - 20
        
    ax = plt.subplot(223, projection=ccrs.PlateCarree())
    ax.set_global()
    plt.title("Native projection")
    plt.contourf(lons, lats, data, numlev, transform=ccrs.PlateCarree())
    ax.coastlines()    
    
    plt.subplot(224, projection=ccrs.PlateCarree())
    plt.title("Non-native projection")
    ax = plt.gca()
    ax.set_global()
    cs = plt.contourf(lons, lats, data, numlev, transform=ccrs.Geodetic())
    
    ax.coastlines()    
    

if __name__=='__main__':
    import nose
    nose.runmodule(argv=['-s','--with-doctest'], exit=False)