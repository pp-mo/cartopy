"""
Function to find the 2d rectangular region in one Iris coordinate system that
covers the whole of a 2d region in a different coordinate system.

Uses numerical search for projected coordinate minima and maxima.

EXPERIMENTAL TEST CODE : not public, needs publishing and some tests.

CURRENT STATUS:
notional, untested.

"""
import numpy as np


_DO_DEBUG = True
#_DO_DEBUG = False


def transform_xy_arrays(crs_from, x, y, crs_to):
    """
    Transform 2d points between cartopy coordinate reference systems.

    NOTE: copied private function from iris.analysis.cartography.

    Args:

    * crs_from, crs_to (:class:`cartopy.crs.Projection`):
        The coordinate reference systems.
    * x, y (arrays):
        point locations defined in 'crs_from'.

    Returns:
        x, y :  Arrays of locations defined in 'crs_to'.

    """
    x = np.asarray(x)
    y = np.asarray(y)
    pts = crs_to.transform_points(crs_from, x, y)
    return pts[..., 0], pts[..., 1]


def meshgrid_linspace_2d(x_min, x_max, y_min, y_max, nx, ny):
    x_pts = np.linspace(x_min, x_max, nx)
    y_pts = np.linspace(y_min, y_max, ny)
    xs, ys = np.meshgrid(x_pts, y_pts)
    return xs, ys


def argmin_2d(array):
    argmin_1d = np.argmin(array)
    j_y, i_x = np.unravel_index(argmin_1d, array.shape)
    return i_x, j_y


def xy_for_min_value_in_xy_region(x0, x1, y0, y1,
                                  values_function_of_xs_ys,
                                  nx_samples=5, ny_samples=5,
                                  depth=0):
    done = False
    while not done:
        if _DO_DEBUG:
            print 'find :  x in ({}, {}), y in ({}, {})'.format(x0, x1, y0, y1)
        # Sample at a fixed number of points to look for an improved minimum.
        x_samples, y_samples = meshgrid_linspace_2d(x0, x1, y0, y1,
                                                    nx_samples, ny_samples)
        samples = values_function_of_xs_ys(x_samples, y_samples)

        # Find location of minimum in the sampled points.
        i_min, j_min = argmin_2d(samples)
        i_max, j_max = argmin_2d(-samples)

        # See if we have achieved any improvement to the minimum value.
        if _DO_DEBUG:
            print 'Difference : ', samples[j_max, i_max] - samples[j_min, i_min]
        # Recurse into a region around the location of the current best.
        # N.B. this is tail recursion, so could be looped ...
        j0, j1 = max(0, j_min - 1), min(ny_samples - 1, j_min + 1)
        i0, i1 = max(0, i_min - 1), min(nx_samples - 1, i_min + 1)
        next_x0, next_x1 = x_samples[0, i0], x_samples[0, i1]
        next_y0, next_y1 = y_samples[j0, 0], y_samples[j1, 0]
        done = (((next_x0 < x0 or next_x1 > x1 or
                  next_y0 < y0 or next_y1 > y1) or
                 (next_x0 <= x0) and next_x1 >= x1 and
                 next_y0 <= y0 and next_y1 >= y1))
        x0, x1, y0, y1 = next_x0, next_x1, next_y0, next_y1
        depth = depth + 1

    result = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    if _DO_DEBUG:
        print 'depth = ', depth

    return result


def projected_extent(source_coord_system,
                     source_min_x, source_max_x, source_min_y, source_max_y,
                     target_coord_system,
                     n_source_x_testpoints=77, n_source_y_testpoints=77):
    """
    Determine the region in the target crs that encloses all of the points from
    the target crs within the given 2d rectangle.

    Optional args:
        n_source_x_testpoints, n_source_y_testpoints (int):
            Control the number of sample points used in the initial search for
            local extrema.

    Returns:
        (target_min_x, target_max_x, target_min_y, target_max_y)

    """
    # Tweak to allow passing Iris coord systems :  Probably need removing.
    src_crs, tgt_crs = source_coord_system, target_coord_system
    if hasattr(src_crs, 'as_cartopy_crs'):
        src_crs = source_coord_system.as_cartopy_crs()
    if hasattr(tgt_crs, 'as_cartopy_crs'):
        tgt_crs = target_coord_system.as_cartopy_crs()

    # Subsidiary functions whose minima are min/max of x/y target coords.
    def target_x_coord_positive_function(src_x_points, src_y_points):
        tgt_x_points, tgt_y_points = transform_xy_arrays(
            src_crs, src_x_points, src_y_points, tgt_crs)
        return tgt_x_points

    def target_x_coord_negative_function(src_x_points, src_y_points):
        tgt_x_points, tgt_y_points = transform_xy_arrays(
            src_crs, src_x_points, src_y_points, tgt_crs)
        return -tgt_x_points

    def target_y_coord_positive_function(src_x_points, src_y_points):
        tgt_x_points, tgt_y_points = transform_xy_arrays(
            src_crs, src_x_points, src_y_points, tgt_crs)
        return tgt_y_points

    def target_y_coord_negative_function(src_x_points, src_y_points):
        tgt_x_points, tgt_y_points = transform_xy_arrays(
            src_crs, src_x_points, src_y_points, tgt_crs)
        return -tgt_y_points

    # Function to get the xy source coords of a target coord extremum.
    def xy_for_min(coordinate_value_function):
        # Get the source xy for the required extreme of a target coordinate.
        src_x, src_y = xy_for_min_value_in_xy_region(
            source_min_x, source_max_x, source_min_y, source_max_y,
            values_function_of_xs_ys=coordinate_value_function,
            nx_samples=n_source_x_testpoints, ny_samples=n_source_y_testpoints)
        return src_x, src_y

    # Calculate minimum target-cs x coordinate value.
    src_x, src_y = xy_for_min(target_x_coord_positive_function)
    tgt_x_min = transform_xy_arrays(
        src_crs, [src_x], [src_y], tgt_crs)[0][0]

    # Calculate maximum target-cs x coordinate value.
    src_x, src_y = xy_for_min(target_x_coord_negative_function)
    tgt_x_max = transform_xy_arrays(src_crs, [src_x], [src_y], tgt_crs)[0][0]

    # Calculate minimum target-cs y coordinate value.
    src_x, src_y = xy_for_min(target_y_coord_positive_function)
    tgt_y_min = transform_xy_arrays(src_crs, [src_x], [src_y], tgt_crs)[1][0]

    # Calculate maximum target-cs y coordinate value.
    src_x, src_y = xy_for_min(target_y_coord_negative_function)
    tgt_y_max = transform_xy_arrays(src_crs, [src_x], [src_y], tgt_crs)[1][0]

    return tgt_x_min, tgt_x_max, tgt_y_min, tgt_y_max


def do_test():
    import cartopy.crs as ccrs
    crs_latlon = ccrs.PlateCarree()

    # Define a source region defined in rotated-pole coords.
    pole_lat = 53.4
    pole_lon = -173.2
    crs_rotated = ccrs.RotatedPole(
        pole_longitude=pole_lon, pole_latitude=pole_lat)
    x0, x1 = -60.0, 60.0
    y0, y1 = -45.0, 20.0
    xx0, xx1, yy0, yy1 = projected_extent(
        source_coord_system=crs_rotated,
        source_min_x=x0, source_max_x=x1,
        source_min_y=y0, source_max_y=y1,
        target_coord_system=crs_latlon)
    print xx0, xx1, yy0, yy1


if __name__ == '__main__':
    do_test()
