"""Shared plotting functions
"""
import json
import os

from collections import namedtuple, OrderedDict

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from osgeo import gdal
import numpy as np


def load_config():
    """Read config.json
    """
    config_path = os.path.join(os.path.dirname(__file__), '..', 'config.json')
    with open(config_path, 'r') as config_fh:
        config = json.load(config_fh)
    return config


def get_axes(extent=[102.2, 109.5, 8.5, 23.3], epsg=None):
    """Get transverse mercator axes (default to Vietnam extent)
    EPSG:4756
    """
    if epsg is not None:
        ax_proj = ccrs.epsg(epsg)
    else:
        x0, x1, y0, y1 = extent
        cx = x0 + ((x1 - x0) / 2)
        cy = y0 + ((y1 - y0) / 2)
        ax_proj = ccrs.LambertConformal(central_longitude=cx, central_latitude=cy)

    plt.figure(figsize=(6, 10), dpi=300)
    ax = plt.axes([0.025, 0.025, 0.95, 0.95], projection=ax_proj)
    proj = ccrs.PlateCarree()
    ax.set_extent(extent, crs=proj)
    set_ax_bg(ax)
    return ax


def save_fig(output_filename):
    plt.savefig(output_filename)


def set_ax_bg(ax, color='#c6e0ff'):
    """Set axis background color
    """
    ax.background_patch.set_facecolor(color)


def plot_basemap(ax, data_path, focus='VNM', neighbours=['VNM', 'CHN', 'LAO', 'KHM', 'THA'],country_border='white',plot_regions=True):
    """Plot countries and regions background
    """
    proj = ccrs.PlateCarree()

    states_filename = os.path.join(
        data_path,
        'Global_boundaries',
        'ne_10m_admin_0_countries_lakes.shp'
    )

    states_over_lakes_filename = os.path.join(
        data_path,
        'Global_boundaries',
        'ne_10m_admin_0_countries.shp'
    )

    provinces_filename = os.path.join(
        data_path,
        'Vietnam_boundaries',
        'who_boundaries',
        'who_provinces.shp'
    )

    lakes_filename = os.path.join(
        data_path,
        'Global_boundaries',
        'ne_10m_lakes.shp'
    )

    # Neighbours
    for record in shpreader.Reader(states_filename).records():
        country_code = record.attributes['ISO_A3']
        if country_code in neighbours:
            geom = record.geometry
            ax.add_geometries(
                [geom],
                crs=proj,
                edgecolor=country_border,
                facecolor='#e0e0e0',
                zorder=1)

    # Regions
    if plot_regions == True:
        for record in shpreader.Reader(provinces_filename).records():
            geom = record.geometry
            ax.add_geometries([geom], crs=proj, edgecolor='#ffffff', facecolor='#d2d2d2')

    # Lakes
    for record in shpreader.Reader(lakes_filename).records():
        name = record.attributes['name']
        geom = record.geometry
        ax.add_geometries(
            [geom],
            crs=proj,
            edgecolor='none',
            facecolor='#c6e0ff',
            zorder=1)

def plot_basemap_labels_large_region(ax, data_path):

        labels = [
            ('Vietnam', 108.633, 13.625, 9),
            ('Myanmar', 97.383, 21.535, 9),
            ('Malaysia', 99.404, 8.624, 9),
            ('Indonesia', 97.822, 3.338, 9),
            ('Singapore', 103.799, 1.472, 9),
            ('Cambodia', 105.25, 12.89, 9),
            ('Laos', 105.64, 16.55, 9),
            ('Thailand', 101.360, 16.257, 9),
            ('China', 108.08, 22.71, 9),
            ('South China Sea', 108.17, 17.37, 7)
        ]
        plot_basemap_labels(ax, data_path, labels, province_zoom=False, plot_regions=False)


def plot_basemap_labels(ax, data_path, labels=None, province_zoom=False, plot_regions=True):
    """Plot countries and regions background
    """
    proj = ccrs.PlateCarree()
    extent = ax.get_extent()

    if labels is None:
        labels = [
            ('Cambodia', 105.25, 12.89, 9),
            ('Laos', 105.64, 16.55, 9),
            ('Thailand', 103.64, 15.25, 9),
            ('China', 108.08, 22.71, 9),
            ('South China Sea', 108.17, 17.37, 7)]

        if plot_regions == True:
            labels = labels + [
            # Provinces
            ('An Giang', 105.182, 10.491, 5),
            ('Ba Ria-Vung Tau', 107.250, 10.510, 5),
            ('Bac Giang', 106.480, 21.357, 5),
            ('Bac Kan', 105.826, 22.261, 5),
            ('Bac Lieu', 105.489, 9.313, 5),
            ('Bac Ninh', 106.106, 21.109, 5),
            ('Ben Tre', 106.469, 10.118, 5),
            ('Binh Dinh', 108.951, 14.121, 5),
            ('Binh Duong', 106.658, 11.216, 5),
            ('Binh Phuoc', 106.907, 11.754, 5),
            ('Binh Thuan', 108.048, 11.117, 5),
            ('Ca Mau', 105.036, 9.046, 5),
            ('Can Tho', 105.530, 10.184, 5),
            ('Cao Bang', 106.087, 22.744, 5),
            # ('Da Nang', 109.634, 16.188, 5),
            ('Dak Lak', 108.212, 12.823, 5),
            ('Dak Nong', 107.688, 12.228, 5),
            ('Dien Bien', 103.022, 21.710, 5),
            ('Dong Nai', 107.185, 11.058, 5),
            ('Dong Thap', 105.608, 10.564, 5),
            ('Gia Lai', 108.241, 13.797, 5),
            ('Ha Giang', 104.979, 22.767, 5),
            ('Ha Nam', 105.966, 20.540, 5),
            ('Ha Noi', 105.700, 20.998, 5),
            ('Ha Tinh', 105.737, 18.290, 5),
            ('Hai Duong', 106.361, 20.930, 5),
            ('Hai Phong', 106.686, 20.798, 5),
            ('Hau Giang', 105.624, 9.784, 5),
            ('Ho Chi Minh', 106.697, 10.743, 5),
            ('Hoa Binh', 105.343, 20.684, 5),
            ('Hung Yen', 106.060, 20.814, 5),
            # ('Khanh Hoa', 111.307, 10.890, 5),
            ('Kien Giang', 104.942, 9.998, 5),
            ('Kon Tum', 107.875, 14.647, 5),
            ('Lai Chau', 103.187, 22.316, 5),
            ('Lam Dong', 108.095, 11.750, 5),
            ('Lang Son', 106.621, 21.838, 5),
            ('Lao Cai', 104.112, 22.365, 5),
            ('Long An', 106.171, 10.730, 5),
            ('Nam Dinh', 106.217, 20.268, 5),
            ('Nghe An', 104.944, 19.236, 5),
            ('Ninh Binh', 105.903, 20.170, 5),
            ('Ninh Thuan', 108.869, 11.705, 5),
            ('Phu Tho', 105.116, 21.308, 5),
            ('Phu Yen', 109.059, 13.171, 5),
            ('Quang Binh', 106.293, 17.532, 5),
            ('Quang Nam', 107.960, 15.589, 5),
            ('Quang Ngai', 108.650, 14.991, 5),
            ('Quang Ninh', 107.278, 21.245, 5),
            ('Quang Tri', 106.929, 16.745, 5),
            ('Soc Trang', 105.928, 9.558, 5),
            ('Son La', 104.070, 21.192, 5),
            ('Tay Ninh', 106.161, 11.404, 5),
            ('Thai Binh', 106.419, 20.450, 5),
            ('Thai Nguyen', 105.823, 21.692, 5),
            ('Thanh Hoa', 105.319, 20.045, 5),
            ('Thua Thien Hue', 107.512, 16.331, 5),
            ('Tien Giang', 106.309, 10.396, 5),
            ('Tra Vinh', 106.318, 9.794, 5),
            ('Tuyen Quang', 105.267, 22.113, 5),
            ('Vinh Long', 105.991, 10.121, 5),
            ('Vinh Phuc', 105.559, 21.371, 5),
            ('Yen Bai', 104.568, 21.776, 5),
        ]

    for text, x, y, size in labels:

        if province_zoom == True:
            size = 18

        if within_extent(x, y, extent):
            ax.text(
                x, y,
                text,
                alpha=0.7,
                size=size,
                horizontalalignment='center',
                transform=proj)



def within_extent(x, y, extent):
    xmin, xmax, ymin, ymax = extent
    return (xmin < x) and (x < xmax) and (ymin < y) and (y < ymax)


def scale_bar(ax, length=100, location=(0.5, 0.05), linewidth=3):
    """Draw a scale bar

    Adapted from https://stackoverflow.com/questions/32333870/how-can-i-show-a-km-ruler-on-a-cartopy-matplotlib-plot/35705477#35705477

    Parameters
    ----------
    ax : axes
    length : int
        length of the scalebar in km.
    location: tuple
        center of the scalebar in axis coordinates (ie. 0.5 is the middle of the plot)
    linewidth: float
        thickness of the scalebar.
    """
    # lat-lon limits
    llx0, llx1, lly0, lly1 = ax.get_extent(ccrs.PlateCarree())

    # Transverse mercator for length
    x = (llx1 + llx0) / 2
    y = lly0 + (lly1 - lly0) * location[1]
    tmc = ccrs.TransverseMercator(x, y)

    # Extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(tmc)

    # Scalebar location coordinates in metres
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]
    bar_xs = [sbx - length * 500, sbx + length * 500]

    # Plot the scalebar and label
    ax.plot(bar_xs, [sby, sby], transform=tmc, color='k', linewidth=linewidth)
    ax.text(sbx, sby + 10*length, str(length) + ' km', transform=tmc,
            horizontalalignment='center', verticalalignment='bottom', size=8)


def generate_weight_bins(weights, n_steps=9, width_step=0.01):
    """Given a list of weight values, generate <n_steps> bins with a width
    value to use for plotting e.g. weighted network flow maps.
    """
    min_weight = min(weights)
    max_weight = max(weights)

    width_by_range = OrderedDict()

    mins = np.linspace(min_weight, max_weight, n_steps)
    maxs = list(mins)
    maxs.append(max_weight*10)
    maxs = maxs[1:]

    assert len(maxs) == len(mins)

    for i, (min_, max_) in enumerate(zip(mins, maxs)):
        width_by_range[(min_, max_)] = (i+1) * width_step

    return width_by_range


Style = namedtuple('Style', ['color', 'zindex', 'label'])
Style.__doc__ += """: class to hold an element's styles

Used to generate legend entries, apply uniform style to groups of map elements
(See network_map.py for example.)
"""

def legend_from_style_spec(ax, styles, loc='lower left'):
    """Plot legend
    """
    handles = [
        mpatches.Patch(color=style.color, label=style.label)
        for style in styles.values()
    ]
    ax.legend(
        handles=handles,
        loc=loc
    )

def round_sf(x, places=1):
    """Round number to significant figures
    """
    if x == 0:
        return 0
    sign = x / abs(x)
    x = abs(x)
    exp = floor(log10(x)) + 1
    shift = 10 ** (exp - places)
    rounded = round(x / shift) * shift
    return rounded * sign

def get_data(filename):
    """Read in data (as array) and extent of each raster
    """
    gdal.UseExceptions()
    ds = gdal.Open(filename)
    data = ds.ReadAsArray()
    data[data < 0] = 0

    gt = ds.GetGeoTransform()

    # get the edge coordinates
    width = ds.RasterXSize
    height = ds.RasterYSize
    xres = gt[1]
    yres = gt[5]

    xmin = gt[0]
    xmax = gt[0] + (xres * width)
    ymin = gt[3] + (yres * height)
    ymax = gt[3]

    lat_lon_extent = (xmin, xmax, ymax, ymin)

    return data, lat_lon_extent