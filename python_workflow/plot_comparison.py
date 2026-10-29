"""Align four exported rasters and save a common-extent comparison."""
from pathlib import Path
import json
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    '--data-dir',
    type=Path,
    default=Path(__file__).resolve().parent.parent / 'outputs'
)
parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent.parent / 'outputs')
args = parser.parse_args()
DATA = args.data_dir
OUTPUT = args.output_dir
OUTPUT.mkdir(parents=True, exist_ok=True)
with rasterio.open(DATA / 'Antakya_NDUI_before.tif') as src:
    grid = dict(crs=src.crs, transform=src.transform,
                width=src.width, height=src.height)
    bounds = src.bounds
extent = (bounds.left / 1000, bounds.right / 1000,
          bounds.bottom / 1000, bounds.top / 1000)

def aligned(filename):
    destination = np.full((grid['height'], grid['width']), np.nan, dtype='float32')
    with rasterio.open(DATA / filename) as src:
        reproject(
            source=rasterio.band(src, 1), destination=destination,
            src_transform=src.transform, src_crs=src.crs, src_nodata=src.nodata,
            dst_transform=grid['transform'], dst_crs=grid['crs'],
            dst_nodata=np.nan, resampling=Resampling.nearest
        )
    return np.ma.masked_invalid(destination)

# Python exports use Antakya names; earlier GEE coarse exports used Hatay.
coarse_prefix = 'Antakya' if (DATA / 'Antakya_BlackMarble_coarse_before.tif').exists() else 'Hatay'
files = [
    f'{coarse_prefix}_BlackMarble_coarse_before.tif', 'Antakya_NDUI_before.tif',
    f'{coarse_prefix}_BlackMarble_coarse_after.tif', 'Antakya_NDUI_after.tif'
]
images = [aligned(name) for name in files]
cmap = LinearSegmentedColormap.from_list(
    'black_marble', ['#000000', '#332200', '#ffaa00', '#ffff66'])
cmap.set_bad('#d6d9df')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
fig, axes = plt.subplots(2, 2, figsize=(10, 11), sharex=True, sharey=True)
fig.subplots_adjust(left=.12, right=.97, top=.88, bottom=.18, wspace=.06, hspace=.10)
plots = []
for i, (ax, image) in enumerate(zip(axes.flat, images)):
    coarse = i % 2 == 0
    plots.append(ax.imshow(image, extent=extent, origin='upper', cmap=cmap,
                           vmin=0 if coarse else -1, vmax=50 if coarse else 1,
                           interpolation='nearest', aspect='equal'))
    ax.tick_params(labelsize=8)
    if i >= 2:
        ax.set_xlabel('UTM easting (km)')
    if coarse:
        ax.set_ylabel('UTM northing (km)')
axes[0, 0].set_title('Black Marble | native ~500 m', fontsize=12, pad=12)
axes[0, 1].set_title('HLS-assisted NDUI | 30 m grid', fontsize=12, pad=12)
for row, name in enumerate(['Before earthquake', 'After earthquake']):
    axes[row, 0].annotate(name, xy=(-.24, .5), xycoords='axes fraction',
                           rotation=90, ha='center', va='center', weight='bold')
for col, label in enumerate(['Black Marble radiance (nW cm$^{-2}$ sr$^{-1}$)',
                              'NDUI (dimensionless)']):
    pos = axes[1, col].get_position()
    cax = fig.add_axes([pos.x0, .105, pos.width, .014])
    fig.colorbar(plots[col], cax=cax, orientation='horizontal', label=label)
fig.suptitle('Antakya | Nighttime-light visualization comparison',
             fontsize=16, weight='bold', y=.965)
fig.text(.5, .925, 'Identical map extent and alignment across all four panels',
         ha='center', color='#444444')
fig.text(.5, .010,
         'Gray: NoData. Coarse rasters use nearest-neighbor reprojection; no visual smoothing.\n'
         '30 m NDUI is a modelled index, not a measured 30 m radiance product.\n'
         'Separate column scales: colors are not numerically equivalent. Sources: NASA VIIRS / HLS.',
         ha='center', va='bottom', fontsize=7, color='#444444', linespacing=1.3)
for suffix in ['png', 'pdf']:
    fig.savefig(OUTPUT / f'antakya_comparison.{suffix}', dpi=300, facecolor='white')
plt.close(fig)
(OUTPUT / 'alignment.json').write_text(json.dumps({
    'target_crs': str(grid['crs']), 'target_transform': list(grid['transform'])[:6],
    'shape': [grid['height'], grid['width']], 'extent': list(bounds),
    'resampling': 'nearest', 'files': files,
    'valid_pixels': {name: int(image.count()) for name, image in zip(files, images)},
    'dates': 'Use dates in the generating GEE script; GeoTIFF names do not encode them.'
}, indent=2))
print('Saved PNG, PDF and alignment metadata to', OUTPUT)
