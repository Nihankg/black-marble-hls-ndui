from pathlib import Path
import math

import ee
import numpy as np
import requests
import rasterio
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from earth_engine_pipeline import load_config, build_products

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT.parent / "outputs"
OUTPUT.mkdir(exist_ok=True)

# Mevcut projenin tarih ve kalite ayarları.
config = load_config(ROOT / "config.local.json")
if not config.get("roi_asset"):
    raise ValueError("Set roi_asset and roi_filter in config.local.json to an accessible Hatay boundary asset.")
ee.Initialize(project=config["project"])

hatay = (
    ee.FeatureCollection(config["roi_asset"])
    .filter(ee.Filter.eq(
        config["roi_filter"]["field"],
        config["roi_filter"]["value"]
    ))
    .geometry()
)

products, _, _, native = build_products(config)
projection = native.getInfo()

# Hatay sınırlarını kapsayan, orijinal VIIRS ağına oturan alan.
coordinates = hatay.bounds().getInfo()["coordinates"][0]
west = min(point[0] for point in coordinates)
east = max(point[0] for point in coordinates)
south = min(point[1] for point in coordinates)
north = max(point[1] for point in coordinates)

a, b, x0, d, e, y0 = projection["transform"]

if projection["crs"] != "EPSG:4326" or b != 0 or d != 0:
    raise ValueError("Beklenmeyen VIIRS projeksiyonu.")

col0 = math.floor((west - x0) / a)
col1 = math.ceil((east - x0) / a)
row0 = math.floor((north - y0) / e)
row1 = math.ceil((south - y0) / e)

width = col1 - col0
height = row1 - row0

transform = [
    a, 0, x0 + col0 * a,
    0, e, y0 + row0 * e,
]

# Küçük doğrudan indirme isteği sınırını kontrol et.
if max(width, height) > 10000 or width * height * 4 > 30 * 1024**2:
    raise ValueError("Alan doğrudan indirmek için fazla büyük.")

print("Deprem öncesi Hatay indiriliyor...", flush=True)

image = (
    products["Antakya_BlackMarble_coarse_before"]
    .clip(hatay)
    .toFloat()
    .unmask(value=-9999, sameFootprint=False)
)

url = image.getDownloadURL({
    "crs": projection["crs"],
    "crs_transform": transform,
    "dimensions": [width, height],
    "format": "GEO_TIFF",
})

response = requests.get(url, timeout=(30, 300))
response.raise_for_status()

path = OUTPUT / "Hatay_BlackMarble_coarse_before.tif"
temporary = path.with_suffix(".tif.part")
temporary.write_bytes(response.content)

with rasterio.open(temporary, "r+") as src:
    src.nodata = -9999
    data = np.ma.masked_invalid(src.read(1, masked=True))
    bounds = src.bounds

    if data.count() == 0:
        raise ValueError("Geçerli gözlem bulunamadı.")

temporary.replace(path)
print(f"TIFF kaydedildi: {path}")

cmap = LinearSegmentedColormap.from_list(
    "black_marble",
    ["#000000", "#332200", "#ffaa00", "#ffff66"]
)
cmap.set_bad("white")

fig, ax = plt.subplots(figsize=(7, 9), layout="constrained")

plot = ax.imshow(
    data,
    extent=[
        bounds.left, bounds.right,
        bounds.bottom, bounds.top
    ],
    cmap=cmap,
    vmin=0,
    vmax=50,
    interpolation="nearest",
    aspect=1 / np.cos(np.deg2rad((south + north) / 2)),
)

ax.set_title(
    "Hatay | Black Marble — before earthquake\n"
    f"{config['pre'][0]} to {config['pre'][1]} (end exclusive)"
)
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")

fig.colorbar(
    plot, ax=ax,
    orientation="horizontal",
    pad=0.07,
    label="Radiance (nW cm$^{-2}$ sr$^{-1}$)"
)

fig.supxlabel(
    "Native VIIRS grid (~500 m). White: NoData or outside Hatay.",
    fontsize=9
)

png_path = OUTPUT / "hatay_coarse_before.png"
fig.savefig(png_path, dpi=300, facecolor="white")
print(f"PNG kaydedildi: {png_path}")

plt.show()