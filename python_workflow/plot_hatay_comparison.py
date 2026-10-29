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

images = []

for phase in ["before", "after"]:
    print(f"Hatay {phase} indiriliyor...", flush=True)

    image = (
        products[f"Antakya_BlackMarble_coarse_{phase}"]
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

    path = OUTPUT / f"Hatay_BlackMarble_coarse_{phase}.tif"
    temporary = path.with_suffix(".tif.part")
    temporary.write_bytes(response.content)

    with rasterio.open(temporary, "r+") as src:
        src.nodata = -9999
        data = np.ma.masked_invalid(src.read(1, masked=True))

        if data.count() == 0:
            raise ValueError(f"{phase}: geçerli gözlem bulunamadı.")

        bounds = src.bounds

    temporary.replace(path)
    images.append(data)
    print(f"Kaydedildi: {path}")

# İki dönemde aynı renk ölçeğini kullan.
cmap = LinearSegmentedColormap.from_list(
    "black_marble",
    ["#000000", "#332200", "#ffaa00", "#ffff66"]
)
cmap.set_bad("white")

fig, axes = plt.subplots(
    1, 2, figsize=(11, 8),
    sharex=True, sharey=True,
    layout="constrained"
)

extent = [
    bounds.left, bounds.right,
    bounds.bottom, bounds.top
]

# Enlem-boylam gösterimindeki oranı yaklaşık düzelt.
aspect = 1 / np.cos(np.deg2rad((south + north) / 2))

for ax, data, period, label in zip(
    axes,
    images,
    [config["pre"], config["post"]],
    ["Before", "After"]
):
    plot = ax.imshow(
        data,
        extent=extent,
        cmap=cmap,
        vmin=0,
        vmax=50,
        interpolation="nearest",
        aspect=aspect,
    )

    ax.set_title(
        f"{label}\n{period[0]} to {period[1]} (end exclusive)"
    )
    ax.set_xlabel("Longitude")

axes[0].set_ylabel("Latitude")

fig.suptitle(
    "Hatay | Black Marble before and after the earthquake\n"
    "Native VIIRS grid (~500 m); white = NoData or outside Hatay",
    fontsize=13
)

fig.colorbar(
    plot,
    ax=axes,
    orientation="horizontal",
    shrink=0.8,
    pad=0.06,
    label="Radiance (nW cm$^{-2}$ sr$^{-1}$)"
)

png_path = OUTPUT / "hatay_coarse_comparison.png"
fig.savefig(png_path, dpi=300, facecolor="white")
print(f"PNG kaydedildi: {png_path}")

plt.show()