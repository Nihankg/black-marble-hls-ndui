"""Port of the supplied GEE NDUI workflow. Computation stays on Earth Engine."""
import argparse
from datetime import date, datetime, timezone
import json
import hashlib
import math
from pathlib import Path

import ee

ROOT = Path(__file__).resolve().parent
BAND = 'DNB_BRDF_Corrected_NTL'
VIIRS = 'NASA/VIIRS/002/VNP46A2'
L30 = 'NASA/HLS/HLSL30/v002'
S30 = 'NASA/HLS/HLSS30/v002'
HD_TRANSFORM = [30, 0, 0, 0, -30, 0]


def load_config(path):
    c = json.loads(Path(path).read_text())
    for period in ('pre', 'post', 'optical', 'reference'):
        start, end = map(date.fromisoformat, c[period])
        if start >= end:
            raise ValueError(f'{period}: start must precede exclusive end')
        if period in ('optical', 'reference') and end > date.fromisoformat(c['post'][0]):
            raise ValueError(f'{period} must end before the post-event window')
    if c['pre'][1] > c['post'][0]:
        raise ValueError('Before and after periods must not overlap')
    if c['min_valid_days'] < 1:
        raise ValueError('min_valid_days must be positive')
    if c['cloud_quality_min'] not in range(4) or c['cloud_class_max'] not in (0, 1):
        raise ValueError('Invalid cloud QA settings; cloudy classes are not supported')
    west, south, east, north = c['export_bbox']
    if not (-180 <= west < east <= 180 and -90 <= south < north <= 90):
        raise ValueError('export_bbox must be [west, south, east, north]')
    return c


def build_products(c):
    export_region = ee.Geometry.Rectangle(c['export_bbox'], 'EPSG:4326', False)
    if c.get('roi_asset'):
        filt = c['roi_filter']
        boundary = ee.FeatureCollection(c['roi_asset']).filter(
            ee.Filter.eq(filt['field'], filt['value']))
        if boundary.size().getInfo() == 0:
            raise ValueError('ROI asset filter returned no features')
        roi = boundary.geometry()
    else:
        # Portable example mode: different normalization region from Hatay asset.
        roi = export_region
    grid30 = ee.Projection('EPSG:32636', HD_TRANSFORM)
    source = ee.ImageCollection(VIIRS).filterBounds(roi)
    ref_source = source.filterDate(*c['reference'])
    if ref_source.size().getInfo() == 0:
        raise ValueError('No reference VIIRS scenes')
    native = ee.Image(ref_source.first()).select(BAND).projection()

    def actual_observation(img):
        q = img.select('QF_Cloud_Mask')
        good = (img.select('Mandatory_Quality_Flag').lte(1)
                .And(img.select('Snow_Flag').eq(0))
                .And(q.bitwiseAnd(1).eq(0))
                .And(q.rightShift(4).bitwiseAnd(3).gte(c['cloud_quality_min']))
                .And(q.rightShift(6).bitwiseAnd(3).lte(c['cloud_class_max']))
                .And(q.bitwiseAnd(1 << 9).eq(0))
                .And(q.bitwiseAnd(1 << 10).eq(0)))
        if c['mask_shadow']:
            good = good.And(q.bitwiseAnd(1 << 8).eq(0))
        return (img.select(BAND).rename('NTL').updateMask(good)
                .copyProperties(img, ['system:time_start'])
                .set('day', img.date().format('YYYY-MM-dd')))

    observations = source.map(actual_observation)
    empty = (ee.Image.constant(0).rename('NTL')
             .updateMask(ee.Image.constant(0)).setDefaultProjection(native))

    def summarize(period):
        subset = observations.filterDate(*period)
        dates = ee.List(subset.aggregate_array('day')).distinct().sort()

        def mosaic_day(day):
            day_date = ee.Date(day)
            return (subset.filterDate(day_date, day_date.advance(1, 'day')).mosaic()
                    .setDefaultProjection(native)
                    .set('system:time_start', day_date.millis()).set('day', day))

        daily = ee.ImageCollection.fromImages(dates.map(mosaic_day))
        safe = daily.merge(ee.ImageCollection.fromImages([empty]))
        count = (safe.count().rename('valid_days').unmask(0)
                 .setDefaultProjection(native).clip(roi))
        median = (safe.median().rename('NTL').setDefaultProjection(native)
                  .updateMask(count.gte(c['min_valid_days'])).clip(roi))
        return {'ntl': median, 'count': count, 'dates': dates}

    before, after, reference = [summarize(c[k]) for k in ('pre', 'post', 'reference')]
    p99 = reference['ntl'].reduceRegion(
        reducer=ee.Reducer.percentile([99]), geometry=roi,
        crs=native, scale=native.nominalScale(), maxPixels=1e13).get('NTL')

    def prepare_hls(img, nir_band):
        clear = img.select('Fmask').bitwiseAnd(30).eq(0)
        bands = img.select(['B3', 'B4', nir_band], ['green', 'red', 'nir'])
        bands = bands.updateMask(clear.And(bands.reduce(ee.Reducer.min()).gte(0)))
        green, red, nir = [bands.select(b) for b in ('green', 'red', 'nir')]
        ndvi = (nir.subtract(red).divide(nir.add(red))
                .updateMask(nir.add(red).gt(0)).rename('NDVI'))
        ndwi = (green.subtract(nir).divide(green.add(nir))
                .updateMask(green.add(nir).gt(0)).rename('NDWI'))
        return ndvi.addBands(ndwi).copyProperties(img, ['system:time_start'])

    landsat = ee.ImageCollection(L30).filterBounds(roi).filterDate(*c['optical'])
    sentinel = ee.ImageCollection(S30).filterBounds(roi).filterDate(*c['optical'])
    hls = landsat.map(lambda img: prepare_hls(img, 'B5')).merge(
        sentinel.map(lambda img: prepare_hls(img, 'B8A')))
    composite = hls.median().setDefaultProjection(grid30).clip(roi)
    ndvi, ndwi = composite.select('NDVI').max(0), composite.select('NDWI')
    common = (before['count'].gte(c['min_valid_days'])
              .And(after['count'].gte(c['min_valid_days'])).reproject(grid30)
              .And(ndwi.lte(0)).And(ndvi.mask()).clip(roi))

    def ndui(ntl):
        n = ntl.divide(ee.Number(p99)).clamp(0, 1).resample('bicubic').reproject(grid30)
        return (n.subtract(ndvi).divide(n.add(ndvi).max(0.001))
                .updateMask(common).rename('NDUI').reproject(grid30).clip(roi))

    products = {'Antakya_NDUI_before': ndui(before['ntl']),
                'Antakya_NDUI_after': ndui(after['ntl']),
                'Antakya_BlackMarble_coarse_before': before['ntl'],
                'Antakya_BlackMarble_coarse_after': after['ntl']}
    diagnostics = ee.Dictionary({
        'before_product_dates': before['dates'], 'after_product_dates': after['dates'],
        'reference_p99': p99, 'landsat_scene_count': landsat.size(),
        'sentinel_scene_count': sentinel.size(),
        'native_grid': native, 'hd_grid': grid30,
        'hd_pixel_size_m': grid30.nominalScale(),
        'before_max_valid_days': before['count'].reduceRegion(
            ee.Reducer.max(), export_region, crs=native,
            scale=native.nominalScale(), maxPixels=1e8).get('valid_days'),
        'after_max_valid_days': after['count'].reduceRegion(
            ee.Reducer.max(), export_region, crs=native,
            scale=native.nominalScale(), maxPixels=1e8).get('valid_days')})
    return products, diagnostics, export_region, native


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT / 'config.local.json')
    parser.add_argument('--authenticate', action='store_true')
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--export', action='store_true', help='Start four Google Drive tasks')
    action.add_argument('--download', action='store_true', help='Download four TIFFs locally')
    parser.add_argument('--output-dir', type=Path, default=ROOT.parent / 'outputs')
    parser.add_argument('--overwrite', action='store_true', help='Replace local downloaded TIFFs')
    args = parser.parse_args()
    c = load_config(args.config)
    if args.authenticate:
        ee.Authenticate()
    ee.Initialize(project=c['project'])
    products, diagnostic_graph, region, native = build_products(c)
    diagnostics = diagnostic_graph.getInfo()
    p99 = diagnostics.get('reference_p99')
    if p99 is None or p99 <= 0:
        raise ValueError('No positive normalization p99; inspect reference coverage')
    if not diagnostics['landsat_scene_count'] + diagnostics['sentinel_scene_count']:
        raise ValueError('No HLS scenes in the optical window')
    if not diagnostics.get('before_max_valid_days') or not diagnostics.get('after_max_valid_days'):
        raise ValueError('No valid VIIRS observations in at least one export period')
    print(json.dumps(diagnostics, indent=2))
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    manifest = {'created_utc': datetime.now(timezone.utc).isoformat(), 'config': c,
                'diagnostics': diagnostics, 'datasets': [VIIRS, L30, S30],
                'tasks': [], 'note': 'NDUI is a modelled index, not measured 30 m radiance.'}
    mode = 'download' if args.download else 'export' if args.export else 'check'
    manifest_path = out / f'{mode}_manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2))
    if args.download:
        import numpy as np
        import requests
        import rasterio
        from rasterio.warp import transform_bounds

        native_info = native.getInfo()
        manifest['files'] = []
        for name, image in products.items():
            hd = 'NDUI' in name
            crs = 'EPSG:32636' if hd else native_info['crs']
            base = HD_TRANSFORM if hd else native_info['transform']
            a, b, x0, d, e, y0 = base
            if b or d or a <= 0 or e >= 0:
                raise ValueError('Local download expects a north-up grid')
            west, south, east, north = transform_bounds(
                'EPSG:4326', crs, *c['export_bbox'], densify_pts=21)
            col0, col1 = math.floor((west-x0)/a), math.ceil((east-x0)/a)
            row0, row1 = math.floor((north-y0)/e), math.ceil((south-y0)/e)
            width, height = col1-col0, row1-row0
            # Explicit dimensions are essential: region is ignored with CRS+transform.
            transform = [a, 0, x0+col0*a, 0, e, y0+row0*e]
            if max(width, height) > 10000 or width*height*4 > 30*1024*1024:
                raise ValueError('Region too large for direct download; use --export')
            target = out / f'{name}.tif'
            if target.exists() and not args.overwrite:
                raise FileExistsError(f'{target} already exists; choose another output folder')
            print(f'Downloading {name} ({width} x {height})...', flush=True)
            url = image.clip(region).toFloat().unmask(
                value=-9999, sameFootprint=False).getDownloadURL({
                    'crs': crs, 'crs_transform': transform,
                    'dimensions': [width, height], 'format': 'GEO_TIFF'})
            partial = target.with_suffix('.tif.part')
            try:
                with requests.get(url, stream=True, timeout=(30, 300)) as response:
                    response.raise_for_status()
                    with partial.open('wb') as stream:
                        for chunk in response.iter_content(1024*1024):
                            stream.write(chunk)
                with rasterio.open(partial, 'r+') as raster:
                    raster.nodata = -9999
                    if raster.shape != (height, width) or raster.count != 1:
                        raise ValueError('Unexpected downloaded raster dimensions')
                    if raster.crs != rasterio.crs.CRS.from_string(crs):
                        raise ValueError('Unexpected downloaded raster CRS')
                    values = raster.read(1, masked=True)
                    valid = values.compressed()
                    if not len(valid) or not np.isfinite(valid).all():
                        raise ValueError('Downloaded raster has no valid finite values')
                    info = {'filename': target.name, 'shape': list(raster.shape),
                            'crs': str(raster.crs), 'transform': list(raster.transform)[:6],
                            'valid_pixels': int(values.count()),
                            'min': float(valid.min()), 'max': float(valid.max())}
                partial.replace(target)
            finally:
                partial.unlink(missing_ok=True)
            info['sha256'] = hashlib.sha256(target.read_bytes()).hexdigest()
            manifest['files'].append(info)
            manifest_path.write_text(json.dumps(manifest, indent=2))
            print(f'Saved {target}', flush=True)
    elif args.export:
        native_info = native.getInfo()
        for name, image in products.items():
            hd = 'NDUI' in name
            task = ee.batch.Export.image.toDrive(
                image=image.toFloat().unmask(value=-9999, sameFootprint=False),
                description=name, fileNamePrefix=name, folder=c['drive_folder'],
                region=region, crs='EPSG:32636' if hd else native_info['crs'],
                crsTransform=HD_TRANSFORM if hd else native_info['transform'],
                maxPixels=1e10, fileFormat='GeoTIFF',
                formatOptions={'cloudOptimized': True, 'noData': -9999})
            task.start()
            manifest['tasks'].append({'name': name, 'id': task.id})
            manifest_path.write_text(json.dumps(manifest, indent=2))
            print(f'Started {name}: {task.id}')
    else:
        print('Check completed. Use --download for local TIFFs or --export for Drive.')


if __name__ == '__main__':
    main()
