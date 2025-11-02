import argparse
import numpy as np
import pandas as pd
from typing import Optional


def _gaussian_kernel_2d(bw_lon: float, bw_lat: float, dlon: float, dlat: float, truncate: float = 3.0) -> np.ndarray:
    sx = max(bw_lon / max(dlon, 1e-12), 1e-6)
    sy = max(bw_lat / max(dlat, 1e-12), 1e-6)
    rx = int(np.ceil(truncate * sx))
    ry = int(np.ceil(truncate * sy))
    x = np.arange(-rx, rx + 1)
    y = np.arange(-ry, ry + 1)
    xx, yy = np.meshgrid(x, y, indexing="xy")
    k = np.exp(-0.5 * ((xx / sx) ** 2 + (yy / sy) ** 2))
    k /= k.sum() if k.sum() > 0 else 1.0
    return k.astype(np.float64, copy=False)


def _fft_convolve2d_same(img: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    h, w = img.shape
    kh, kw = kernel.shape
    pad_h = h + kh - 1
    pad_w = w + kw - 1
    fh = int(2 ** np.ceil(np.log2(pad_h)))
    fw = int(2 ** np.ceil(np.log2(pad_w)))
    F_img = np.fft.rfftn(img, s=(fh, fw))
    F_ker = np.fft.rfftn(kernel, s=(fh, fw))
    conv = np.fft.irfftn(F_img * F_ker, s=(fh, fw)).real
    sh = (kh - 1) // 2
    sw = (kw - 1) // 2
    out = conv[sh:sh + h, sw:sw + w]
    return out


def kde_grid(
    lons: np.ndarray,
    lats: np.ndarray,
    lon_range=None,
    lat_range=None,
    grid_res_deg: float = 0.01,
    bandwidth_lon_deg: float = 0.05,
    bandwidth_lat_deg: float = 0.05,
):
    lons = np.asarray(lons, dtype=float)
    lats = np.asarray(lats, dtype=float)
    m = np.isfinite(lons) & np.isfinite(lats)
    lons = lons[m]
    lats = lats[m]
    if lons.size == 0:
        raise ValueError("No valid coordinates provided")
    if lon_range is None:
        lon_min = float(np.nanmin(lons))
        lon_max = float(np.nanmax(lons))
    else:
        lon_min, lon_max = float(lon_range[0]), float(lon_range[1])
    if lat_range is None:
        lat_min = float(np.nanmin(lats))
        lat_max = float(np.nanmax(lats))
    else:
        lat_min, lat_max = float(lat_range[0]), float(lat_range[1])
    d = max(grid_res_deg, 1e-6)
    lon_edges = np.arange(lon_min, lon_max + d, d, dtype=float)
    lat_edges = np.arange(lat_min, lat_max + d, d, dtype=float)
    H, xedges, yedges = np.histogram2d(
        lons, lats, bins=[lon_edges, lat_edges]
    )
    dlon = float(np.diff(xedges).mean()) if xedges.size > 1 else d
    dlat = float(np.diff(yedges).mean()) if yedges.size > 1 else d
    K = _gaussian_kernel_2d(bandwidth_lon_deg, bandwidth_lat_deg, dlon, dlat)
    smoothed = _fft_convolve2d_same(H, K)
    n = float(lons.size)
    cell_area = dlon * dlat
    density = smoothed / max(n * cell_area, 1e-12)
    return density, xedges, yedges


def run_from_csv(
    in_csv: str,
    lon_col: str = "Longitude",
    lat_col: str = "Latitude",
    grid_res_deg: float = 0.01,
    bandwidth_deg: float = 0.05,
    lon_range=None,
    lat_range=None,
    out_npz: str = "kde_grid.npz",
    out_png: Optional[str] = None,
): 
    df = pd.read_csv(in_csv)
    # normalize headers: trim whitespace
    df.columns = [c.strip() for c in df.columns]
    if lon_col not in df.columns or lat_col not in df.columns:
        raise ValueError(f"Missing required columns: {lon_col}, {lat_col}")
    lons = df[lon_col].to_numpy(dtype=float, copy=False)
    lats = df[lat_col].to_numpy(dtype=float, copy=False)
    density, lon_edges, lat_edges = kde_grid(
        lons,
        lats,
        lon_range=lon_range,
        lat_range=lat_range,
        grid_res_deg=grid_res_deg,
        bandwidth_lon_deg=bandwidth_deg,
        bandwidth_lat_deg=bandwidth_deg,
    )
    np.savez_compressed(
        out_npz,
        density=density.astype(np.float32),
        lon_edges=lon_edges.astype(np.float32),
        lat_edges=lat_edges.astype(np.float32),
    )
    if out_png is not None:
        try:
            import matplotlib.pyplot as plt
            X, Y = np.meshgrid(lon_edges, lat_edges, indexing="xy")
            plt.figure(figsize=(6, 6))
            plt.pcolormesh(X, Y, density.T, shading="auto")
            plt.xlabel("Longitude")
            plt.ylabel("Latitude")
            plt.title("Spatial KDE")
            plt.colorbar(label="density (1/deg^2)")
            plt.tight_layout()
            plt.savefig(out_png, dpi=200)
            plt.close()
        except Exception:
            pass
    return density, lon_edges, lat_edges


def _parse_range(s: Optional[str]):
    if not s:
        return None
    parts = [p.strip() for p in s.split(",")]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("Range must be 'min,max'")
    return float(parts[0]), float(parts[1])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", required=True)
    p.add_argument("--lon-col", default="Longitude")
    p.add_argument("--lat-col", default="Latitude")
    p.add_argument("--grid-res", type=float, default=0.01)
    p.add_argument("--bandwidth", type=float, default=0.05)
    p.add_argument("--lon-range", type=str, default=None)
    p.add_argument("--lat-range", type=str, default=None)
    p.add_argument("--out-npz", default="kde_grid.npz")
    p.add_argument("--out-png", default=None)
    args = p.parse_args()
    lon_range = _parse_range(args.lon_range)
    lat_range = _parse_range(args.lat_range)
    run_from_csv(
        args.csv,
        lon_col=args.lon_col,
        lat_col=args.lat_col,
        grid_res_deg=args.grid_res,
        bandwidth_deg=args.bandwidth,
        lon_range=lon_range,
        lat_range=lat_range,
        out_npz=args.out_npz,
        out_png=args.out_png,
    )


if __name__ == "__main__":
    main()

