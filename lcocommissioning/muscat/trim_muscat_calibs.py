#!/usr/bin/env python3

import argparse
import os
import re

from astropy.io import fits


# Hardcoded output date settings
OUTPUT_DATE_TAG = "20260602"  # YYYYMMDD for filename
OUTPUT_DATE_OBS = "2026-06-02T15:00:00.000"  # DATE-OBS header value

TRIM_LEFT= []


def find_image_hdu_index(hdul):
    for i, hdu in enumerate(hdul):
        if getattr(hdu, "data", None) is not None:
            if hasattr(hdu.data, "ndim") and hdu.data.ndim >= 2:
                return i
    raise RuntimeError("No image HDU with 2D data found.")


def get_instrume(hdul, image_hdu_index):
    for hdr in (hdul[0].header, hdul[image_hdu_index].header):
        val = hdr.get("INSTRUME")
        if val:
            return str(val).strip()
    raise RuntimeError("INSTRUME keyword not found.")


def trim_image(data, side):
    if data.shape[-1] != 2080 or data.shape[-2] != 2048:
        raise RuntimeError(f"Expected image shape (2048, 2080), got {data.shape[-2:]}.")
    if side == "left":
        return data[..., :, 32:]
    if side == "right":
        return data[..., :, :2048]
    raise RuntimeError(f"Unknown trim side: {side}")


def build_output_name(input_path):
    base = os.path.basename(input_path)
    # Replace first YYYYMMDD token if present; otherwise prepend it.
    if re.search(r"\d{8}", base):
        out = re.sub(r"\d{8}", OUTPUT_DATE_TAG, base, count=1)
    else:
        out = f"{OUTPUT_DATE_TAG}_{base}"
    if not out.endswith(".fits.fz"):
        out = re.sub(r"(\.fz)?$", ".fits.fz", out)
    return out


def main():
    parser = argparse.ArgumentParser(description="Trim MUSCAT compressed FITS calibrations.")
    parser.add_argument("input_fits_fz", help="Input compressed FITS file (.fits.fz)")
    parser.add_argument("-o", "--output", help="Output filename (.fits.fz). If omitted, auto-generated.")
    args = parser.parse_args()

    with fits.open(args.input_fits_fz) as hdul:
        img_idx = find_image_hdu_index(hdul)
        instrume = get_instrume(hdul, img_idx)

        side = "left" if instrume in TRIM_LEFT else "right"
        print ("instrume:", instrume, "-> trimming", side, "side")
        trimmed = trim_image(hdul[img_idx].data, side)

        primary_header = hdul[0].header.copy()
        image_header = hdul[img_idx].header.copy()

    #primary_header["DATE-OBS"] = OUTPUT_DATE_OBS
    image_header["DATE-OBS"] = OUTPUT_DATE_OBS

    output_path = args.output if args.output else os.path.join(
        os.path.dirname(args.input_fits_fz), build_output_name(args.input_fits_fz)
    )
    print ("file names: \n input:", args.input_fits_fz, "\noutput:", output_path)

    phdu = fits.PrimaryHDU(header=primary_header)
    chdu = fits.CompImageHDU(data=trimmed, header=image_header, name="BPM")
    fits.HDUList([phdu, chdu]).writeto(output_path, overwrite=True)


if __name__ == "__main__":
    main()
