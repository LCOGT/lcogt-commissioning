from pathlib import Path
import subprocess

import numpy as np
from astropy.io import ascii


def _assert_table_matches_reference(generated_path, reference_path, rtol=1e-5, atol=1e-5):
    generated = ascii.read(generated_path)
    reference = ascii.read(reference_path)

    assert generated.colnames == reference.colnames, (
        f"Column mismatch for {generated_path.name}: {generated.colnames} != {reference.colnames}"
    )

    for column in generated.colnames:
        generated_values = np.asarray(generated[column])
        reference_values = np.asarray(reference[column])

        if np.issubdtype(generated_values.dtype, np.number) or np.issubdtype(reference_values.dtype, np.number):
            assert np.allclose(generated_values, reference_values, rtol=rtol, atol=atol, equal_nan=True), (
                f"Column '{column}' in {generated_path.name} differs from reference beyond tolerance"
            )
        else:
            assert np.array_equal(generated_values, reference_values), (
                f"Column '{column}' in {generated_path.name} differs from reference"
            )


def test_noisegainmef_runs_on_ep60_data(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "testdata" / "ep60"
    files = sorted(data_dir.glob("*.fits.fz"))

    assert files, f"No FITS files found in {data_dir}"
    assert len(files) >= 4, "Expected at least a bias/flat set for a meaningful noisegain run"

    cmd = [
        "noisegainmef",
        "--quadrant",
        "--makepng",
        "--readmode",
        "full_frame",
        "--sortby",
        "filterlevel",
        *[str(path) for path in files],
    ]

    result = subprocess.run(
        cmd,
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        "noisegainmef failed on testdata/ep60\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    expected_outputs = {
        "ptc_full_frame_dateobs_flux.png",
        "ptc_full_frame_level_flux.png",
        "ptc_full_frame_levelgain.png",
        "ptc_full_frame_ptc.png",
        "ptc_full_frame_texplevel.png",
        "ptc_data_0_ll.dat",
        "ptc_data_0_lr.dat",
        "ptc_data_0_ul.dat",
        "ptc_data_0_ur.dat",
    }

    produced = {path.name for path in tmp_path.iterdir() if path.is_file()}
    missing = sorted(expected_outputs - produced)
    assert not missing, f"Missing expected noisegain outputs in {tmp_path}: {missing}"

    for filename in sorted({name for name in expected_outputs if name.endswith(".dat") }):
        generated_path = tmp_path / filename
        reference_path = data_dir / filename
        assert reference_path.exists(), f"Reference file missing: {reference_path}"
        _assert_table_matches_reference(generated_path, reference_path)



