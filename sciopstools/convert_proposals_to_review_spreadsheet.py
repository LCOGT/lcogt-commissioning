#!/usr/bin/env python3
"""Convert Proposal TSV data into a review spreadsheet TSV.

The input file uses tab-delimited ScienceApplication-style proposal columns.
The output groups requested time by semester and telescope size, with review
adjustment and comment columns for each semester/size group.

Example:
    python3 convert_proposals_to_review_spreadsheet.py input.tsv output.tsv
"""

import csv
from io import StringIO
import re
import sys
from collections import Counter


COLUMN_RE = re.compile(
    r'^(?P<semester>\d{4}[AB])\s+'
    r'(?P<aperture>0M4|1M0|2M0)-'
    r'(?P<family>[A-Z0-9]+)-'
    r'(?P<camera>[A-Z0-9]+)\s+'
    r'(?P<mode>Queue|RR|TC)$'
)


def canonical_instrument(aperture, family, camera):
    family = family.upper()
    camera = camera.upper()

    if aperture == "0M4":
        if camera == "QHY600":
            return "QHY"
        return None  # SBIG ignored

    if aperture == "1M0":
        if family == "NRES":
            return "NRES"
        if camera == "SINISTRO":
            return "Sinistro"
        return None

    if aperture == "2M0":
        if family == "FLOYDS":
            return "FLOYDS"
        if camera == "MUSCAT":
            return "MuSCAT"
        return None  # SPECTRAL ignored

    return None


def checkbox_value(flag):
    return "TRUE" if flag else ""


def parse_tags(raw_tags):
    return {t.strip().lower() for t in (raw_tags or "").split("|") if t.strip()}


def safe_float(value):
    if value is None:
        return 0.0
    s = str(value).strip()
    if not s:
        return 0.0
    try:
        return float(s)
    except ValueError:
        return 0.0


def number_or_blank(value):
    if value == 0:
        return ""
    return int(value) if float(value).is_integer() else value


def blank_row(row):
    return not any((value or "").strip() for value in row.values())


def validate_unique_headers(headers):
    duplicates = [header for header, count in Counter(headers).items() if count > 1]
    if duplicates:
        print("Duplicate input headers found:", file=sys.stderr)
        for header in duplicates:
            print(f"  {header}", file=sys.stderr)
        sys.exit(1)


def semester_sort_key(semester):
    return (int(semester[:4]), 0 if semester[4] == "A" else 1)


def telescope_group(instrument):
    if instrument in {"FLOYDS", "MuSCAT"}:
        return "2m"
    if instrument in {"NRES", "Sinistro"}:
        return "1m"
    if instrument == "QHY":
        return "0m4"
    return None


def convertproposals(input_tsv):
    
    reader = csv.DictReader(StringIO(input_tsv), delimiter="\t")
    input_headers = reader.fieldnames or []
    validate_unique_headers(input_headers)
    input_rows = [row for row in reader if not blank_row(row)]

    semester_columns = []
    seen_semester_cols = set()

    for col in input_headers:
        m = COLUMN_RE.match(col)
        if not m:
            continue

        semester = m.group("semester")
        aperture = m.group("aperture")
        family = m.group("family")
        camera = m.group("camera")
        mode = m.group("mode")

        instrument = canonical_instrument(aperture, family, camera)
        if instrument is None:
            continue

        key = (semester, instrument, mode)
        if key not in seen_semester_cols:
            seen_semester_cols.add(key)
            semester_columns.append(key)

    instrument_order = {"FLOYDS": 0, "MuSCAT": 1, "Sinistro": 2, "QHY": 3}
    mode_order = {"Queue": 0, "RR": 1, "TC": 2}

    semester_columns.sort(
        key=lambda x: (
            *semester_sort_key(x[0]),
            instrument_order.get(x[1], 99),
            mode_order.get(x[2], 99),
        )
    )
    semesters = sorted({semester for semester, _instrument, _mode in semester_columns}, key=semester_sort_key)

    output_headers = [
        "Category",
        "Prop ID",
        "Rank",
        "Title",
        "PI name",
        "PI Institution",
        "PI email",
        "Open Access",
        "long term?",
        "Student?",
    ]

    telescope_order = ["2m", "1m", "0m4"]

    for semester in semesters:
        for group in telescope_order:
            output_headers.append(f"Hours ({semester} {group})")
            output_headers += [
                f"{col_semester} {instrument} {mode.lower()}"
                for col_semester, instrument, mode in semester_columns
                if col_semester == semester and telescope_group(instrument) == group
            ]
            output_headers += [
                f"Adjusted Hours ({semester} {group})",
                f"{semester} {group} comment",
            ]

    output_headers += [
        "TAC grade",
        "TAC priority",
        "Priority algorithm",
        "Final Priority",
        "Manual Change",
        "Rank (final)",
        "Proposal Code",
        "RR/TC requests valid?",
        "Comment",
        "Memberships",
    ]

    output_rows = []

    for row in input_rows:
        tags = parse_tags(row.get("Tags", ""))

        out = {h: "" for h in output_headers}

        out["Category"] = ""
        out["Prop ID"] = row.get("Proposal ID", "")
        out["Rank"] = row.get("Rank", "")
        out["Title"] = row.get("Title", "")
        out["PI name"] = row.get("PI Name", "")
        out["PI Institution"] = row.get("PI Institution", "")
        out["PI email"] = row.get("PI Email", "")
        out["Open Access"] = checkbox_value("open-access" in tags)
        out["Student?"] = checkbox_value("student" in tags)

        semester_totals = {
            (semester, group): 0.0
            for semester in semesters
            for group in telescope_order
        }

        for col in input_headers:
            m = COLUMN_RE.match(col)
            if not m:
                continue

            semester = m.group("semester")
            aperture = m.group("aperture")
            family = m.group("family")
            camera = m.group("camera")
            mode = m.group("mode")

            instrument = canonical_instrument(aperture, family, camera)
            if instrument is None:
                continue

            value = safe_float(row.get(col, ""))
            out_col = f"{semester} {instrument} {mode.lower()}"

            previous = safe_float(out.get(out_col, ""))
            new_total = previous + value
            out[out_col] = number_or_blank(new_total)

            group = telescope_group(instrument)
            if group is not None:
                semester_totals[(semester, group)] += value

        for (semester, group), total in semester_totals.items():
            out[f"Hours ({semester} {group})"] = number_or_blank(total)

        last_active_semester = ""
        for semester in semesters:
            if any(semester_totals[(semester, group)] for group in telescope_order):
                last_active_semester = semester
        out["long term?"] = last_active_semester if "long-term" in tags else ""

        out["Proposal Code"] = ""
        out["Memberships"] = ""

        output_rows.append(out)

    

    data = StringIO()
    
    writer = csv.DictWriter(data, fieldnames=output_headers, dialect='excel-tab', lineterminator="\n")
    #writer.writeheader()
    writer.writerows(output_rows)
    
    return data.getvalue()


if __name__ == "__main__":
    main()
