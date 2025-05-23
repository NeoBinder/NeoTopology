import re
from collections import defaultdict, deque

import numpy as np

from ttk.math import RotationMatrix


def split_pdb_content(content):
    if not isinstance(content, list):
        content = content.split("\n")
    content = deque(content)
    headers = deque([])
    while len(content):
        line = content.popleft()
        if line[0:6] in ("ATOM  ", "HETATM", "MODEL "):
            content.appendleft(line)
            break
        else:
            headers.append(line)
    return list(headers), list(content)


def get_references(inl):
    # REMARK   1 REFERENCE 1                                                  1CSE  11
    # REMARK   1  AUTH   W.BODE,E.PAPAMOKOS,D.MUSIL                           1CSE  12
    references = []
    actref = ""
    for l in inl:
        if re.search(r"\AREMARK   1", l):
            if re.search(r"\AREMARK   1 REFERENCE", l):
                if actref != "":
                    actref = re.sub(r"\s\s+", " ", actref)
                    if actref != " ":
                        references.append(actref)
                    actref = ""
            else:
                actref += l[19:72].lower()

    if actref != "":
        actref = re.sub(r"\s\s+", " ", actref)
        if actref != " ":
            references.append(actref)
    return references


def get_symmetry_matrices(header):
    matrices = defaultdict(list)
    for line in header:
        if re.search("SMTRY", line):
            line = line.strip().split()
            matrices[int(line[3])].append(line[4:])
    for k, v in matrices.items():
        rmatrix = np.asarray(v, dtype=np.float64)
        # convert to nanometer
        rmatrix[:3, 3] /= 10
        matrices[k] = RotationMatrix.from_matrix(rmatrix)
    return matrices


def get_journal(inl):
    # JRNL        AUTH   L.CHEN,M.DOI,F.S.MATHEWS,A.Y.CHISTOSERDOV,           2BBK   7
    journal = ""
    for l in inl:
        if re.search(r"\AJRNL", l):
            journal += l[19:72].lower()
    journal = re.sub(r"\s\s+", " ", journal)
    return journal


def chop_end_codes(line):
    """Chops lines ending with  '     1CSA  14' and the like (PRIVATE)."""
    return re.sub(r"\s\s\s\s+[\w]{4}.\s+\d*\Z", "", line)


def chop_end_misc(line):
    """Chops lines ending with  '     14-JUL-97  1CSA' and the like (PRIVATE)."""
    return re.sub(r"\s+\d\d-\w\w\w-\d\d\s+[1-9][0-9A-Z]{3}\s*\Z", "", line)


def nice_case(line):
    """Make A Lowercase String With Capitals (PRIVATE)."""
    line_lower = line.lower()
    s = ""
    i = 0
    nextCap = 1
    while i < len(line_lower):
        c = line_lower[i]
        if c >= "a" and c <= "z" and nextCap:
            c = c.upper()
            nextCap = 0
        elif c in " .,;:\t-_":
            nextCap = 1
        s += c
        i += 1
    return s


# bring dates to format: 1909-01-08
def format_date(pdb_date):
    """Convert dates from DD-Mon-YY to YYYY-MM-DD format (PRIVATE)."""
    date = ""
    year = int(pdb_date[7:])
    if year < 50:
        century = 2000
    else:
        century = 1900
    date = str(century + year) + "-"
    all_months = [
        "xxx",
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ]
    month = str(all_months.index(pdb_date[3:6]))
    if len(month) == 1:
        month = "0" + month
    date = date + month + "-" + pdb_date[:2]
    return date


def parse_remark_465(line):
    """Parse missing residue remarks.
    Returns a dictionary describing the missing residue.
    The specification for REMARK 465 at
    http://www.wwpdb.org/documentation/file-format-content/format33/remarks2.html#REMARK%20465
    only gives templates, but does not say they have to be followed.
    So we assume that not all pdb-files with a REMARK 465 can be understood.
    Returns a dictionary with the following keys:
    "model", "res_name", "chain", "ssseq", "insertion"
    """
    if line:
        # Note that line has been stripped.
        assert line[0] != " " and line[-1] not in "\n ", "line has to be stripped"
    pattern = re.compile(
        r"""
        (\d+\s[\sA-Z][\sA-Z][A-Z] |   # Either model number + residue name
            [A-Z]{1,3})               # Or only residue name with 1 (RNA) to 3 letters
        \s ([A-Za-z0-9])              # A single character chain
        \s+(-?\d+[A-Za-z]?)$          # Residue number: A digit followed by an optional
                                      # insertion code (Hetero-flags make no sense in
                                      # context with missing res)
        """,
        re.VERBOSE,
    )
    match = pattern.match(line)
    if match is None:
        return None
    residue = {}
    if " " in match.group(1):
        model, residue["res_name"] = match.group(1).split()
        residue["model"] = int(model)
    else:
        residue["model"] = None
        residue["res_name"] = match.group(1)
    residue["chain"] = match.group(2)
    try:
        residue["ssseq"] = int(match.group(3))
    except ValueError:
        residue["insertion"] = match.group(3)[-1]
        residue["ssseq"] = int(match.group(3)[:-1])
    else:
        residue["insertion"] = None
    return residue


def parse_pdb_header_list(header):
    # database fields
    pdbh_dict = {
        "name": "",
        "head": "",
        "idcode": "",
        "deposition_date": "1909-01-08",
        "release_date": "1909-01-08",
        "structure_method": "unknown",
        "resolution": None,
        "structure_reference": "unknown",
        "journal_reference": "unknown",
        "periodic_box": None,
        "author": "",
        "compound": {"1": {"misc": ""}},
        "source": {"1": {"misc": ""}},
        "has_missing_residues": False,
        "missing_residues": [],
    }

    pdbh_dict["structure_reference"] = get_references(header)
    pdbh_dict["journal_reference"] = get_journal(header)
    pdbh_dict["symmetry_rotation_matrix"] = get_symmetry_matrices(header)
    comp_molid = "1"
    last_comp_key = "misc"
    last_src_key = "misc"
    for hh in header:
        h = re.sub(r"[\s\n\r]*\Z", "", hh)  # chop linebreaks off
        # key=re.sub("\s.+\s*","",h)
        key = h[:6].strip()
        # tail=re.sub("\A\w+\s+\d*\s*","",h)
        tail = h[10:].strip()
        # print("%s:%s" % (key, tail)

        # From here, all the keys from the header are being parsed
        if key == "TITLE":
            name = chop_end_codes(tail).lower()
            pdbh_dict["name"] = " ".join([pdbh_dict["name"], name]).strip()
        elif key == "CRYST1":
            pdbh_dict["periodic_box"] = hh
        elif key == "HEADER":
            rr = re.search(r"\d\d-\w\w\w-\d\d", tail)
            if rr is not None:
                pdbh_dict["deposition_date"] = format_date(nice_case(rr.group()))
            rr = re.search(r"\s+([1-9][0-9A-Z]{3})\s*\Z", tail)
            if rr is not None:
                pdbh_dict["idcode"] = rr.group(1)
            head = chop_end_misc(tail).lower()
            pdbh_dict["head"] = head
        elif key == "COMPND":
            tt = re.sub(r"\;\s*\Z", "", chop_end_codes(tail)).lower()
            # look for E.C. numbers in COMPND lines
            rec = re.search(r"\d+\.\d+\.\d+\.\d+", tt)
            if rec:
                pdbh_dict["compound"][comp_molid]["ec_number"] = rec.group()
                tt = re.sub(r"\((e\.c\.)*\d+\.\d+\.\d+\.\d+\)", "", tt)
            tok = tt.split(":")
            if len(tok) >= 2:
                ckey = tok[0]
                cval = re.sub(r"\A\s*", "", tok[1])
                if ckey == "mol_id":
                    pdbh_dict["compound"][cval] = {"misc": ""}
                    comp_molid = cval
                    last_comp_key = "misc"
                else:
                    pdbh_dict["compound"][comp_molid][ckey] = cval
                    last_comp_key = ckey
            else:
                pdbh_dict["compound"][comp_molid][last_comp_key] += tok[0] + " "
        elif key == "SOURCE":
            tt = re.sub(r"\;\s*\Z", "", chop_end_codes(tail)).lower()
            tok = tt.split(":")
            # print(tok)
            if len(tok) >= 2:
                ckey = tok[0]
                cval = re.sub(r"\A\s*", "", tok[1])
                if ckey == "mol_id":
                    pdbh_dict["source"][cval] = {"misc": ""}
                    comp_molid = cval
                    last_src_key = "misc"
                else:
                    pdbh_dict["source"][comp_molid][ckey] = cval
                    last_src_key = ckey
            else:
                pdbh_dict["source"][comp_molid][last_src_key] += tok[0] + " "
        elif key == "KEYWDS":
            kwd = chop_end_codes(tail).lower()
            if "keywords" in pdbh_dict:
                pdbh_dict["keywords"] += " " + kwd
            else:
                pdbh_dict["keywords"] = kwd
        elif key == "EXPDTA":
            expd = chop_end_codes(tail)
            # chop junk at end of lines for some structures
            expd = re.sub(r"\s\s\s\s\s\s\s.*\Z", "", expd)
            # if re.search('\Anmr',expd,re.IGNORECASE): expd='nmr'
            # if re.search('x-ray diffraction',expd,re.IGNORECASE): expd='x-ray diffraction'
            pdbh_dict["structure_method"] = expd.lower()
        elif key == "CAVEAT":
            # make Annotation entries out of these!!!
            pass
        elif key == "REVDAT":
            rr = re.search(r"\d\d-\w\w\w-\d\d", tail)
            if rr is not None:
                pdbh_dict["release_date"] = format_date(nice_case(rr.group()))
        elif key == "JRNL":
            # print("%s:%s" % (key, tail))
            if "journal" in pdbh_dict:
                pdbh_dict["journal"] += tail
            else:
                pdbh_dict["journal"] = tail
        elif key == "AUTHOR":
            auth = nice_case(chop_end_codes(tail))
            if "author" in pdbh_dict:
                pdbh_dict["author"] += auth
            else:
                pdbh_dict["author"] = auth
        elif key == "REMARK":
            if re.search("REMARK   2 RESOLUTION.", hh):
                r = chop_end_codes(re.sub("REMARK   2 RESOLUTION.", "", hh))
                r = re.sub(r"\s+ANGSTROM.*", "", r)
                try:
                    pdbh_dict["resolution"] = float(r)
                except ValueError:
                    # print('nonstandard resolution %r' % r)
                    pdbh_dict["resolution"] = None
            elif hh.startswith("REMARK 465"):
                if tail:
                    pdbh_dict["has_missing_residues"] = True
                    missing_res_info = parse_remark_465(tail)
                    if missing_res_info:
                        pdbh_dict["missing_residues"].append(missing_res_info)
            elif hh.startswith("REMARK  99 ASTRAL"):
                if tail:
                    remark_99_keyval = tail.replace("ASTRAL ", "").split(": ")
                    if type(remark_99_keyval) == list and len(remark_99_keyval) == 2:
                        if "astral" not in pdbh_dict:
                            pdbh_dict["astral"] = {
                                remark_99_keyval[0]: remark_99_keyval[1]
                            }
                        else:
                            pdbh_dict["astral"][remark_99_keyval[0]] = remark_99_keyval[
                                1
                            ]
        else:
            # print(key)
            pass
    if pdbh_dict["structure_method"] == "unknown":
        res = pdbh_dict["resolution"]
        if res is not None and res > 0.0:
            pdbh_dict["structure_method"] = "x-ray diffraction"
    return pdbh_dict
