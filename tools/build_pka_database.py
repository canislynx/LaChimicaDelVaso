#!/usr/bin/env python3
"""Genera pka-database.js a partire dal foglio "Database" di CurTiPot (curtipot_i.xlsm).

Uso (dalla cartella del progetto):
    python3 -I tools/build_pka_database.py Source/curtipot_i.xlsm pka-database.js

Legge SOLO i dati: il file .xlsm è uno zip di XML e qui si usa la libreria standard,
non si apre in Excel e non si esegue nessuna macro.

Struttura del foglio (colonne): B = numero, C = nome, D = carica della forma
completamente deprotonata, E..J = pKa1..pKa6, K = °C, L = forza ionica, M = formula.
Sezione indicatori (righe 280-288): O = pKa, P/R = colore (formato Excel BGR) sotto/sopra il pKa.

Interventi sui dati (tutti elencati qui, niente è modificato di nascosto):
  * si scartano le righe senza nome o senza pKa (segnaposto del foglio);
  * le voci della sezione "frequently used" ripetute in ordine alfabetico si fondono
    (stessa formula, carica e pKa); il nome alternativo resta cercabile;
  * "A = B = C" nei nomi: il primo è il nome, gli altri sono sinonimi cercabili;
  * pKa ordinati in modo crescente (Hexamethylenediamine li ha invertiti);
  * formula con "CI" scritto con la I maiuscola al posto di "Cl" (es. C6H6CIN);
  * formule con carica esplicita ("HCrO4-", "HSeO4-"): la carica indicata nel foglio come
    "forma deprotonata" è incoerente, la si ricava dalla formula;
  * refusi evidenti nei nomi (TYPO_FIXES).
"""
import re
import sys
import json
import zipfile
import xml.etree.ElementTree as ET

M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
SHEET = "xl/worksheets/sheet12.xml"  # foglio "Database" (xl/workbook.xml, rId12)

TYPO_FIXES = {
    "4-Animobenzoic acid": "4-Aminobenzoic acid",
    "2,5-Dmethylpyridine": "2,5-Dimethylpyridine",
    "3.5-Dimethylpyridine": "3,5-Dimethylpyridine",
    "Hydroazoic": "Hydrazoic acid",
    "Tris(hydroxymethyl)- aminomethane (tris)": "Tris(hydroxymethyl)aminomethane (tris)",
}


def read_sheet(path):
    z = zipfile.ZipFile(path)
    shared = []
    for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(M + "si"):
        shared.append("".join(t.text or "" for t in si.iter(M + "t")))
    rows = {}
    for c in ET.fromstring(z.read(SHEET)).iter(M + "c"):
        v = c.find(M + "v")
        if v is None:
            continue
        val = shared[int(v.text)] if c.get("t") == "s" else v.text
        m = re.match(r"([A-Z]+)(\d+)", c.get("r"))
        rows.setdefault(int(m.group(2)), {})[m.group(1)] = val
    return rows


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def clean(s):
    return re.sub(r"\s+", " ", s).strip()


def sort_key(name):
    s = name
    while True:
        t = re.sub(r"^\d[\d,.'\s-]*", "", s)
        t = re.sub(r"^(sec|tert|meso|[A-Za-z])-", "", t)
        t = re.sub(r"^\([^)]*\)\s*", "", t)
        if t == s:
            break
        s = t
    return s.lower()


def bgr(n):
    n = int(float(n))
    return [n & 255, (n >> 8) & 255, (n >> 16) & 255]


def main(src, dst):
    rows = read_sheet(src)
    subs, seen = [], {}
    for r in sorted(rows):
        d = rows[r]
        if r < 22 or r >= 279 or not d.get("C") or "E" not in d:
            continue
        pks = [num(d[k]) for k in "EFGHIJ" if k in d]
        if None in pks:  # riga di note (es. "Insert new lines anywhere...")
            continue
        name = clean(d["C"])
        parts = [clean(p) for p in name.split(" = ")]
        main_name = TYPO_FIXES.get(parts[0], parts[0])
        alts = parts[1:]
        formula = clean(d.get("M", "")).replace("CI", "Cl")
        z = int(float(d["D"])) if d.get("D") not in (None, "") else 0
        m_charge = re.search(r"([+-])$", formula)
        if m_charge:  # formula di una specie carica: la più protonata ha quella carica
            z = (1 if m_charge.group(1) == "+" else -1) - len(pks)
        pks = sorted(round(p, 3) for p in pks)
        t = num(d.get("K"))
        ion = num(d.get("L"))
        entry = {"n": main_name, "f": formula, "z": z, "pk": pks}
        if alts:
            entry["a"] = alts
        if t is not None:
            entry["t"] = int(t) if t == int(t) else t
        if ion is not None:
            entry["i"] = int(ion) if ion == int(ion) else ion
        key = (re.sub(r"\s", "", formula).lower(), z, tuple(pks))
        if key in seen:
            prev = seen[key]
            if prev["n"] != main_name and main_name not in prev.get("a", []):
                prev.setdefault("a", []).append(main_name)
            continue
        seen[key] = entry
        subs.append(entry)

    # ordine alfabetico ignorando locanti e prefissi (2-, m-, sec-...)
    subs.sort(key=lambda e: (sort_key(e["n"]), e["n"]))

    inds = []
    for r in range(280, 289):
        d = rows.get(r)
        if not d or "O" not in d:
            continue
        inds.append({
            "n": clean(d["C"]),
            "pk": round(float(d["O"]), 2),
            "z": int(float(d["D"])),
            "ac": bgr(d["P"]),
            "bc": bgr(d["R"]),
        })
    inds.sort(key=lambda e: e["pk"])

    out = []
    out.append("/* pKa e sostanze: dal foglio \"Database\" di CurTiPot (versione 4.4.0, gennaio 2024),")
    out.append("   di Ivano G. R. Gutz, Instituto de Quimica, Universidade de Sao Paulo, Brasile.")
    out.append("   https://www.iq.usp.br/gutz")
    out.append("   File GENERATO da tools/build_pka_database.py: non modificarlo a mano.")
    out.append("   Chiavi: n nome, a sinonimi, f formula (specie neutra), z carica della forma completamente")
    out.append("   deprotonata, pk pKa crescenti, t temperatura in C, i forza ionica; indicatori: ac/bc colore RGB")
    out.append("   della forma acida/basica. La maggior parte dei valori e' a 25 C e forza ionica 0. */")
    out.append("window.CURTIPOT_DB = {")
    out.append('  meta: ' + json.dumps({
        "name": "CurTiPot", "version": "4.4.0", "date": "2024-01",
        "author": "Ivano G. R. Gutz", "affiliation": "Instituto de Química, Universidade de São Paulo, Brasile",
        "url": "https://www.iq.usp.br/gutz"}, ensure_ascii=False) + ",")
    out.append("  substances: [")
    out.extend("    " + json.dumps(e, ensure_ascii=False, separators=(",", ":")) + "," for e in subs)
    out.append("  ],")
    out.append("  indicators: [")
    out.extend("    " + json.dumps(e, ensure_ascii=False, separators=(",", ":")) + "," for e in inds)
    out.append("  ]")
    out.append("};")
    with open(dst, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")
    print("sostanze:", len(subs), " indicatori:", len(inds), " ->", dst)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
