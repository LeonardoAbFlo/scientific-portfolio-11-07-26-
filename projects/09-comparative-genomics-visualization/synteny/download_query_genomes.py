#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

BASE = Path.home() / "mosquitos/results_p2/data/ref_gb"

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def http_get(url: str) -> str:
    with urllib.request.urlopen(url) as r:
        return r.read().decode("utf-8", errors="replace")

def esearch_ids(term: str, retmax: int = 20):
    q = {
        "db": "nuccore",
        "retmode": "json",
        "retmax": str(retmax),
        "term": term,
    }
    url = ESEARCH + "?" + urllib.parse.urlencode(q)
    data = json.loads(http_get(url))
    return data.get("esearchresult", {}).get("idlist", [])

def efetch_gb_by_id(nuccore_id: str) -> str:
    q = {
        "db": "nuccore",
        "id": nuccore_id,
        "rettype": "gbwithparts",
        "retmode": "text",
    }
    url = EFETCH + "?" + urllib.parse.urlencode(q)
    return http_get(url)

def parse_accession_and_len(gb_text: str):
    # Read the accession version
    acc = None
    m = re.search(r"^VERSION\s+(\S+)", gb_text, flags=re.M)
    if m:
        acc = m.group(1)
    # Read sequence length from the LOCUS line
    length = None
    m2 = re.search(r"^LOCUS\s+\S+\s+(\d+)\s+bp", gb_text, flags=re.M)
    if m2:
        length = int(m2.group(1))
    return acc, length

def is_complete_mitogenome(gb_text: str, min_bp: int = 14000) -> bool:
    acc, length = parse_accession_and_len(gb_text)
    if length is None:
        return False
    # Check sequence length and annotation keywords
    if length < min_bp:
        return False
    head = gb_text[:4000].lower()
    return ("mitochond" in head) and ("complete genome" in head or "complete mitochondrial genome" in head)

def write_record(acc: str, genus: str, species: str, gb_text: str):
    folder = BASE / f"{acc.split('.')[0]}-{genus}_{species}"
    folder.mkdir(parents=True, exist_ok=True)
    out = folder / f"{acc.split('.')[0]}.gb"
    out.write_text(gb_text, encoding="utf-8")
    print(f"[OK] {acc} -> {out}")

def fetch_by_accession(acc: str) -> str:
    # Download the record
    q = {
        "db": "nuccore",
        "id": acc,
        "rettype": "gbwithparts",
        "retmode": "text",
    }
    url = EFETCH + "?" + urllib.parse.urlencode(q)
    return http_get(url)

def find_best_complete_mitogenome(term_variants):
    # Return accession and GenBank text when found
    for term in term_variants:
        ids = esearch_ids(term, retmax=50)
        for nuccore_id in ids:
            gb = efetch_gb_by_id(nuccore_id)
            if is_complete_mitogenome(gb):
                acc, _ = parse_accession_and_len(gb)
                if acc:
                    return acc, gb
    return None, None

def main():
    BASE.mkdir(parents=True, exist_ok=True)

    # Download known accessions
    direct = [
        # Existing files can be downloaded again
        ("OK662580", "Coquillettidia", "nigricans"),
        ("MK575482", "Limatus", "flavisetosus"),
        ("OM275430", "Ochlerotatus", "serratus"),   # stored as Ochlerotatus_* in your convention
        ("OM275429", "Ochlerotatus", "scapularis"), # proxy nearest to fulvus
        ("MK575483", "Mansonia", "amazonensis"),    # proxy nearest to humeralis
        ("MN342085", "Mansonia", "uniformis"),
    ]

    for acc, genus, sp in direct:
        gb = fetch_by_accession(acc)
        accv, _ = parse_accession_and_len(gb)
        if not accv:
            print(f"[WARN] Could not parse accession for {acc}")
            continue
        write_record(accv, genus, sp, gb)

    # Search for Coquillettidia venezuelensis
    cq_terms = [
        'Coquillettidia venezuelensis[Organism] mitochondrion[filter] "complete genome"',
        '"Coquillettidia venezuelensis" mitochondrion[filter] "complete genome"',
        'Coquillettidia venezuelensis[Organism] "complete mitochondrial genome"',
    ]
    acc, gb = find_best_complete_mitogenome(cq_terms)
    if acc:
        write_record(acc, "Coquillettidia", "venezuelensis", gb)
    else:
        print("[MISS] No complete mitogenome found for Coquillettidia venezuelensis via NCBI search terms.")

    # Search both Ochlerotatus and Aedes names for O. fulvus
    of_terms = [
        'Ochlerotatus fulvus[Organism] mitochondrion[filter] "complete genome"',
        'Aedes fulvus[Organism] mitochondrion[filter] "complete genome"',
        '"Ochlerotatus fulvus" "complete mitochondrial genome"',
        '"Aedes fulvus" "complete mitochondrial genome"',
    ]
    acc, gb = find_best_complete_mitogenome(of_terms)
    if acc:
        write_record(acc, "Ochlerotatus", "fulvus", gb)
    else:
        print("[MISS] No complete mitogenome found for Ochlerotatus/Aedes fulvus (keeping scapularis+serratus as nearest proxies).")

if __name__ == "__main__":
    main()
