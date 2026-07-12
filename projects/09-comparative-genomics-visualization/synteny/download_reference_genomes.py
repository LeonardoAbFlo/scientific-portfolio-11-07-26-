#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import re
import time
import random
import urllib.parse
import urllib.request
from pathlib import Path
from urllib.error import HTTPError, URLError

BASE = Path.home() / "mosquitos/results_p2/data/ref_gb"

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# NCBI request settings
# Set these values if available:
#   export NCBI_EMAIL="you@domain.com"
#   export NCBI_API_KEY="YOUR_NCBI_KEY"
NCBI_EMAIL   = os.getenv("NCBI_EMAIL", "").strip()
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "").strip()
NCBI_TOOL    = os.getenv("NCBI_TOOL", "mosquito_mitogenome_fetcher").strip()

# Request delay depends on API-key availability
_MIN_DELAY = 0.11 if NCBI_API_KEY else 0.35
_last_req_ts = 0.0

def _rate_limit_sleep():
    global _last_req_ts
    now = time.time()
    wait = _MIN_DELAY - (now - _last_req_ts)
    if wait > 0:
        time.sleep(wait)
    _last_req_ts = time.time()

def http_get(url: str, max_retries: int = 6) -> str:
    """Fetch a URL with rate limiting and retries."""
    headers = {
        "User-Agent": f"{NCBI_TOOL}/1.0 ({NCBI_EMAIL or 'no-email-set'})",
        "Accept": "*/*",
    }
    req = urllib.request.Request(url, headers=headers, method="GET")

    backoff = 1.0
    for attempt in range(1, max_retries + 1):
        try:
            _rate_limit_sleep()
            with urllib.request.urlopen(req) as r:
                return r.read().decode("utf-8", errors="replace")

        except HTTPError as e:
            if e.code == 429:
                retry_after = e.headers.get("Retry-After")
                if retry_after:
                    try:
                        wait_s = float(retry_after)
                    except ValueError:
                        wait_s = backoff
                else:
                    wait_s = backoff

                wait_s = wait_s + random.uniform(0, 0.5)
                print(f"[WARN] HTTP 429. Sleeping {wait_s:.2f}s (attempt {attempt}/{max_retries})")
                time.sleep(wait_s)
                backoff = min(backoff * 2.0, 30.0)
                continue

            if e.code in (500, 502, 503, 504):
                wait_s = backoff + random.uniform(0, 0.5)
                print(f"[WARN] HTTP {e.code}. Sleeping {wait_s:.2f}s (attempt {attempt}/{max_retries})")
                time.sleep(wait_s)
                backoff = min(backoff * 2.0, 30.0)
                continue

            raise

        except URLError as e:
            wait_s = backoff + random.uniform(0, 0.5)
            print(f"[WARN] URLError: {e}. Sleeping {wait_s:.2f}s (attempt {attempt}/{max_retries})")
            time.sleep(wait_s)
            backoff = min(backoff * 2.0, 30.0)
            continue

    raise RuntimeError(f"Failed to fetch after {max_retries} retries: {url}")

def _ncbi_params(extra: dict) -> dict:
    p = dict(extra)
    if NCBI_TOOL:
        p["tool"] = NCBI_TOOL
    if NCBI_EMAIL:
        p["email"] = NCBI_EMAIL
    if NCBI_API_KEY:
        p["api_key"] = NCBI_API_KEY
    return p

def esearch_ids(term: str, retmax: int = 20):
    q = _ncbi_params({
        "db": "nuccore",
        "retmode": "json",
        "retmax": str(retmax),
        "term": term,
    })
    url = ESEARCH + "?" + urllib.parse.urlencode(q)
    data = json.loads(http_get(url))
    return data.get("esearchresult", {}).get("idlist", [])

def efetch_gb_by_id(nuccore_id: str) -> str:
    q = _ncbi_params({
        "db": "nuccore",
        "id": nuccore_id,
        "rettype": "gbwithparts",
        "retmode": "text",
    })
    url = EFETCH + "?" + urllib.parse.urlencode(q)
    return http_get(url)

def parse_accession_and_len(gb_text: str):
    acc = None
    m = re.search(r"^VERSION\s+(\S+)", gb_text, flags=re.M)
    if m:
        acc = m.group(1)

    length = None
    m2 = re.search(r"^LOCUS\s+\S+\s+(\d+)\s+bp", gb_text, flags=re.M)
    if m2:
        length = int(m2.group(1))
    return acc, length

def parse_organism(gb_text: str):
    """
    Extract ORGANISM line, return (genus, species_joined_with_underscores).
    If organism has more than 2 tokens (subgenus/subspecies), keep them joined.
    """
    m = re.search(r"^\s{2}ORGANISM\s+(.+)$", gb_text, flags=re.M)
    if not m:
        return None, None
    org = m.group(1).strip()

    # Remove parentheses before building the directory name
    org_clean = re.sub(r"[()]", " ", org)
    org_clean = re.sub(r"\s+", " ", org_clean).strip()
    toks = org_clean.split(" ")
    if not toks:
        return None, None
    genus = toks[0]
    species = "_".join(toks[1:]) if len(toks) > 1 else "unknown"
    return genus, species

def sanitize_name(s: str) -> str:
    s = s.strip().replace(" ", "_")
    s = re.sub(r"[^A-Za-z0-9_]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "unknown"

def is_complete_mitogenome(gb_text: str, min_bp: int = 14000) -> bool:
    acc, length = parse_accession_and_len(gb_text)
    if length is None or length < min_bp:
        return False
    head = gb_text[:5000].lower()
    if "mitochond" not in head:
        return False
    return ("complete genome" in head) or ("complete mitochondrial genome" in head)

def write_record(acc: str, genus: str, species: str, gb_text: str):
    acc0 = acc.split(".")[0]
    genus_s = sanitize_name(genus)
    sp_s = sanitize_name(species)
    folder = BASE / f"{acc0}-{genus_s}_{sp_s}"
    folder.mkdir(parents=True, exist_ok=True)
    out = folder / f"{acc0}.gb"
    out.write_text(gb_text, encoding="utf-8")
    print(f"[OK] {acc} -> {out}")

def write_record_auto(acc: str, gb_text: str):
    g, sp = parse_organism(gb_text)
    if not g:
        g, sp = "Unknown", "unknown"
    write_record(acc, g, sp, gb_text)

def find_complete_mitogenomes(term_variants, n: int = 1, exclude_acc0=None,
                              retmax: int = 80, max_candidates_per_term: int = 20):
    """
    Return list of (acc_version, gb_text) up to n items.
    - exclude_acc0: set of acc without version (e.g., {"OK662580"})
    """
    if exclude_acc0 is None:
        exclude_acc0 = set()

    found = []
    seen = set()

    for term in term_variants:
        ids = esearch_ids(term, retmax=retmax)
        if not ids:
            continue

        for nuccore_id in ids[:max_candidates_per_term]:
            gb = efetch_gb_by_id(nuccore_id)
            if not is_complete_mitogenome(gb):
                continue

            acc, _ = parse_accession_and_len(gb)
            if not acc:
                continue

            acc0 = acc.split(".")[0]
            if acc0 in exclude_acc0:
                continue
            if acc0 in seen:
                continue

            found.append((acc, gb))
            seen.add(acc0)

            if len(found) >= n:
                return found

    return found

def genus_fallback_terms(genus: str):
    return [
        f'{genus}[Organism] mitochondrion[filter] "complete genome"',
        f'"{genus}" mitochondrion[filter] "complete genome"',
        f'{genus}[Organism] "complete mitochondrial genome"',
        f'"{genus}" "complete mitochondrial genome"',
    ]

def main():
    BASE.mkdir(parents=True, exist_ok=True)

    # Target species
    targets = [
        # Accepted Anopheles name variants
        ("Anopheles", "nyssorhynchus_triannulatus", [
            'Anopheles triannulatus[Organism] mitochondrion[filter] "complete genome"',
            '"Anopheles triannulatus" mitochondrion[filter] "complete genome"',
            '"Anopheles (Nyssorhynchus) triannulatus"[Organism] mitochondrion[filter] "complete genome"',
            '"Nyssorhynchus triannulatus" mitochondrion[filter] "complete genome"',
            'Anopheles triannulatus[Organism] "complete mitochondrial genome"',
        ]),

        # Culex quinquefasciatus
        ("Culex", "quinquefasciatus", [
            'Culex quinquefasciatus[Organism] mitochondrion[filter] "complete genome"',
            '"Culex quinquefasciatus" mitochondrion[filter] "complete genome"',
            'Culex quinquefasciatus[Organism] "complete mitochondrial genome"',
        ]),

        # Culex iolambdis
        ("Culex", "iolambdis", [
            'Culex iolambdis[Organism] mitochondrion[filter] "complete genome"',
            '"Culex iolambdis" mitochondrion[filter] "complete genome"',
            'Culex iolambdis[Organism] "complete mitochondrial genome"',
        ]),

        # Psorophora ferox
        ("Psorophora", "ferox", [
            'Psorophora ferox[Organism] mitochondrion[filter] "complete genome"',
            '"Psorophora ferox" mitochondrion[filter] "complete genome"',
            'Psorophora ferox[Organism] "complete mitochondrial genome"',
        ]),
    ]

    downloaded_by_genus = {}   # genus -> set(acc0)
    genus_had_miss = {}        # genus -> bool

    def _mark_download(genus: str, acc: str):
        acc0 = acc.split(".")[0]
        downloaded_by_genus.setdefault(genus, set()).add(acc0)

    for genus, species, terms in targets:
        label = f"{genus} {species.replace('_', ' ')}"
        print(f"[INFO] Searching complete mitogenome for: {label}")

        exclude = downloaded_by_genus.get(genus, set())
        hit = find_complete_mitogenomes(terms, n=1, exclude_acc0=exclude, retmax=80, max_candidates_per_term=20)

        if hit:
            acc, gb = hit[0]
            write_record(acc, genus, species, gb)
            _mark_download(genus, acc)
            genus_had_miss.setdefault(genus, False)
        else:
            print(f"[MISS] No complete mitogenome found for {label} via NCBI search terms.")
            genus_had_miss[genus] = True

    # If a target is missing, download two complete mitogenomes for that genus
    for genus, had_miss in genus_had_miss.items():
        if not had_miss:
            continue

        current = downloaded_by_genus.get(genus, set())
        need = 2 - len(current)
        if need <= 0:
            continue

        print(f"[INFO] Fallback: need {need} more complete mitogenome(s) for genus {genus}")
        extra_terms = genus_fallback_terms(genus)
        extras = find_complete_mitogenomes(extra_terms, n=need, exclude_acc0=current, retmax=120, max_candidates_per_term=35)

        if not extras:
            print(f"[MISS] Fallback failed: no additional complete mitogenomes found for genus {genus}.")
            continue

        for acc, gb in extras:
            # Use the organism name from the GenBank record
            write_record_auto(acc, gb)
            _mark_download(genus, acc)

if __name__ == "__main__":
    main()
