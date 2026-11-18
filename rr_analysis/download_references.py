"""Download reference material for precise RR typing.

The RR3 GenBank record KR856283 is the primary structural reference because it
is a finished and annotated sequence. The Wang et al. 2023 representative runs
are listed for completeness. PacBio RS II runs from 2013 are stored as tar.gz
archives of legacy HDF5 files and are not downloaded by default. Illumina runs
can be downloaded and assembled separately.
"""
import argparse
import os
import sys
import time
import urllib.request

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

RR3_GENBANK = "KR856283"

# Representative strains read from Figure 1A and verified against the ENA
# Portal API. PacBio entries point to submitted tar.gz archives; Illumina
# entries point to paired FASTQ files.
REPRESENTATIVES = [
    ("RR1", "YZU0623", "ERR11897013", "pacbio", ""),
    ("RR2", "YZU1053", "ERR11011713", "illumina",
     "vol1/fastq/ERR110/013/ERR11011713/ERR11011713"),
    ("RR3", "YZU0496", "ERR11897011", "pacbio", ""),
    ("RR4", "YZU3081", "ERR11897014", "pacbio", ""),
    ("RR5", "YZU0530", "ERR11897012", "pacbio", ""),
    ("RR6", "YZU0651", "ERR11011855", "pacbio", ""),
    ("RR7", "YZU0432", "ERR11897007", "pacbio", ""),
    ("RR8", "YZU0484", "ERR12019130", "illumina",
     "vol1/fastq/ERR120/013/ERR12019130/ERR12019130"),
]


def download(url, destination):
    if os.path.exists(destination) and os.path.getsize(destination) > 0:
        print("exists", destination)
        return
    tmp = destination + ".tmp"
    print("downloading", url)
    urllib.request.urlretrieve(url, tmp)
    os.replace(tmp, destination)


def fetch_rr3(out_dir):
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           "?db=nuccore&id=" + RR3_GENBANK + "&rettype=gb&retmode=text")
    download(url, os.path.join(out_dir, RR3_GENBANK + ".gb"))


def fetch_illumina(out_dir):
    base = "https://ftp.sra.ebi.ac.uk/"
    for rr_type, strain, run, platform, path in REPRESENTATIVES:
        if platform != "illumina":
            continue
        target = os.path.join(out_dir, run)
        os.makedirs(target, exist_ok=True)
        for mate in ("1", "2"):
            url = f"{base}{path}_{mate}.fastq.gz"
            download(url, os.path.join(target, f"{run}_{mate}.fastq.gz"))
        time.sleep(1)


def fetch_pacbio(out_dir):
    for rr_type, strain, run, platform, path in REPRESENTATIVES:
        if platform != "pacbio":
            continue
        target = os.path.join(out_dir, run)
        os.makedirs(target, exist_ok=True)
        url = (f"https://ftp.sra.ebi.ac.uk/vol1/run/{run[:7]}/{run}/"
               f"{strain}.tar.gz")
        download(url, os.path.join(target, f"{strain}.tar.gz"))
        time.sleep(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--with-illumina", action="store_true")
    ap.add_argument("--with-pacbio", action="store_true",
                    help="legacy 2013 PacBio RS II tar.gz archives")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    fetch_rr3(args.out_dir)
    if args.with_illumina:
        fetch_illumina(args.out_dir)
    if args.with_pacbio:
        fetch_pacbio(args.out_dir)


if __name__ == "__main__":
    main()
