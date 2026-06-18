import requests
from collections import defaultdict
import json
import time
from pathlib import Path
import pandas as pd
from collections import Counter

def check_relationship(df, gse1, gse2):
    """
    Pull rows associated with two GSE accessions from a dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with a 'gse' column.
    gse1, gse2 : str
        GSE accessions to retrieve.

    Returns
    -------
    pandas.DataFrame
        Rows matching gse1 or gse2.
    """
    wanted = [gse1, gse2]

    out = df[df["gse"].isin(wanted)].copy()

    print(out)

    return out

def get_gse_soft_metadata(gse: str) -> dict:
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    params = {
        "acc": gse,
        "targ": "self",
        "form": "text",
        "view": "full",
    }

    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()

    metadata = defaultdict(list)

    for line in r.text.splitlines():
        if line.startswith("!Series_"):
            key, value = line.split(" = ", 1)
            key = key.removeprefix("!Series_")
            metadata[key].append(value)

    # collapse singleton fields, keep repeated fields as lists
    return {
        key: values[0] if len(values) == 1 else values
        for key, values in metadata.items()
    }

def read_completed_gses(jsonl_path) -> set[str]:
    jsonl_path = Path(jsonl_path)

    if not jsonl_path.exists():
        return set()

    completed = set()

    with open(jsonl_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                record = json.loads(line)
                completed.add(record["gse"])

    return completed

def get_data(gse_metadata, completed, failed_jsonl, output_jsonl):
    for i, gse in enumerate(gse_metadata, start=1):
        if gse in completed:
            continue

        print(f"[{i}/{len(gse_metadata)}] Fetching {gse}...")

        try:
            soft_metadata = get_gse_soft_metadata(gse)

            record = {
                "gse": gse,
                "original_metadata": gse_metadata[gse],
                "geo_soft": soft_metadata,
            }

            with open(output_jsonl, "a", encoding="utf-8") as f:
                f.write(json.dumps(record, ensure_ascii=False) + "\n")
                f.flush()

        except Exception as e:
            print(f"  Failed: {e}")

            fail_record = {
                "gse": gse,
                "error": str(e),
            }

            with open(failed_jsonl, "a", encoding="utf-8") as f:
                f.write(json.dumps(fail_record, ensure_ascii=False) + "\n")
                f.flush()

        time.sleep(0.34)

def load_original_metadata(input_json):
    with open(input_json, "r", encoding="utf-8") as f:
        return json.load(f)

def collapse_list(x, sep="; "):
    """Turn list-valued cells into strings so the dataframe is CSV-friendly."""
    if isinstance(x, list):
        return sep.join(map(str, x))
    return x


def geo_jsonl_to_dataframe(jsonl_path):
    jsonl_path = Path(jsonl_path)

    rows = []

    with open(jsonl_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue

            record = json.loads(line)

            gse = record["gse"]
            original = record.get("original_metadata", {})
            geo_soft = record.get("geo_soft", {})

            row = {
                "gse": gse,

                # useful original fields
                "original_title": original.get("Title"),
                "original_experiment_type": original.get("Experiment type"),
                "original_organism": original.get("Organism"),
                "original_sra": original.get("SRA"),
                "original_n_samples": len(original.get("Samples", [])),

                # useful fetched GEO fields
                "title": geo_soft.get("title"),
                "geo_accession": geo_soft.get("geo_accession"),
                "status": geo_soft.get("status"),
                "submission_date": geo_soft.get("submission_date"),
                "last_update_date": geo_soft.get("last_update_date"),
                "pubmed_id": geo_soft.get("pubmed_id"),
                "contributor": geo_soft.get("contributor"),
                "contact_name": geo_soft.get("contact_name"),
                "contact_email": geo_soft.get("contact_email"),
                "contact_institute": geo_soft.get("contact_institute"),
                "contact_country": geo_soft.get("contact_country"),
                "type": geo_soft.get("type"),
                "platform_id": geo_soft.get("platform_id"),
                "sample_id": geo_soft.get("sample_id"),
                "n_geo_samples": len(geo_soft.get("sample_id", [])),
                "relation": geo_soft.get("relation"),
                "supplementary_file": geo_soft.get("supplementary_file"),
            }

            rows.append(row)

    df = pd.DataFrame(rows)

    # Collapse list-valued columns into semicolon-separated strings
    for col in df.columns:
        df[col] = df[col].map(collapse_list)

    return df

def main():
    input_json = Path("../metadata/metadata_human.json")
    gse_metadata = load_original_metadata(input_json)
    total_entries = len(gse_metadata)

    output_jsonl = Path("../metadata/geo_soft_metadata_human.jsonl")
    failed_jsonl = Path("../metadata/failed_gses_human.jsonl")

    completed = read_completed_gses(output_jsonl)
    failed_before = read_completed_gses(failed_jsonl)

    print(f"Already completed: {len(completed)}")
    print(f"Previously failed: {len(failed_before)}")

    if len(completed) < total_entries:
        get_data(gse_metadata, completed, failed_jsonl, output_jsonl)

    df = geo_jsonl_to_dataframe("../metadata/geo_soft_metadata_human.jsonl")
    # df.to_csv("../metadata/human_geo_soft.csv", index=False)

    info_dict = dict()
    all_pmids = list()
    all_keys = set()

    for _, row in df.iterrows():
        gse = row["gse"]
        pmid = row.get("pubmed_id")

        info_dict[gse] = dict()

        if pd.notna(pmid):
            info_dict[gse]["pmid"] = pmid
            all_pmids.append(pmid)

        try:
            identifiers = row["relation"].split("; ")
        except AttributeError:
            print(f"{gse} has no relation")
            continue
        for identifier in identifiers:
            name, contents = [i.strip() for i in identifier.split(": ")]
            entries = info_dict[gse].get(name, [])
            entries.append(contents)
            info_dict[gse][name] = entries

    print(all_keys)
    # Write the data to a file
    rows = []

    for gse, inner_dict in info_dict.items():
        rows.append({
            "gse": gse,
            "pmid": inner_dict.get("pmid"),
            "subseries": inner_dict.get('SubSeries of'),
            "superseries": inner_dict.get('SuperSeries of'),
            "affiliation": inner_dict.get('Affiliated with'),
            "BioProject": inner_dict.get("BioProject"),
            "SRA": inner_dict.get("SRA"),
        })

    outdf = pd.DataFrame(rows)
    outdf.to_csv("../metadata/gse_pmid_subseries_of.tsv", sep="\t", index=False)



if __name__ == "__main__":
    #main()
    data = pd.read_csv("../metadata/gse_pmid_subseries_of.tsv", delimiter='\t')
    check_relationship(data, "GSE101157", "GSE101167")