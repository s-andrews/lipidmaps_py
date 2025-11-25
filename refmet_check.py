import csv
import requests
import time
import logging
from typing import List, Dict, Tuple
from io import StringIO

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

REFMET_API_URL = "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_minID.php"
BATCH_SIZE = 1000
DELAY_BETWEEN_BATCHES = 1  # seconds

def load_lm_data(filename: str) -> List[Dict[str, str]]:
    """Load all LM IDs and their associated names as a list of dicts."""
    logger.info(f"Loading data from {filename}")
    lm_entries = []
    with open(filename, mode='r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        if 'name' not in reader.fieldnames or 'lm_id' not in reader.fieldnames:
            raise ValueError(f"CSV must have 'name' and 'lm_id' columns. Found: {reader.fieldnames}")
        for row_num, row in enumerate(reader, start=2):
            name = (row.get('name') or '').strip()
            lm_id = (row.get('lm_id') or '').strip()
            if not name or not lm_id:
                continue
            lm_entries.append({"lm_id": lm_id, "name": name})
    logger.info(f"Loaded {len(lm_entries)} LM ID/name pairs")
    return lm_entries

def query_refmet(names: List[str]) -> List[Dict]:
    """Query RefMet API (TSV response)."""
    if not names:
        return []
    payload = {"metabolite_name": "\n".join(names)}
    try:
        resp = requests.post(REFMET_API_URL, data=payload, timeout=40)
        resp.raise_for_status()
        text = resp.text.strip()
        if not text:
            logger.warning("Empty RefMet response")
            return []
        reader = csv.DictReader(StringIO(text), delimiter='\t')
        results = []
        for row in reader:
            cleaned = {k.strip(): (v if v != '-' else '') for k, v in row.items()}
            if any(cleaned.values()):
                results.append(cleaned)
        return results
    except requests.RequestException as e:
        logger.error(f"RefMet API error: {e}")
        return []
    except Exception as e:
        logger.error(f"Parse error: {e}")
        return []

def process_batches(lm_entries: List[Dict[str, str]]) -> List[Dict]:
    """Batch process LM IDs; submit all names (one per LM ID) to RefMet and join results back to LM IDs."""
    # Limit to first 20 for test run
    lm_entries = lm_entries
    all_results: List[Dict] = []
    total_batches = (len(lm_entries) + BATCH_SIZE - 1) // BATCH_SIZE
    logger.info(f"Processing {len(lm_entries)} LM IDs/names in {total_batches} batches (batch size {BATCH_SIZE})")

    for i in range(0, len(lm_entries), BATCH_SIZE):
        batch_num = i // BATCH_SIZE + 1
        batch_entries = lm_entries[i:i + BATCH_SIZE]
        batch_names = [entry["name"] for entry in batch_entries]
        logger.info(f"Batch {batch_num}/{total_batches} ({len(batch_names)} names)")
        refmet_rows = query_refmet(batch_names)

        # Map each LM ID to its name and join with RefMet results by name
        name_to_lm_ids = {}
        for entry in batch_entries:
            name_to_lm_ids.setdefault(entry["name"], []).append(entry["lm_id"])

        for row in refmet_rows:
            input_name = row.get("Input name", "")
            lm_ids = name_to_lm_ids.get(input_name, [])
            for lm_id in lm_ids:
                refmet_lm_id = row.get("LM_ID", "")
                lm_id_mismatch = bool(refmet_lm_id) and refmet_lm_id != lm_id
                combined = {
                    "lm_id": lm_id,
                    "lm_name": input_name,
                    "refmet_standardized_name": row.get("Standardized name", ""),
                    "refmet_formula": row.get("Formula", ""),
                    "refmet_exact_mass": row.get("Exact mass", ""),
                    "refmet_super_class": row.get("Super class", ""),
                    "refmet_main_class": row.get("Main class", ""),
                    "refmet_sub_class": row.get("Sub class", ""),
                    "refmet_pubchem_cid": row.get("PubChem_CID", ""),
                    "refmet_chebi_id": row.get("ChEBI_ID", ""),
                    "refmet_hmdb_id": row.get("HMDB_ID", ""),
                    "refmet_lm_id": refmet_lm_id,
                    "refmet_kegg_id": row.get("KEGG_ID", ""),
                    "refmet_inchi_key": row.get("INCHI_KEY", ""),
                    "refmet_id": row.get("RefMet_ID", ""),
                    "lm_id_mismatch": lm_id_mismatch,
                }
                all_results.append(combined)

        if batch_num < total_batches:
            time.sleep(DELAY_BETWEEN_BATCHES)

    logger.info(f"Completed; collected {len(all_results)} joined rows")
    return all_results

def write_output_csv(results: List[Dict], output_filename: str) -> None:
    if not results:
        logger.warning("No results to write")
        return
    fieldnames = [
        "lm_id",
        "lm_name",
        "refmet_standardized_name",
        "refmet_formula",
        "refmet_exact_mass",
        "refmet_super_class",
        "refmet_main_class",
        "refmet_sub_class",
        "refmet_pubchem_cid",
        "refmet_chebi_id",
        "refmet_hmdb_id",
        "refmet_lm_id",
        "refmet_kegg_id",
        "refmet_inchi_key",
        "refmet_id",
        "lm_id_mismatch",
    ]
    logger.info(f"Writing {len(results)} rows to {output_filename}")
    with open(output_filename, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(results)
    logger.info("Write complete")

def analyze_results(results: List[Dict]) -> None:
    total = len(results)
    logger.info("\n" + "="*60)
    logger.info("SUMMARY STATISTICS")
    logger.info("="*60)
    if total == 0:
        logger.info("No results to analyze.")
        logger.info("="*60 + "\n")
        return

    def pct(n: int) -> str:
        return f"{(n/total*100):.1f}%"

    matched_std = sum(1 for r in results if r.get("refmet_standardized_name"))
    matched_lm = sum(1 for r in results if r.get("refmet_lm_id"))
    mismatches = sum(1 for r in results if r.get("lm_id_mismatch"))

    logger.info(f"Total joined rows: {total}")
    logger.info(f"Rows with standardized name: {matched_std} ({pct(matched_std)})")
    logger.info(f"Rows with RefMet LM_ID: {matched_lm} ({pct(matched_lm)})")
    logger.info(f"LM ID mismatches (RefMet LM_ID != our LM ID): {mismatches} ({pct(mismatches)})")

    # Example problematic LM IDs (limit 10)
    if mismatches:
        logger.info("Examples of LM ID mismatches:")
        count = 0
        for r in results:
            if r.get("lm_id_mismatch"):
                logger.info(f"  {r['lm_id']} | name={r['lm_name']} | refmet={r['refmet_lm_id']}")
                count += 1
                if count >= 10:
                    break
    logger.info("="*60 + "\n")

def main():
    input_file = "lm_id_names.csv"
    output_file = "lm_refmet_joined.csv"
    try:
        lm_entries = load_lm_data(input_file)
        results = process_batches(lm_entries)
        write_output_csv(results, output_file)
        analyze_results(results)
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
        raise

if __name__ == "__main__":
    main()