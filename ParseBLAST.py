import re
from collections import defaultdict
from dataclasses import dataclass
from ..pandas.panda_tricks import df

# Define a dataclass to hold alignment block information, similar to the C++ class.
END_OF_QUERY = ("Effective search",">") #,"S2:", "BLAST", "TBLASTN")

@dataclass
class AlignmentBlock:
    """Represents a single alignment block in a BLAST hit."""
    q_start: int = 0
    h_start: int = 0
    q_end: int = 0
    h_end: int = 0
    positives: int = 0
    gaps: int = 0
    e_val: float = 0.0

    @property
    def similarity(self) -> int:
        return (2 * self.positives - self.gaps) // 2

    @property
    def q_len(self) -> int:
        return self.q_end - self.q_start + 1

    @property
    def h_len(self) -> int:
        return self.h_end - self.h_start + 1


def parse_block(block_text: str) -> AlignmentBlock:
    """Parses a string containing a single BLAST alignment block."""
    bl = AlignmentBlock()

    # Use regex to find values more robustly.
    e_val_match = re.search(r"Expect = ([\de.-]+)", block_text)
    if e_val_match:
        bl.e_val = float(e_val_match.group(1))

    positives_match = re.search(r"Positives = (\d+)/\d+", block_text)
    if positives_match:
        bl.positives = int(positives_match.group(1))

    gaps_match = re.search(r"Gaps = (\d+)/\d+", block_text)
    bl.gaps = int(gaps_match.group(1)) if gaps_match else 0

    # Extract all query and sbjct coordinates
    coords = re.findall(r"(?:Query|Sbjct)\s+(\d+)\s+.*\s+(\d+)", block_text)
    if coords:
        q_coords = [c for i, c in enumerate(coords) if i % 2 == 0]
        h_coords = [c for i, c in enumerate(coords) if i % 2 != 0]

        bl.q_start = int(q_coords[0][0])
        bl.q_end = int(q_coords[-1][1])
        bl.h_start = int(h_coords[0][0])
        bl.h_end = int(h_coords[-1][1])

        # For TBLASTN, the hit can be reversed.
        if bl.h_end < bl.h_start:
            bl.h_start, bl.h_end = bl.h_end, bl.h_start
            
    return bl


def get_blocks(alignment_text: str) -> list[AlignmentBlock]:
    """Extracts all alignment blocks from a hit's alignment section."""
    blocks = alignment_text.split('\n\n\n')
    return [parse_block(block) for block in blocks if block]


def get_hit_descr(hit_alignment: str) -> tuple[str, int,str]:
    """Parses the description and length of a BLAST hit.
    output:
      Returns a tuple of (description, length, alignment).
    """
    hit_descr = ""
    hit_len = 0
    
    # Find the end of the description line(s)
    length_match = re.search(r"Length\s*=\s*([\d,]+)", hit_alignment)
    if not length_match:
        return "", 0

    header_end_pos = length_match.start()
    header = hit_alignment[:header_end_pos].strip()
    
    # Clean up the description string
    hit_descr = header.lstrip('>').replace('\n', ' ').strip()
    hit_descr = re.sub(r'\s+', ' ', hit_descr)
    
    hit_len = int(length_match.group(1).replace(",", ""))
    
    return hit_descr, hit_len, hit_alignment[length_match.end():]


def _calculate_merged_score(blocks: list[AlignmentBlock], coord_type: str) -> tuple[int, int]:
    """Helper function to calculate merged score for query or hit."""
    if not blocks:
        return 0, 0

    # Determine which attributes to use based on coord_type ('q' or 'h')
    start_attr = f"{coord_type}_start"
    end_attr = f"{coord_type}_end"
    len_attr = f"{coord_type}_len"
    
    # Sort blocks by their start coordinate
    blocks.sort(key=lambda b: (getattr(b, start_attr), getattr(b, end_attr)))

    merged_segments = []
    
    # Initial segment from the first block
    first_block = blocks[0]
    merged_segments.append([
        first_block.similarity,
        getattr(first_block, start_attr),
        getattr(first_block, len_attr)
    ])

    for i in range(1, len(blocks)):
        current_block = blocks[i]
        last_segm = merged_segments[-1]
        
        previous_similarity, previous_start, previous_len = last_segm
        previous_end = previous_start + previous_len - 1

        current_similarity = current_block.similarity
        current_start = getattr(current_block, start_attr)
        current_len = getattr(current_block, len_attr)
        current_end = getattr(current_block, end_attr)

        if current_start > previous_end: # No overlap
            merged_segments.append([current_similarity, current_start, current_len])
        else: # Overlap exists
            overlap_length = previous_end - current_start + 1
            if overlap_length >= current_len:
                continue # Current segment is fully contained and not better

            # Approximate similarity contributions
            overlap_sim1 = (previous_similarity * overlap_length / previous_len) if previous_len else 0
            overlap_sim2 = (current_similarity * overlap_length / current_len) if current_len else 0
            max_overlap_sim = max(overlap_sim1, overlap_sim2)

            prev_sim_before_ol = (previous_similarity * (previous_len - overlap_length) / previous_len) if previous_len else 0
            cur_sim_after_ol = (current_similarity * (current_len - overlap_length) / current_len) if current_len else 0
            
            average_sim = prev_sim_before_ol + max_overlap_sim + cur_sim_after_ol
            
            # Update the last segment with merged values
            last_segm[2] += current_end - previous_end # new length
            last_segm[0] = min(average_sim, last_segm[2]) # new similarity

    total_score = sum(seg[0] for seg in merged_segments)
    align_start = merged_segments[0][1]
    align_end = max(seg[1] + seg[2] for seg in merged_segments)
    total_align_length = align_end - align_start

    return int(total_score), total_align_length


def similarity_score(q_len: int, hit_alignment: str, is_tblastn: bool = False) -> dict:
    """Calculates a composite similarity score from a BLAST hit alignment."""
    results = {
        "hit_descr": "", "best_eval": 10.0, "align_score": 0.0,
        "q_align_len": 0, "q_score": 0, "h_align_len": 0,
        "hit_len": 0, "h_score": 0
    }
    
    hit_descr, hit_len,alignments = get_hit_descr(hit_alignment)
    if not hit_len: return results
    
    results["hit_descr"] = hit_descr
    results["hit_len"] = hit_len
    
    blocks = get_blocks(alignments)
    if not blocks: return results

    results["best_eval"] = blocks[0].e_val

    # TBLASTN requires hit coordinates to be adjusted for codon space
    if is_tblastn:
        for block in blocks:
            block.h_start = (block.h_start // 3)
            block.h_end = (block.h_end // 3)

    q_score, q_align_len = _calculate_merged_score(list(blocks), 'q')
    h_score, h_align_len = _calculate_merged_score(list(blocks), 'h')
    
    results["q_score"] = q_score
    results["q_align_len"] = q_align_len
    results["h_score"] = h_score
    results["h_align_len"] = h_align_len

    # Final score calculation
    score1 = (h_score / hit_len) if hit_len else 0
    score2 = (q_score / q_len) if q_len else 0
    
    if is_tblastn:
        score1 *= 3

    results["align_score"] = (score1 + score2) / 2
    return results

import re
from typing import Iterator, Dict

# Assume the 'similarity_score' function is defined as before.

def generate_scores_fromGermini(blast_results_path: str, blast_type: int) -> Iterator[Dict]:
    """
    Parses a BLAST output file and yields a score dictionary for each hit as it's found.
    This optimized version processes the file as a stream, making it highly memory-efficient.
    """
    is_tblastn = (blast_type != 0)
    q_len = 0
    q_descr = ""
    hit_lines = []

    with open(blast_results_path, 'r') as f:
      for line in f:
        if line.startswith((END_OF_QUERY)):
          if hit_lines:
            scores = similarity_score(q_len, hit_alignment="".join(hit_lines), is_tblastn=is_tblastn)
            scores['q_descr'] = q_descr
            scores['q_len'] = q_len
            yield scores
            hit_lines = []  # Reset for the next block

          if line.startswith("Query="):
            # A new query section begins.
            q_descr_parts = [line.split("=", 1)[1].strip()]           
            while line: # Read ahead to find the query length
              line = str(next(f,default=""))
              len_match = re.search(r"(\d{1,3}(?:,\d{3})*|\d+)\s+letters|Length=(\d+)", line)
              if len_match:
                num_str = next(s for s in len_match.groups() if s is not None)
                q_len = int(num_str.replace(",", ""))
                break  # Length found
              stripped_line = line.strip()
              if stripped_line:
                q_descr_parts.append(stripped_line)
            q_descr = " ".join(q_descr_parts)

          elif line.startswith('>') or hit_lines: # Start of a new hit. Begin collecting its lines or keep collecting
            hit_lines.append(line) 
            

def generate_scores(blast_results_path: str,blast_type: int):
    """Parses a BLAST output file and generates a tab-separated map of hits."""
    is_tblastn = (blast_type != 0)
    def _qlen(l: str) -> bool:
      len_match = re.search(r"(\d{1,3}(?:,\d{3})*|\d+)\s+letters", l) or re.search(r"Length=(\d+)", l)     
      return int(len_match.group(1).replace(",", "")) if len_match else 0

    with open(blast_results_path, 'r') as f_in:
      query_alignments = []
      for line in f_in:
        line = line.strip()
        query_alignments.append(line)
        if line.startswith(END_OF_QUERY):
          i = 0
          while not query_alignments[i].startswith("Query="):
            i += 1
          q_descr = [query_alignments[i].split("=", 1)[1].strip()]

          next_line = i+1
          for q_line in query_alignments[next_line:]:
            i += 1
            q_len = _qlen(q_line)
            if q_len:
              q_descr = " ".join(q_descr)
              break
            else:
              if q_line:
                q_descr.append(q_line)

          next_line = i+1
          for q_line in query_alignments[next_line:]:
            i += 1
            if query_alignments[i].startswith(">"):
              break

          hit_alignment = query_alignments[i]
          for q_line in query_alignments[i+1:]:
            if q_line.startswith(">"):
              query_scores = similarity_score(q_len, hit_alignment,is_tblastn)
              query_scores['q_descr'] = q_descr
              query_scores['q_len'] = q_len
              yield query_scores
              hit_alignment = q_line
            else:
              hit_alignment += q_line + '\n'

          query_alignments = []


def generate_blast_map(blast_results_path: str, output_path: str, alscore_cutoff: float, blast_type: int):
    """Parses a BLAST output file and generates a tab-separated map of hits."""
    with (open(output_path, 'w')) as f_out:
      for scores in generate_scores_fromGermini(blast_results_path, blast_type):
        q_descr = scores['q_descr']
        h_descr = scores['hit_descr']
        if scores["align_score"] > alscore_cutoff and q_descr != scores["hit_descr"]:
          f_out.write(
            f"{q_descr}\t{scores['hit_descr']}\t{scores['align_score']:.3f}\t"
            f"{scores['best_eval']:.1E}\t{scores['q_align_len']}\t{scores['q_len']}\t"
            f"{scores['q_score']}\t{scores['h_align_len']}\t{scores['hit_len']}\t"
            f"{scores['h_score']}\n"
          )

                

def get_best_pair(pairs: list[str]) -> tuple[str, float, float, int]:
    """Finds the best hit from a list of formatted pair strings."""
    best_score = -1.0
    best_eval = 1e10
    best_pair_descr = ""
    best_hit_index = -1

    for i, pair_str in enumerate(pairs):
        parts = pair_str.split('\t')
        hit_descr, score, ev = parts[0], float(parts[1]), float(parts[2])

        if score > best_score:
            best_score = score
            best_eval = ev
            best_pair_descr = hit_descr
            best_hit_index = i
        elif score == best_score and ev < best_eval:
            best_eval = ev
            best_pair_descr = hit_descr
            best_hit_index = i
            
    return best_pair_descr, best_score, best_eval, best_hit_index


def _is_unique_pair(q_descr: str, best_query_of_bh: str, q_score: float, bh_score: float, q_eval: float, bh_eval: float) -> bool:
    """
    Determines if a query-hit pair is unique based on reciprocality and scores,
    mirroring the logic from the C++ IsUniquePair function.
    
    A pair is considered unique if:
    1. It's a true reciprocal best hit.
    2. The current query has a definitively better alignment score for the hit.
    3. Scores are tied, but the current query has a better (lower) E-value.
    """
    if q_descr == best_query_of_bh:
        return True
    if q_score > bh_score:
        return True
    if q_score == bh_score and q_eval < bh_eval:
        return True
    return False


def _find_best_pairs_pass(q_to_h: dict[str, set[str]], h_to_q: dict[str, list[str]], best_unique_pairs: dict[str, str]):
    """
    Executes a single pass of the reciprocal best hit algorithm.
    This function depletes `h_to_q` from found hits and adds them to `best_unique_pairs`.
    """
    # Sorting keys ensures deterministic behavior when scores are tied.
    for q_descr in sorted(q_to_h.keys()):
        # Skip queries that have already been paired in a previous pass.
        if q_descr in best_unique_pairs:
            continue

        candidate_pairs = list(q_to_h.get(q_descr, []))
        # Continuously check the best available hit for the current query.
        while candidate_pairs:
            best_hit_descr, q_score, q_eval, bh_index = get_best_pair(candidate_pairs)
            
            best_hit_pairs = h_to_q.get(best_hit_descr)
            
            # If the best hit is no longer available (it was claimed by another query),
            # remove it from this query's list and try again.
            if not best_hit_pairs:
                candidate_pairs = [p for i, p in enumerate(candidate_pairs) if i != bh_index]
                continue

            best_query_of_bh, bh_score, bh_eval, _ = get_best_pair(best_hit_pairs)

            if _is_unique_pair(q_descr, best_query_of_bh, q_score, bh_score, q_eval, bh_eval):
                # Success! A unique pair was found.
                best_unique_pairs[q_descr] = candidate_pairs[bh_index]
                # Claim this hit by removing it from the pool of available hits.
                del h_to_q[best_hit_descr]
                # Stop processing this query and move to the next one.
                break 
            else:
                # This hit is not a unique partner. Remove it from consideration for this
                # query and check the next-best hit in the next loop iteration.
                candidate_pairs = [p for i, p in enumerate(candidate_pairs) if i != bh_index]


def get_unique_hits(q_to_h: dict[str, set[str]], h_to_q: dict[str, set[str]]) -> dict[str, str]:
    """
    Identifies reciprocal best hits from query-to-hit and hit-to-query maps.

    This function uses a two-pass approach to resolve conflicts:
    1. First pass finds all the "easy" reciprocal best hits.
    2. Second pass runs on the remaining queries to find pairs among them,
       now that the initial high-confidence pairs have been removed.
    """
    best_unique_pairs = {}
    original_bup_size = 0
    
    # --- Pass 1 ---
    # Create copies of the dictionaries to avoid modifying the originals directly in the first pass.
    h_to_q_pass1 = {k: list(v) for k, v in h_to_q.items()}
    _find_best_pairs_pass(q_to_h, h_to_q_pass1, best_unique_pairs)

    # If the first pass didn't assign a pair to every query, run a multiple passes until .
    # The pool of available hits (h_to_q_pass1) has already been reduced by the first pass.
    while len(best_unique_pairs) - original_bup_size > 0:
      original_bup_size = len(best_unique_pairs)
      _find_best_pairs_pass(q_to_h, h_to_q_pass1, best_unique_pairs)
        
    return best_unique_pairs


def blast_reciprocal_map(blast_results_path:str, out_path:str,align_score_cutoff:float=0.3, blast_type:int=0) -> dict[str, str]:
    """
    Parses a BLAST output file to find reciprocal best hits.

    This function reads through a BLAST report, calculates a custom similarity
    score for each query-hit pair, and then applies a two-pass algorithm to
    identify unique, reciprocal best-hit pairs that satisfy the score cutoff.

    Args:
        blast_results_path: Path to the BLAST output file.
        align_score_cutoff: The minimum similarity score for a pair to be considered.
        blast_type: An integer indicating the BLAST type (0 for BLASTP, non-zero for TBLASTN).

    Returns:
        A dictionary mapping each query description to its best unique hit's
        description and associated scores.
    """
    q_to_h = defaultdict(set)
    h_to_q = defaultdict(set)
    
    for scores in generate_scores_fromGermini(blast_results_path, blast_type):
      q_descr = scores['q_descr']
      h_descr = scores['hit_descr']
      # Apply the score cutoff and avoiding self hits in case of paralog calculations.
      if scores["align_score"] < align_score_cutoff or q_descr == scores["hit_descr"]:
        continue

      # Format the output string for storage.
      add_descr = [f'{scores['align_score']:.3f}',f'{scores['best_eval']:.3e}']
      add_descr += list(map(str, [scores['q_align_len'],scores['q_len'],scores['q_score'],
                    scores['h_align_len'],scores['hit_len'],scores['h_score']])
                    )

      # Store the forward relationship (q_to_h:Query->Hit, h_to_q:Hit->Query).
      q_to_h[q_descr].add('\t'.join([h_descr]+add_descr))
      h_to_q[h_descr].add('\t'.join([q_descr]+add_descr))
      
    # --- 3. Resolve the Maps to Find Unique Best Hits ---
    best_unique_pairs = get_unique_hits(q_to_h, h_to_q)
    rows = []
    for q_descr, pair_str in best_unique_pairs.items():
        info = pair_str.split('\t')
        rows.append([q_descr]+info)
      
    df_header = ["Query","Hit","Average alignment Score","Best Eval",
                 "Query alignment length","Query Length","Query alignment Score",
                 "Hit alignment length","Hit Length","Hit alignment Score"]
    BRHmap = df.from_rows(rows, df_header)

    for c in ['Average alignment Score','Best Eval']: 
       BRHmap[c] = BRHmap[c].astype(float)
    for c in ['Query alignment Score','Hit alignment Score','Query alignment length','Query Length','Hit alignment length','Hit Length']: 
       BRHmap[c] = BRHmap[c].astype(int)

    BRHmap = BRHmap.sortrows(by=["Average alignment Score","Best Eval"], ascending=[False, True])
    BRHmap.to_csv(out_path,sep='\t', index=False)

    return best_unique_pairs