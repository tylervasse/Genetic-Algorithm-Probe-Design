import matplotlib.pyplot as plt
import matplotlib.patches as patches

def annotate_complementary_pairs(dot_bracket):
    """
    Annotate the entire dot-bracket structure by mapping each base to its complementary base.
    Returns a dictionary where each base index is mapped to its complementary index.
    """
    n = len(dot_bracket)
    complement_map = {}
    stack = []

    for i, char in enumerate(dot_bracket):
        if char == '(':
            stack.append(i)
        elif char == ')':
            # Complementary base found; map the positions
            opening_idx = stack.pop()
            complement_map[opening_idx] = i
            complement_map[i] = opening_idx

    return complement_map

def is_balanced(dot_bracket):
    """Check if a dot-bracket string is balanced (equal number of '(' and ')')."""
    return dot_bracket.count('(') == dot_bracket.count(')')

def is_fully_closed(dot_bracket_segment):
    """Ensure no unpaired brackets left and no long-range pairing beyond this segment."""
    stack = []
    for c in dot_bracket_segment:
        if c == '(':
            stack.append(1)
        elif c == ')':
            if not stack:
                return False
            stack.pop()
    return not stack

def find_fully_enclosed_segments(sequence, dot_bracket, max_len=500, min_len=15):
    """
    Find the largest possible fully enclosed segments of length less than `max_len`.
    These segments must be fully closed (balanced and no open pairs).
    """
    n = len(dot_bracket)
    used = [False] * n
    segments = []

    # Annotate the dot-bracket structure for complementary base pairs
    complement_map = annotate_complementary_pairs(dot_bracket)

    i = 0
    while i < n:
        if used[i]:
            i += 1
            continue

        found = False
        # Try to grow as large as possible
        for j in range(i + min_len, min(i + max_len, n) + 1):
            if any(used[i:j]):
                break
            segment_db = dot_bracket[i:j]

            if is_balanced(segment_db) and is_fully_closed(segment_db):
                # Found a perfect region
                segments.append((i, j))
                for k in range(i, j):
                    used[k] = True
                i = j  # move past this segment
                found = True
                break  # move to next i

        if not found:
            i += 1  # Move forward if no segment found

    # Merge adjacent segments that are less than 500bp in total length and have no gap between them
    merged_segments = []
    i = 0
    while i < len(segments):
        start, end = segments[i]
        # Check if the next segment starts immediately after this one (i.e., no gap)
        while (i + 1 < len(segments) and 
               segments[i + 1][0] == end and  # No gap between segments
               (end - start + segments[i + 1][1] - segments[i + 1][0]) <= max_len):  # Check if merged length <= max_len
            # Merge the current segment with the next one
            end = segments[i + 1][1]
            i += 1
        
        merged_segments.append((start, end))
        i += 1

    return merged_segments

def plot_segment_distribution(segments, sequence_length):
    """
    Plot the distribution of fully enclosed segments in the RNA sequence.
    """
    fig, ax = plt.subplots(figsize=(15, 2))
    ax.set_xlim(0, sequence_length)
    ax.set_ylim(0, 1)
    for start, end in segments:
        rect = patches.Rectangle((start, 0.2), end - start, 0.6, linewidth=1, edgecolor='black', facecolor='skyblue')
        ax.add_patch(rect)
    ax.set_yticks([])
    ax.set_xlabel("Nucleotide Position")
    ax.set_title("RNA Segment Distribution")
    plt.tight_layout()
    plt.savefig("rna_segment_distribution.png")
    plt.show()
    plt.close()


def identify_gap_regions(segments, sequence_length):
    """
    Identify the gaps between the fully closed segments.
    Gaps are defined as regions not included in any segment.
    """
    gap_regions = []
    last_end = 0  # Starting point for the first gap

    for start, end in segments:
        if start > last_end:
            gap_regions.append((last_end, start-1))  # Add gap region
        last_end = end

    # If there's a remaining gap at the end of the sequence
    if last_end < sequence_length:
        gap_regions.append((last_end+1, sequence_length-1))

    return gap_regions


def find_complementary_gap_pairs(gap_portions, complement_map):
    """
    For each base in each gap, find complementary bases in other gaps based on the complement_map.
    """
    complementary_pairs = []

    # Iterate through each gap portion to compare with other gap portions
    for i, (gap_start, gap_end) in enumerate(gap_portions):
        for j, (other_gap_start, other_gap_end) in enumerate(gap_portions):
            if i != j:  # Don't compare a gap with itself
                # Check complementarity between bases in gap i and gap j
                for base_i in range(gap_start, gap_end+1):
                    for base_j in range(other_gap_start, other_gap_end+1):
                        # Check if base_i is complementary to base_j
                        if complement_map.get(base_i) == base_j:
                            complementary_pairs.append((base_i, base_j))

    return complementary_pairs