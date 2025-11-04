# Project 3 - Multiple Sequence Alignment & Database Search

A comprehensive implementation of **Star Alignment** algorithm with **block-based iterative refinement** for multiple sequence alignment (MSA), combined with theoretical exploration of database search algorithms (FASTA, BLAST) and progressive alignment methods (ClustalW, T-Coffee).

[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![Algorithm](https://img.shields.io/badge/Algorithm-Star%20Alignment-green.svg)](#)
[![Course](https://img.shields.io/badge/Course-Bioinformatics-blue.svg)](#)
[![University](https://img.shields.io/badge/University-AUT-red.svg)](https://aut.ac.ir/en)

<details> <summary><h2>📚 Table of Contents</h2></summary>

- [Overview](#-overview)
- [Multiple Sequence Alignment Fundamentals](#-multiple-sequence-alignment-fundamentals)
- [Programming Assignment](#-programming-assignment)
  - [Star Alignment Algorithm](#star-alignment-algorithm)
  - [Block-Based Improvement](#block-based-improvement)
  - [Implementation Details](#implementation-details)
  - [Input/Output Format](#inputoutput-format)
  - [Usage Examples](#usage-examples)
- [Theoretical Assignment](#-theoretical-assignment)
  - [Question 1: FASTA vs BLAST vs Dynamic Programming](#question-1-fasta-vs-blast-vs-dynamic-programming)
  - [Question 2: BLAST Word Search Tree](#question-2-blast-word-search-tree)
  - [Question 3: FASTA Longest Common Substring](#question-3-fasta-longest-common-substring)
  - [Question 4: Star Alignment Manual Calculation](#question-4-star-alignment-manual-calculation)
  - [Question 5: Star vs ClustalW vs T-Coffee](#question-5-star-vs-clustalw-vs-t-coffee)
- [Installation & Setup](#%EF%B8%8F-installation--setup)
- [Algorithm Complexity](#%EF%B8%8F-algorithm-complexity)
- [Results Summary](#-results-summary)
- [Project Structure](#%EF%B8%8F-project-structure)
- [Key Concepts Covered](#-key-concepts-covered)
- [Learning Outcomes](#-learning-outcomes)
- [Project Information](#ℹ%EF%B8%8F-project-information)
- [Contact](#-contact)

</details>

---

## 📋 Overview

This project explores **Multiple Sequence Alignment (MSA)**, a cornerstone technique in bioinformatics for comparing three or more biological sequences simultaneously. MSA is essential for:

- **Phylogenetic Analysis**: Constructing evolutionary trees
- **Conserved Region Detection**: Identifying functional domains
- **Protein Family Classification**: Grouping related proteins
- **Structure Prediction**: Leveraging homology modeling
- **Database Searching**: Finding similar sequences efficiently

The project consists of two main components:

### 1. Programming Component
**Implementation**: Star Alignment with iterative block-based refinement
- Initial MSA using center-star heuristic
- Iterative improvement through block realignment
- Scoring with match=+3, mismatch=-1, gap=-2
- Outputs alignment score and sequences in input order

### 2. Theoretical Component
**Topics Covered**:
- Database search algorithms (FASTA, BLAST, Dynamic Programming)
- BLAST word indexing and search tree construction
- FASTA k-tuple matching and longest common substring
- Manual star alignment calculation with DNAfull matrix
- Comparative analysis: Star vs ClustalW vs T-Coffee

---

## 🧬 Multiple Sequence Alignment Fundamentals

### Why MSA?

**Pairwise vs Multiple Alignment:**
```
Pairwise (2 sequences):
  Seq1: ACTG
  Seq2: AC-G

Multiple (3+ sequences):
  Seq1: ACTG--
  Seq2: AC-G--
  Seq3: A--GTC
  
Advantage: Reveals patterns invisible in pairwise comparisons
```

**Key Applications:**

1. **Evolutionary Studies**
   - Trace evolutionary relationships
   - Identify orthologs and paralogs
   - Construct phylogenetic trees

2. **Functional Analysis**
   - Detect conserved functional regions
   - Predict active sites
   - Identify structural motifs

3. **Comparative Genomics**
   - Compare genomes across species
   - Find regulatory elements
   - Understand gene families

### MSA Challenges

**Computational Complexity:**
```
Pairwise:     O(n²) - tractable with dynamic programming
3 sequences:  O(n³) - still manageable
k sequences:  O(nᵏ) - exponential explosion!

For 10 sequences of length 100:
Exact solution requires 10¹⁰⁰ operations (impossible!)
```

**Solution**: Use **heuristic algorithms** (Star, ClustalW, T-Coffee)
- Trade optimality for speed
- Provide good (not necessarily optimal) alignments
- Run in reasonable time (polynomial complexity)

### MSA Scoring

**Sum-of-Pairs (SP) Score:**
```
Total Score = Σ (all pairwise column scores)

Example column: [A, A, C]
Pairs: (A,A)=+3, (A,C)=-1, (A,C)=-1
Column score: 3 + (-1) + (-1) = 1
```

**Implementation Note**: This project uses SP scoring with:
- Match: +3
- Mismatch: -1
- Gap: -2
- (Gap, Gap): 0 (no penalty for double gaps)

---

## 💻 Programming Assignment

### Star Alignment Algorithm

**Overview**: Center-star heuristic for fast approximate MSA

**Core Idea**:
1. Find a "center" sequence most similar to all others
2. Align all sequences to this center
3. Merge alignments using "once a gap, always a gap" rule

**Why "Star"?**
```
Visual representation:

    Seq2
     |
Seq3-Center-Seq4
     |
    Seq1
    
Center is the "hub", others are "spokes"
```

**Algorithm Steps:**

#### Step 1: Build Pairwise Score Matrix

For each pair of sequences (i, j):
- Perform global alignment using provided `global_align()` function
- Store score in matrix S[i][j]

```python
S = [[0 for _ in range(n)] for _ in range(n)]

for i in range(n):
    for j in range(n):
        if i != j:
            _, _, score = global_align(seq[i], seq[j], 
                                       s_match=3, s_mismatch=-1, s_gap=-2)
            S[i][j] = score
```

**Score Matrix Example:**
```
      Seq1  Seq2  Seq3  Seq4
Seq1   -     20    15    25    → Total: 60
Seq2   20    -     18    22    → Total: 60
Seq3   15    18    -     30    → Total: 63  ← Center!
Seq4   25    22    30    -     → Total: 77
```

#### Step 2: Select Center Sequence

**Center Selection Rule:**
- Choose sequence with **maximum sum of pairwise scores**
- This sequence is most similar to all others on average

```python
center_seq = max(sequences, 
                 key=lambda s: sum(S[sequences.index(s)].values()))
```

**Why this works:**
- Center sequence has highest average similarity
- Minimizes total alignment errors
- Provides good "representative" for the family

#### Step 3: Progressive Alignment to Center

**Order**: Align sequences in **descending order** of their score with center
```
Best match first → Worst match last
```

**Process for each sequence**:
1. Align new sequence with current center
2. Update center alignment (add gaps)
3. **Propagate gaps** to all previously aligned sequences
4. Apply "once a gap, always a gap" rule

**"Once a Gap, Always a Gap" Rule:**
```
If center gets a gap at position i:
→ ALL previously aligned sequences must get gap at position i

Example:
Previous state:
  Center: ACTG
  Seq1:   AC-G

Aligning Seq2 introduces gap:
  Center: A-CTG
  Seq2:   AGCTG

Propagate gap to Seq1:
  Center: A-CTG
  Seq1:   A-C-G  ← Gap added at position 1
  Seq2:   AGCTG
```

**Implementation Function**: `always_a_gap()`
```python
def always_a_gap(center, updated_center, last_seqs):
    """
    Synchronize gaps between old and new center alignments
    Propagate new gaps to all previously aligned sequences
    """
    loc1 = 0  # Position in old center
    loc2 = 0  # Position in updated center
    loc3 = 0  # Position in final alignment
    
    while loc1 < len(center) and loc2 < len(updated_center):
        if center[loc1] == updated_center[loc2]:
            # Characters match, move forward
            loc1 += 1
            loc2 += 1
            loc3 += 1
        elif updated_center[loc2] == '-':
            # New gap in updated center
            # Insert gap in all previous sequences at position loc3
            for seq in last_seqs:
                last_seqs[...] = seq[:loc3] + '-' + seq[loc3:]
            loc2 += 1
            loc3 += 1
```

#### Step 4: Output Alignment

**Important**: Output sequences in **original input order**, not alignment order!

```python
# Find original order
output_seqs = []
for original_seq in input_sequences:
    for aligned_seq in aligned_seqs:
        if aligned_seq.replace('-', '') == original_seq:
            output_seqs.append(aligned_seq)
```

---

### Block-Based Improvement

**Motivation**: Star alignment is greedy and may be suboptimal
- Early alignment errors propagate
- Local regions might benefit from realignment

**Solution**: Iterative block refinement

#### Block Identification Rules

**Valid Block Criteria:**
1. ✅ **Minimum length**: ≥3 columns (at least 3 positions)
2. ✅ **Contains variation**: Has mismatches or gaps (not all identical)
3. ❌ **No all-gap columns**: Cannot have columns with only gaps

**Boundary Definition:**
- Block boundaries are all-gap or all-identical columns
- These columns are excluded from blocks

**Example:**
```
Alignment:
ACGT---TGCA
A-GT---TGCA
ACGT---TGCA
↑↑↑↑   ↑↑↑↑
│││└─→ All T (boundary, exclude)
││└──→ All G (boundary, exclude)
│└───→ Has variation (C/-/C) - VALID BLOCK
└────→ All A (boundary, exclude)

Blocks: [1:3] (positions with variation)
```

#### Improvement Process

**For each identified block:**

1. **Extract block sequences**
```python
block_seqs = []
for seq in aligned_seqs:
    # Remove gaps to get original subsequence
    block_seqs.append(seq[start:end].replace('-', ''))
```

2. **Realign block using Star algorithm**
```python
# Full MSA on block subsequences
new_block_alignment = star_alignment(block_seqs)
```

3. **Calculate scores**
```python
old_score = calculate_alignment_score(old_block_columns)
new_score = calculate_alignment_score(new_block_columns)
```

4. **Replace if improved**
```python
if new_score > old_score:
    # Substitute new alignment for old block
    for i, seq in enumerate(aligned_seqs):
        aligned_seqs[i] = (
            seq[:start] + 
            new_block_alignment[i] + 
            seq[end:]
        )
```

#### Iterative Refinement Loop

```python
improved = True
while improved:
    blocks = identify_blocks(alignment)
    improved = False
    
    for block in blocks:
        if realign_block(block) improves score:
            improved = True
            
    if not improved:
        break  # Convergence reached
```

**Convergence**: Algorithm stops when no block can be improved

---

### Implementation Details

#### Core Functions

**`global_align(x, y, s_match, s_mismatch, s_gap)`**
- Provided global alignment function (from Project 2)
- Returns: (aligned_x, aligned_y, score)
- Used for all pairwise alignments

**`fill_matrix_and_find_center(seqs)`**
```python
def fill_matrix_and_find_center(seqs):
    """
    Build pairwise score matrix
    Identify center sequence
    """
    max_score = float('-inf')
    center_seq = ''
    score_matrix = {}
    
    for seq1 in seqs:
        total_score = 0
        for seq2 in seqs:
            if seq1 != seq2:
                _, _, score = global_align(seq1, seq2, 3, -1, -2)
                score_matrix[...][...] = score
                total_score += score
        
        if total_score > max_score:
            max_score = total_score
            center_seq = seq1
    
    return score_matrix, center_seq
```

**`star_alignment(center_seq, score_matrix, seqs)`**
```python
def star_alignment(center_seq, score_matrix, seqs):
    """
    Progressive alignment to center sequence
    """
    # Sort sequences by score with center (descending)
    sorted_seqs = sorted(score_matrix[center_index].items(), 
                        key=lambda x: x[1], reverse=True)
    
    aligned_seqs = []
    new_center = center_seq
    
    for seq_index, score in sorted_seqs:
        last_center = new_center
        
        # Align new sequence with current center
        new_seq, new_center, _ = global_align(
            seqs[seq_index], new_center, 3, -1, -2
        )
        
        # Propagate gaps to previous sequences
        always_a_gap(last_center, new_center, aligned_seqs)
        
        aligned_seqs.append(new_seq)
    
    # Add center at beginning
    aligned_seqs = [new_center] + aligned_seqs
    
    return aligned_seqs
```

**`calculate_alignment_score(columns)`**
```python
def calculate_alignment_score(columns):
    """
    Sum-of-pairs scoring for MSA
    """
    total_score = 0
    
    for column in columns:
        for i in range(len(column)):
            for j in range(i + 1, len(column)):
                if column[i] == '-' and column[j] == '-':
                    total_score += 0  # No penalty for (gap, gap)
                elif column[i] == '-' or column[j] == '-':
                    total_score += -2  # Gap penalty
                elif column[i] != column[j]:
                    total_score += -1  # Mismatch
                elif column[i] == column[j]:
                    total_score += 3   # Match
    
    return total_score
```

**`block_improvement(seqs_columns, aligned_seqs)`**
```python
def block_improvement(seqs_columns, aligned_seqs):
    """
    Identify and realign suboptimal blocks
    """
    # Find all-gap or all-identical columns (boundaries)
    not_valid_column_indexes = []
    for i, column in enumerate(seqs_columns):
        if all(column[0] == char for char in column):
            not_valid_column_indexes.append(i)
    
    # Extract blocks between boundaries (min length 3)
    blocks = extract_blocks(not_valid_column_indexes, len(seqs_columns))
    
    improved = False
    for block in blocks:
        # Realign block and check improvement
        if realign_improves(block, aligned_seqs):
            improved = True
    
    return improved, aligned_seqs
```

---

### Input/Output Format

#### Input Format
```
Line 1: n (number of sequences)
Line 2 to n+1: Protein sequences (uppercase letters)
```

**Example Input:**
```
4
HEAGAWGHEA
TCIVMREAYE
YIMQEVQQER
WRYIAMREQYES
```

**Constraints:**
- Time limit: 5 seconds
- Memory limit: 512 MB
- n ≥ 2 sequences
- Protein sequences (20 standard amino acids)

#### Output Format
```
Line 1: Initial alignment score (before improvement)
Lines 2 to n+1: Initial aligned sequences (in INPUT order)
Line n+2: Empty line
Line n+3: Final alignment score (after improvement)
Lines n+4 to 2n+3: Final aligned sequences (in INPUT order)
```

**Example Output:**
```
51
-TYI-MREAQYESAQ
-TCIVMREA-YE---
--YI-MQEVQQER--
WRYIAMRE-QYES--

60
-TYI-MREAQYESAQ
-TCIVMREA-YE---
--YIMQEV-QQER--
WRYIAMRE-QYES--
```

**Output Requirements:**
- ✅ Score on first line
- ✅ Sequences in **original input order**
- ✅ Use `-` for gaps
- ✅ Maintain alignment columns
- ✅ Empty line between initial and final
- ❌ Do NOT sort sequences alphabetically

---

### Usage Examples

#### Running the Program

**Method 1: Interactive Input**
```bash
cd "3- Multiple Sequence Alignment - DB Search/src"
python main.py
# Then enter:
4
HEAGAWGHEA
TCIVMREAYE
YIMQEVQQER
WRYIAMREQYES
```

**Method 2: Input File**
```bash
# Create input.txt:
echo -e "4\nHEAGAWGHEA\nTCIVMREAYE\nYIMQEVQQER\nWRYIAMREQYES" > input.txt

# Run:
python main.py < input.txt
```

#### Sample Test Cases

**Test 1: Simple 4-sequence alignment**

**Input:**
```
4
TYIMREAQYESAQ
TCIVMREAYE
YIMQEVQQER
WRYIAMREQYES
```

**Expected Output:**
```
51
-TYI-MREAQYESAQ
-TCIVMREA-YE---
--YI-MQEVQQER--
WRYIAMRE-QYES--

51
-TYI-MREAQYESAQ
-TCIVMREA-YE---
--YI-MQEVQQER--
WRYIAMRE-QYES--
```

**Analysis:**
- Initial star alignment already optimal
- No improvement from block refinement
- Score: 51 (conserved regions well-aligned)

---

**Test 2: 5-sequence DNA alignment**

**Input:**
```
5
TAGCTACCAGGA
CAGCTACCAGG
TAGCTACCAGT
CAGCTATCGCGGC
CAGCTACCAGGA
```

**Expected Output:**
```
240
TAGCTA-C-CAGGA
CAGCTA-C-CAGG-
TAGCTA-C-CA-GT
CAGCTATCGC-GGC
CAGCTA-C-CAGGA

240
TAGCTA-C-CAGGA
CAGCTA-C-CAGG-
TAGCTA-C-CA-GT
CAGCTATCGC-GGC
CAGCTA-C-CAGGA
```

**Analysis:**
- High conservation in "CAGCTA" region
- Gap introduced to align "TCGC" insertion in seq4
- Alignment converged immediately

---

## 📝 Theoretical Assignment

### Question 1: FASTA vs BLAST vs Dynamic Programming

**Task**: Compare three sequence database search methods

#### Dynamic Programming (Smith-Waterman)

**Method**: Exhaustive search using optimal local alignment

**Algorithm**:
```python
# For each database sequence:
for db_seq in database:
    score[query][db_seq] = smith_waterman(query, db_seq)
    
# Return all scores above threshold
```

**Characteristics:**
- ✅ **Guaranteed optimal**: Finds best local alignment
- ✅ **Sensitive**: Detects distant homologs
- ✅ **Precise**: No false positives from algorithm
- ❌ **Slow**: O(mn) per comparison
- ❌ **Impractical**: For large databases (millions of sequences)

**Best For:**
- Small databases
- High-quality alignments needed
- Regulatory/critical applications
- Benchmark comparisons

---

#### FASTA Algorithm

**Method**: Heuristic search using k-tuple matching

**Algorithm Steps:**

**Step 1: K-tuple Hashing**
```
Query: CACGTTGACAT (k=2)
Tuples: CA, AC, CG, GT, TT, TG, GA, AC, CA, AT

Hash table:
AC → [2, 8]
CA → [1, 9]
CG → [3]
GA → [7]
GT → [4]
TG → [6]
TT → [5]
AT → [10]
```

**Step 2: Exact Match Search**
```
For each database sequence:
  Find all exact k-tuple matches with query
  Calculate diagonal offsets
```

**Step 3: Identify Regions**
```
Cluster nearby matches on same diagonal
Example:
  Match at (query_pos=2, db_pos=5) → diagonal = 3
  Match at (query_pos=8, db_pos=11) → diagonal = 3
  → Same diagonal = potential alignment region
```

**Step 4: Extend & Score**
```
Extend high-density diagonal regions
Calculate ungapped alignment score
Refine with gapped alignment (Smith-Waterman on region)
```

**Characteristics:**
- ✅ **More sensitive than BLAST**: Smaller k-tuples (k=1 or 2)
- ✅ **Finds distant homologs**: Lower threshold for matches
- ✅ **Single best result**: Returns top hit per database sequence
- ❌ **Slower than BLAST**: More exact matches to process
- ⚖️ **Medium speed**: Faster than DP, slower than BLAST

**Best For:**
- Short query sequences
- Distant homology searches
- Comprehensive sensitivity needed
- Moderate-sized databases

---

#### BLAST Algorithm

**Method**: Heuristic search using high-scoring word pairs

**Algorithm Steps:**

**Step 1: Word Generation (w=3 for proteins)**
```
Query: ACDEFG
Words: ACD, CDE, DEF, EFG

For each word, find "neighborhood":
  ACD → ACD (score=X), ACE (score>T), ADD (score>T), ...
  Keep only words scoring above threshold T with query word
```

**Step 2: Exact Word Search**
```
Search database for exact matches to:
  - Query words
  - High-scoring neighborhood words
  
Use indexed lookup (very fast)
```

**Step 3: Two-Hit Extension (BLAST2)**
```
Require TWO close hits on same diagonal:
  Hit1: Query[10] matches DB[50]
  Hit2: Query[15] matches DB[55]
  → Distance = 5, within window (e.g., 40) → Extend

Why: Reduces false positive extensions
```

**Step 4: Gapped Alignment**
```
Extend ungapped alignment
If score exceeds threshold:
  Perform gapped extension (dynamic programming)
Calculate E-value (statistical significance)
```

**Characteristics:**
- ✅ **Very fast**: Pre-indexed word lookup
- ✅ **Scalable**: Handles massive databases (nr, nt)
- ✅ **Multiple results**: Returns ranked list with E-values
- ✅ **Family variants**: BLASTn, BLASTp, BLASTx, tBLASTn
- ❌ **Less sensitive**: May miss distant homologs
- ⚖️ **Good specificity**: Two-hit requirement reduces noise

**Best For:**
- Large-scale database searches
- Recent homologs (moderate similarity)
- High-throughput applications
- NCBI nr/nt databases

---

#### Comparison Summary

| Feature | Dynamic Programming | FASTA | BLAST |
|---------|-------------------|-------|-------|
| **Method** | Exhaustive search | K-tuple hashing | Word neighborhood |
| **Speed** | Slowest (O(mn) per seq) | Medium | Fastest |
| **Sensitivity** | Highest (optimal) | High | Moderate |
| **Specificity** | Highest | High | Good |
| **K-tuple size** | N/A | k=1-2 (small) | w=3-11 (large) |
| **Initial match** | N/A | Exact k-tuple | Exact word + neighbors |
| **Extension** | Full DP | DP on regions | Two-hit, then gapped |
| **Database size** | Small (<1000) | Medium | Large (millions) |
| **Results** | All above threshold | Top hit per sequence | Ranked list with E-values |
| **Use case** | Gold standard | Sensitivity > speed | Speed > sensitivity |

**Memory Usage:**
- **DP**: O(mn) per comparison (can optimize to O(m))
- **FASTA**: O(m + database_size) for hash table
- **BLAST**: O(database_index) - largest upfront, but reused

**Accuracy vs Speed Trade-off:**
```
DP ────────→ FASTA ────────→ BLAST
↑                                  ↑
Slow, optimal           Fast, good approximation
```

**From Report Summary:**
- DP guarantees best answer but impractical for large searches
- FASTA: Better sensitivity (smaller k), slower, good for short queries
- BLAST: Speed optimized (indexed, two-hit), family variants, industry standard

---

### Question 2: BLAST Word Search Tree

**Task**: Find all 3-letter words aligning with "QIV" with score ≥ 10

#### Part A: Word Identification

**Given Scoring Matrix** (example values):
```
    Q   I   V   L   M   E
Q   5   -   -   -   -   2
I   -   4   3   2   1   -
V   -   3   4   1   1   -
L   -   2   1   -   -   -
M   -   1   1   -   -   -
E   2   -   -   -   -   -
```

**Calculate scores for all 3-letter combinations:**

**QIV** (query word):
```
QIV = Q-Q + I-I + V-V = 5 + 4 + 4 = 13 ✓
```

**Q** in position 1:
```
QII = 5 + 4 + 3 = 12 ✓
QIL = 5 + 4 + 1 = 10 ✓
QIM = 5 + 4 + 1 = 10 ✓
QVV = 5 + 3 + 4 = 12 ✓
QVI = 5 + 3 + 3 = 11 ✓
QLV = 5 + 2 + 4 = 11 ✓
QLI = 5 + 2 + 3 = 10 ✓
QMV = 5 + 1 + 4 = 10 ✓
```

**E** in position 1:
```
EIV = 2 + 4 + 4 = 10 ✓
```

**Final Word List** (10 words with score ≥ 10):
```
QIV (13), QII (12), QVV (12), QVI (11), QLV (11), 
QIL (10), QIM (10), QLI (10), QMV (10), EIV (10)
```

#### Part B: Word Search Tree Construction

**Tree Structure** (from report):
```
                    Root
                   /    \
                  Q      E
                 /        \
                I          I
            /   |   \      |
           V    I  L  M    V
          [13] [12][10][10][10]
          /  \
         V    I
        [12] [11]
        /     |
       V      I
      [11]   [10]
```

**Explanation:**
- First letter branches: Q, E
- Second letter branches: I (most common), L, M, V
- Third letter: Terminal nodes with scores
- Paths from root to leaves represent valid words

**Search Efficiency:**
- Pre-compute tree once for query
- Database search: Traverse tree for each 3-mer
- O(1) lookup instead of O(n) calculation

---

### Question 3: FASTA Longest Common Substring

**Task**: Use FASTA algorithm to find longest matching regions

#### Part A: K-tuple = 1

**Query**: `CACGTTGACAT`

**Database**:
1. `ATGACATTCGAA`
2. `CGATTCGGACA`

**Solution Process:**

**Step 1: Build Query Hash Table (k=1)**
```
Query positions:
1:C, 2:A, 3:C, 4:G, 5:T, 6:T, 7:G, 8:A, 9:C, 10:A, 11:T

Hash table:
A → [2, 8, 10]
C → [1, 3, 9]
G → [4, 7]
T → [5, 6, 11]
```

**Step 2: Search Target 1** `ATGACATTCGAA`

**Match Table** (diagonal offset for each match):
```
Position:  1  2  3  4  5  6  7  8  9 10 11 12
Target:    A  T  G  A  C  A  T  T  C  G  A  A
Query_A:   1  3  1 -2 -4 -4 -2 -3 -8 -6 -9-10  (pos 2)
Query_A:   7  4  4  4 -2  2 -1 -2  6 -3 -3 -4  (pos 8)
Query_A:   9  9  6  4  4  4  3  0 -1 -2       (pos 10)
```

**Diagonal Analysis:**
- Diagonal 4: Matches at (2,6), (3,7), (4,8), (5,9)
  - Corresponds to query "TGAC" matching target "ACAT" at offset 4
  - **Longest match: TGACAT (6 characters)** with shift 4

**Step 3: Search Target 2** `CGATTCGGACA`

**Match Table:**
```
Position:  1  2  3  4  5  6  7  8  9 10 11
Target:    C  G  A  T  T  C  G  G  A  C  A
Query_C:   0  2 -1  1  0 -5 -3 -4 -7 -9 -9  (pos 1)
Query_C:   2  5  5  2  1 -3  0 -1 -1 -7 -3  (pos 3)
Query_C:   8  7  7  6  3  1 -1 -1           (pos 9)
```

**Diagonal Analysis:**
- Diagonal 1: Matches at (3,4), (9,10), (10,11)
  - **Longest match: GACA (4 characters)** with shift 1

**Answer:**
- Target 1: `TGACAT` (6 chars, shift=4)
- Target 2: `GACA` (4 chars, shift=1)
- **Overall winner: TGACAT from target 1**

---

#### Part B: K-tuple = 2

**Query**: `GTTACCACG`

**Database**: `TACGTCGT`

**Solution Process:**

**Step 1: Build Query Hash Table (k=2)**
```
Query 2-tuples:
Position 1: GT
Position 2: TT
Position 3: TA
Position 4: AC
Position 5: CC
Position 6: CA
Position 7: AC
Position 8: CG

Hash table:
AC → [4, 7]
CA → [6]
CC → [5]
CG → [8]
GT → [1]
TA → [3]
TT → [2]
```

**Step 2: Search Target** `TACGTCGT`

**2-tuples in target:**
```
Position 1: TA
Position 2: AC
Position 3: CG
Position 4: GT
Position 5: TC
Position 6: CG
Position 7: GT
```

**Match Table:**
```
Position:  1   2   3   4   5   6   7
Target:    TA  AC  CG  GT  TC  CG  GT
Diagonal:  2   2   5  -3   -   2   6
          (3) (4) (8) (1)       (8) (1)
```

**Diagonal Analysis:**
- Diagonal 2: 
  - TA (pos 1 in target, pos 3 in query): offset = 1-3 = -2
  - AC (pos 2 in target, pos 4 in query): offset = 2-4 = -2
  - CG (pos 6 in target, pos 8 in query): offset = 6-8 = -2
  - **Cluster: TAC (3 chars) at position 1**

- Diagonal -6:
  - GT (pos 4 in target, pos 1 in query): offset = 4-1 = 3 (not -6, correction)

**Recalculation with proper offsets:**
```
TA at target[1] matches query[3]: shift = 1-3 = -2
AC at target[2] matches query[4]: shift = 2-4 = -2  ← Same diagonal!
CG at target[3] matches query[8]: shift = 3-8 = -5

Result: TAC aligns (3 chars, shift 2)
```

**Additionally:**
```
CG at target[3] and [6] match query[8]
GT at target[4] and [7] match query[1]

ACG region: target[2-4] vs query[4-6] (with shift)
```

**Answer:**
- **TAC** (3 characters) with shift=2
- **ACG** (3 characters) with shift=3
- Both are longest matches with k=2

---

### Question 4: Star Alignment Manual Calculation

**Task**: Manually compute MSA using Star algorithm

**Sequences:**
```
A: ACGCTAAC
B: TTGCACATC
C: TCGGTAGATC
D: TCACTGGAC
```

**Scoring**: DNAfull matrix (https://rosalind.info/glossary/dnafull/)
- Match: +5
- Mismatch: -4 (typically)
- Gap open: -10
- Gap extend: -1

**Tool**: EMBOSS Needle (https://www.ebi.ac.uk/Tools/psa/emboss_needle/)

#### Step 1: Pairwise Alignments

**A vs B:**
```
----ACGCTAAC
TTGCAC---ATC

Identity: 3/12 (25%)
Score: 5.0
```

**A vs C:**
```
---ACGCTAAC--
TCGGTAGATC---

Identity: 2/13 (15%)
Score: 4.0
```

**A vs D:**
```
-ACGCTAAC-
TCACTGGAC-

Identity: 5/10 (50%)
Score: 11.0
```

**B vs C:**
```
--TTGCACATC
TCGGTAGATC-

Identity: 4/11 (36%)
Score: 11.5
```

**B vs D:**
```
TTGCACATC
TCACTGGAC

Identity: 5/9 (56%)
Score: 12.0
```

**C vs D:**
```
TCGGTAGATC
TC-ACTGGAC

Identity: 6/10 (60%)
Score: 10.0
```

#### Step 2: Score Matrix

```
     A     B     C     D   | Total
A    -     5     4    11   |  20
B    5     -    11.5  12   |  28.5
C    4    11.5   -    10   |  25.5
D   11    12    10     -   |  33   ← Highest!
```

**Center Selection**: **Sequence D** (total score = 33)

#### Step 3: Progressive Alignment

**Order by score with D:**
1. B (score 12 with D)
2. A (score 11 with D)
3. C (score 10 with D)

**Round 1: Align D with B**
```
D: TCACTGGAC---
B: ---TTGCACATC
```

**Round 2: Add A to alignment**
```
Align A with current D:
D: ----TCACTGGAC---
A: ACGCTAAC--------

Propagate gaps to B:
D: ----TCACTGGAC---
B: -------TTGCACATC
A: ACGCTAAC--------
```

**Round 3: Add C to alignment**
```
Align C with current D:
C needs 8 gaps at start to align with D's position

Initial attempt:
D: --------TCACTGGAC---
C: ----TCGGTAGATC----------

After propagation:
D: --------TCACTGGAC---
B: -----------TTGCACATC
A: ----ACGCTAAC--------
C: ----TCGGTAGATC----------
```

**Remove all-gap columns** (first 4 positions are all gaps):

#### Step 4: Final Alignment

```
D: ----TCACTGGAC---
B: -------TTGCACATC
A: ACGCTAAC--------
C: TCGGTAGATC------
```

**Output in original order (A, B, C, D):**
```
A: ACGCTAAC--------
B: -------TTGCACATC
C: TCGGTAGATC------
D: ----TCACTGGAC---
```

**Analysis:**
- Center (D) has fewest total gaps
- Conserved regions: "TC", "AC", "C" positions
- High divergence reflects distant sequences

---

### Question 5: Star vs ClustalW vs T-Coffee

**Task**: Compare three progressive MSA algorithms

#### Star Alignment

**Method**: Center-star heuristic

**Algorithm:**
1. Compute all pairwise alignments
2. Select center sequence (highest sum of scores)
3. Progressively align sequences to center
4. Apply "once a gap, always a gap" rule

**Advantages:**
✅ **Fast**: Only n pairwise alignments after center selection  
✅ **Simple**: Easy to understand and implement  
✅ **Predictable**: Deterministic given score matrix  
✅ **Good approximation**: If metric satisfies triangle inequality  

**Disadvantages:**
❌ **Greedy**: Early mistakes cannot be corrected  
❌ **Center-dependent**: Poor center choice = poor alignment  
❌ **No phylogeny**: Doesn't consider evolutionary relationships  
❌ **Suboptimal**: Not guaranteed optimal MSA  

**Limitations:**
- Sensitive to scoring scheme choice
- Assumes center represents entire family
- Poorly handles highly divergent sequences
- No iterative refinement in basic version

**Best For:**
- Quick MSA of moderate sequence sets
- Teaching/learning MSA concepts
- Initial alignment for refinement
- Sequences with clear "representative"

---

#### ClustalW

**Method**: Progressive alignment with guide tree and weighted scoring

**Algorithm:**
1. **Pairwise alignment**: Compute all-vs-all distances
2. **Build guide tree**: Neighbor-joining from distance matrix
3. **Progressive alignment**: 
   - Start from tree leaves (most similar pairs)
   - Align pairs/groups moving toward root
   - Use sequence weights to down-weight over-represented groups

**Key Features:**
- **Position-specific gap penalties**: Lower penalties in flexible regions
- **Dynamic substitution matrices**: BLOSUM for distant, PAM for close
- **Sequence weighting**: Prevents bias from redundant sequences
- **Delay divergent sequences**: Align similar sequences first

**Advantages:**
✅ **Phylogeny-aware**: Uses evolutionary relationships  
✅ **Widely used**: Industry standard, well-tested  
✅ **Adaptive scoring**: Position and distance-dependent parameters  
✅ **Handles redundancy**: Sequence weighting reduces bias  
✅ **Global alignment**: Good for sequences of similar length  

**Disadvantages:**
❌ **Global only**: Poor for sequences with very different lengths  
❌ **Error propagation**: Early mistakes cannot be fixed  
❌ **Greedy progressive**: No optimization after initial alignment  
❌ **Moderate speed**: Slower than Star  

**Limitations:**
- Guide tree quality affects final alignment
- Not ideal for divergent sequences with local similarity
- Computational cost: O(n² + nL²) for n sequences of length L
- Cannot realign after initial pass

**Best For:**
- Protein family alignments
- Sequences of similar length
- Well-defined homologs
- Standard bioinformatics workflows

---

#### T-Coffee (Tree-based Consistency Objective Function)

**Method**: Progressive alignment with library of pairwise alignments

**Algorithm:**
1. **Build library**: 
   - Global pairwise alignments (Needleman-Wunsch)
   - Local pairwise alignments (Smith-Waterman)
   - Combine both into consistency library
2. **Extend library**: Add transitivity information
   - If A[i]=B[j] and B[j]=C[k], then support A[i]=C[k]
3. **Build guide tree**: From library scores
4. **Progressive alignment**: Use library scores as weights

**Key Innovation**: **Consistency-based scoring**
```
Standard: Align A-B based only on A-B similarity
T-Coffee: Align A-B considering A-C, C-B information

Example:
A: ACT
B: AGT  
C: ACT

A-B local alignment: ACT vs AGT (C vs G uncertain)
But: A-C has ACT=ACT (C matched)
     C-B has ACT vs AGT (C in C matches T in B position)
     → Support C alignment in A-B
```

**Advantages:**
✅ **Very accurate**: Best accuracy for divergent sequences  
✅ **Handles length variation**: Combines global + local  
✅ **Consistency check**: Uses transitive alignment information  
✅ **Robust**: Less sensitive to guide tree errors  
✅ **Gold standard**: Often used as benchmark  

**Disadvantages:**
❌ **Slow**: O(n³L²) complexity  
❌ **Memory intensive**: Stores all pairwise alignments  
❌ **Still progressive**: No true global optimization  
❌ **Complex**: Harder to understand and modify  

**Limitations:**
- Computationally expensive for large datasets (>100 sequences)
- Memory requirements scale poorly
- Still cannot correct early alignment errors (progressive)
- Overkill for closely related sequences

**Best For:**
- Difficult alignments (divergent sequences)
- Different-length sequences
- Gold-standard alignments for benchmarking
- Publication-quality MSA (accuracy critical)
- Structural biology applications

---

#### Comparison Summary

| Feature | Star | ClustalW | T-Coffee |
|---------|------|----------|----------|
| **Speed** | Fast | Moderate | Slow |
| **Accuracy** | Low-Moderate | Good | Excellent |
| **Method** | Center-star | Phylo progressive | Library + progressive |
| **Complexity** | O(n²L) | O(n²L + nL²) | O(n³L²) |
| **Memory** | Low | Moderate | High |
| **Best for** | Quick/teaching | Standard proteins | Divergent seqs |
| **Phylogeny** | No | Yes (guide tree) | Yes (guide tree) |
| **Refinement** | Optional | No | No (but consistent) |
| **Gap penalty** | Fixed | Position-specific | Library-weighted |
| **Alignment type** | Global | Global | Global + Local |

**Evolutionary Progression:**
```
Star (1986)
  ↓ (add phylogeny)
ClustalW (1994)
  ↓ (add consistency)
T-Coffee (2000)
```

**Use Case Decision Tree:**
```
Need MSA?
├─ Few sequences (<10), speed important?
│  └─ Star ✓
├─ Standard protein family, proven method?
│  └─ ClustalW ✓
├─ Divergent sequences, accuracy critical?
│  └─ T-Coffee ✓
└─ >1000 sequences, memory limited?
   └─ ClustalW or specialized tool
```

**From Report Summary:**
- **ClustalW**: Phylogeny-aware, position-specific penalties, widely used standard
- **Star**: Fast heuristic, good with triangle inequality, greedy
- **T-Coffee**: Combines global+local, consistency checking, most accurate but slowest

---

## 🛠️ Installation & Setup

### Prerequisites
```bash
Python 3.x (no external libraries required)
```

### Project Setup

**Clone Repository:**
```bash
git clone https://github.com/zamirmehdi/Bioinformatics-Course.git
cd "Bioinformatics-Course/3- Multiple Sequence Alignment - DB Search"
```

**Run Implementation:**
```bash
cd src
python main.py
```

### Testing

**Create Test Input:**
```bash
cat > test_input.txt << EOF
4
TYIMREAQYESAQ
TCIVMREAYE
YIMQEVQQER
WRYIAMREQYES
EOF
```

**Run Test:**
```bash
python main.py < test_input.txt
```

**Expected Output:**
```
51
-TYI-MREAQYESAQ
-TCIVMREA-YE---
--YI-MQEVQQER--
WRYIAMRE-QYES--

51
-TYI-MREAQYESAQ
-TCIVMREA-YE---
--YI-MQEVQQER--
WRYIAMRE-QYES--
```

---

## ⚙️ Algorithm Complexity

### Time Complexity

**Star Alignment:**
```
Step 1: Pairwise alignments
  - n × (n-1) / 2 alignments
  - Each O(L²) dynamic programming
  - Total: O(n²L²)

Step 2: Find center
  - O(n²) to sum scores

Step 3: Progressive alignment
  - n-1 alignments to center
  - Each O(L²)
  - Gap propagation: O(nL) per sequence
  - Total: O(nL² + n²L)

Overall: O(n²L²) dominated by pairwise alignments
```

**Block Improvement:**
```
Per iteration:
  - Identify blocks: O(nL)
  - For each block (average length b):
    - Extract: O(nb)
    - Realign: O(n²b²)
    - Calculate score: O(n²b)
  - Total per iteration: O(k × n²b²) where k = number of blocks

Iterations: Usually converges in 2-5 iterations

Overall: O(k × iterations × n²b²)
Typically much less than initial alignment since b << L
```

**Complete Algorithm:**
```
Total: O(n²L²) + O(k × i × n²b²)
     ≈ O(n²L²) for most practical cases
```

### Space Complexity

**Memory Usage:**
```
score_matrix:      O(n²) integers
direction_matrix:  O(L²) per alignment (reused)
aligned_sequences: O(nL) characters
columns:           O(nL) for scoring

Total: O(n²L + nL) = O(n²L)
```

**For n=100 sequences, L=1000:**
```
Score matrix: 100² × 4 bytes = 40 KB
Sequences: 100 × 1000 = 100 KB
Columns: 100 × 1000 = 100 KB

Total: ~250 KB << 512 MB limit ✓
```

### Practical Constraints

**Given Limits:**
- Time: 5 seconds
- Memory: 512 MB

**Feasible Problem Size:**
```
With O(n²L²):
  n=10, L=1000: ~100M operations → <1 sec ✓
  n=20, L=500:  ~100M operations → <1 sec ✓
  n=50, L=100:  ~25M operations  → <0.5 sec ✓
  
Block improvement adds ~20-50% overhead
Still well within time limits for typical inputs
```

---

## 📈 Results Summary

*(From submitted programming and theoretical assignments)*

### Programming Results

**Test Case 1: 4 Protein Sequences**
```
Input:
TYIMREAQYESAQ
TCIVMREAYE
YIMQEVQQER
WRYIAMREQYES

Initial Score: 51
Final Score: 51 (no improvement possible)

Alignment:
-TYI-MREAQYESAQ
-TCIVMREA-YE---
--YI-MQEVQQER--
WRYIAMRE-QYES--
```

**Analysis:**
- Star alignment already optimal
- Conserved "IMRE" and "YE" regions
- Center sequence: WRYIAMREQYES (highest total score)
- Block refinement found no improvements

---

**Test Case 2: 5 DNA Sequences**
```
Input:
TAGCTACCAGGA
CAGCTACCAGG
TAGCTACCAGT
CAGCTATCGCGGC
CAGCTACCAGGA

Initial Score: 240
Final Score: 240 (converged immediately)

Alignment:
TAGCTA-C-CAGGA
CAGCTA-C-CAGG-
TAGCTA-C-CA-GT
CAGCTATCGC-GGC
CAGCTA-C-CAGGA
```

**Analysis:**
- High conservation in "CAGCTA" prefix
- Single insertion in sequence 4 ("TCGC")
- Gap coordination maintains alignment quality
- No further improvement possible

---

### Theoretical Results

**Question 1: Algorithm Comparison**
- Successfully compared DP, FASTA, BLAST
- Key insight: Speed vs sensitivity trade-off
- Understanding of when to use each method

**Question 2: BLAST Word Search**
- Identified 10 words scoring ≥10 with "QIV"
- Constructed search tree for efficient lookup
- Demonstrated understanding of BLAST indexing

**Question 3: FASTA K-tuple Matching**
- **ktup=1**: Found TGACAT (6 chars) as longest match
- **ktup=2**: Found TAC and ACG (3 chars each)
- Correct diagonal analysis and offset calculations

**Question 4: Manual Star Alignment**
- Correctly identified D as center (score=33)
- Performed progressive alignment with gap propagation
- Final 4-sequence MSA with proper column structure

**Question 5: Algorithm Comparison**
- Comprehensive analysis of Star, ClustalW, T-Coffee
- Understanding of progressive alignment evolution
- Recognition of accuracy vs speed trade-offs

---

## 🗂️ Project Structure

```
3- Multiple Sequence Alignment - DB Search/
├── docs/
│   ├── Programming Instruction MSA.pdf  # Coding assignment (Persian)
│   ├── Theoretical/
│   │   ├── Instruction.pdf              # Theory questions (Persian)
│   │   ├── Report.pdf                   # Completed solutions (Persian)
│   │   └── cstar.pdf                    # Star alignment diagram
├── src/
│   └── main.py                          # Star alignment implementation
└── README.md                            # This documentation
```

**File Descriptions:**

**`main.py`** (343 lines):
- Global alignment function (from Project 2)
- Star alignment core algorithm
- Block identification and refinement
- Score calculation (sum-of-pairs)
- Main execution with input/output handling

**`Programming Instruction MSA.pdf`**:
- Star alignment specification
- Block-based improvement rules
- Scoring scheme details
- Input/output format
- Test cases and examples

**`Theoretical/Instruction.pdf`**:
- 5 comprehensive theory questions
- FASTA/BLAST/DP comparison
- Manual calculation exercises
- Algorithm analysis tasks

**`Theoretical/Report.pdf`**:
- Complete solutions with calculations
- FASTA diagonal tables
- Star alignment step-by-step
- Algorithm comparison essays

**`cstar.pdf`**:
- Visual diagram of center-star method
- Illustrates progressive alignment concept

---

## 🎓 Key Concepts Covered

### Multiple Sequence Alignment
- ✅ Sum-of-pairs scoring
- ✅ Progressive alignment strategies
- ✅ Center-star heuristic
- ✅ "Once a gap, always a gap" rule
- ✅ Block-based iterative refinement

### Database Search Algorithms
- ✅ Dynamic programming (Smith-Waterman)
- ✅ FASTA k-tuple hashing and extension
- ✅ BLAST word neighborhoods and two-hit
- ✅ Speed vs sensitivity trade-offs
- ✅ Indexing and search tree construction

### Alignment Algorithms
- ✅ Star alignment (center-star)
- ✅ ClustalW (phylogeny-guided progressive)
- ✅ T-Coffee (consistency-based library)
- ✅ Progressive vs iterative methods
- ✅ Guide tree construction

### Computational Techniques
- ✅ Dynamic programming reuse
- ✅ Hash table lookup
- ✅ Diagonal analysis
- ✅ Iterative refinement
- ✅ Convergence detection

---

## 🎯 Learning Outcomes

After completing this project, students can:

### Programming Skills
✅ Implement center-star MSA algorithm  
✅ Apply iterative refinement techniques  
✅ Calculate sum-of-pairs alignment scores  
✅ Handle gap propagation across sequences  
✅ Detect and optimize suboptimal regions  

### Bioinformatics Knowledge
✅ Understand MSA fundamentals and applications  
✅ Compare database search algorithms (FASTA, BLAST, DP)  
✅ Recognize heuristic vs exhaustive methods  
✅ Apply appropriate algorithm for given problem  
✅ Interpret MSA quality and biological significance  

### Analytical Thinking
✅ Analyze algorithm complexity (time and space)  
✅ Evaluate trade-offs (speed vs accuracy)  
✅ Manually trace algorithm execution  
✅ Optimize alignment through block refinement  
✅ Compare progressive alignment strategies  

### Practical Skills
✅ Use bioinformatics web tools (EMBOSS)  
✅ Interpret alignment scores and statistics  
✅ Construct search trees and hash tables  
✅ Perform diagonal analysis for substring matching  
✅ Validate implementations with manual calculations  

---

## ℹ️ Project Information

**Project:** Multiple Sequence Alignment & Database Search  
**Author:** Amirmehdi Zarrinnezhad  
**Course:** Bioinformatics  
**University:** Amirkabir University of Technology (Tehran Polytechnic) - Fall 2022  
**Language:** Python 3.x (Programming), English (README), Persian (Reports)  
**GitHub Link:** [3- Multiple Sequence Alignment - DB Search](https://github.com/zamirmehdi/Bioinformatics-Course/tree/main/3-%20Multiple%20Sequence%20Alignment%20-%20DB%20Search)

<div align="center">

**Part of Bioinformatics Course Projects**

[1: Basic Biology](../1-%20Basic%20biology) | [2: Sequence Alignment](../2-%20Pairwise%20Sequence%20Alignment) | [3: MSA & DB Search](.) | [4: Profile HMM](../4-%20Profile%20-%20Hidden%20Markov%20model) | [5: Phylogenetic Trees](../5-%20Phylogenetic%20Trees) | [Final: Virus Classification](../Virus%20Classification%20(Final%20Project))

</div>

---

## 📧 Contact

Questions or collaborations? Feel free to reach out!  
📧 Email: amzarrinnezhad@gmail.com  
💬 Open an [Issue](https://github.com/zamirmehdi/Bioinformatics-Course/issues)  
🌐 GitHub: [@zamirmehdi](https://github.com/zamirmehdi)

---

<div align="center">

[⬆ Back to Main Repository](https://github.com/zamirmehdi/Bioinformatics-Course)

</div>

<p align="right">(<a href="#top">back to top</a>)</p>

<div align="center">

⭐ **If you found this project helpful, please consider giving it a star!** ⭐

*Amirmehdi Zarrinnezhad*

</div>
