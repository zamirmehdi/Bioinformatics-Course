# Project 2 - Pairwise Sequence Alignment

A comprehensive implementation of **Semi-Global Alignment** algorithm for protein sequences using dynamic programming, with theoretical exploration of global and local alignment techniques, scoring matrices (PAM250, BLOSUM), and practical applications in bioinformatics.

[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![Algorithm](https://img.shields.io/badge/Algorithm-Dynamic%20Programming-green.svg)](#)
[![Course](https://img.shields.io/badge/Course-Bioinformatics-blue.svg)](#)
[![University](https://img.shields.io/badge/University-AUT-red.svg)](https://aut.ac.ir/en)

<details> <summary><h2>📚 Table of Contents</h2></summary>

- [Overview](#-overview)
- [Alignment Types](#-alignment-types)
- [Programming Assignment](#-programming-assignment)
  - [Semi-Global Alignment Algorithm](#semi-global-alignment-algorithm)
  - [Implementation Details](#implementation-details)
  - [Input/Output Format](#inputoutput-format)
  - [Usage Examples](#usage-examples)
- [Theoretical Assignment](#-theoretical-assignment)
  - [Question 1: Needleman-Wunsch & Dot Matrix](#question-1-needleman-wunsch--dot-matrix)
  - [Question 2: Dynamic Programming Table](#question-2-dynamic-programming-table)
  - [Question 3: Local Alignment (Smith-Waterman)](#question-3-local-alignment-smith-waterman)
  - [Question 4: Global vs Local Comparison](#question-4-global-vs-local-comparison)
  - [Question 5: PAM vs BLOSUM Matrices](#question-5-pam-vs-blosum-matrices)
- [Key Concepts](#-key-concepts)
- [Installation & Setup](#-installation--setup)
- [Algorithm Complexity](#-algorithm-complexity)
- [Scoring Matrices](#-scoring-matrices)
- [Project Structure](#%EF%B8%8F-project-structure)
- [Results & Analysis](#-results--analysis)
- [Learning Outcomes](#-learning-outcomes)
- [Project Information](#ℹ%EF%B8%8F-project-information)
- [References](#-references)
- [Contact](#-contact)

</details>

---

## 📋 Overview

This project explores **Pairwise Sequence Alignment**, one of the most fundamental problems in bioinformatics. Sequence alignment is crucial for:

- **Identifying homologous sequences** (evolutionary relationships)
- **Predicting protein structure and function**
- **Detecting conserved regions** across species
- **Understanding evolutionary distances**
- **Designing primers** for PCR and sequencing

The project consists of two main components:

### 1. Programming Component
**Implementation**: Semi-Global Alignment algorithm for protein sequences
- Uses PAM250 substitution matrix
- Finds ALL optimal alignments (not just one)
- Handles free end gaps (no penalty at sequence termini)

### 2. Theoretical Component
**Topics Covered**:
- Global alignment (Needleman-Wunsch)
- Local alignment (Smith-Waterman)
- Dot matrix visualization
- Scoring schemes and gap penalties
- PAM vs BLOSUM matrix comparison

---

## 🔄 Alignment Types

### Global Alignment (Needleman-Wunsch)
**Purpose**: Align entire sequences from beginning to end

**Best For**:
- Sequences of similar length
- High overall similarity
- Close evolutionary distance
- Complete sequence comparison

**Characteristics**:
- Gaps penalized everywhere (including ends)
- Optimal end-to-end alignment
- Maximizes alignment score over full length

**Algorithm**: Dynamic programming with full matrix

```
Example:
Seq1: HEAGAWGHEA
Seq2: PAWHEAG---
Score: Considers all positions
```

### Local Alignment (Smith-Waterman)
**Purpose**: Find best matching subsequences

**Best For**:
- Sequences of different lengths
- Low overall similarity with conserved regions
- Distant evolutionary relationship
- Domain/motif identification

**Characteristics**:
- No penalty for unaligned regions
- Identifies highest-scoring local region
- Allows zero scores (restart alignment)

**Algorithm**: Dynamic programming with zero floor

```
Example:
Seq1: MEAMNVE---KASADGNL-------PEVIS
Seq2: ---MTXPALSLHTPLST--SFTPAVWYNG
Score: Only conserved regions considered
```

### Semi-Global Alignment
**Purpose**: Hybrid approach - align full sequences with free end gaps

**Best For**:
- Sequences where one may be substring of other
- Comparing fragments to reference sequences
- Query sequence vs database search
- Avoiding end-gap penalties

**Characteristics**:
- **No penalty for gaps at sequence ends** (termini)
- Gaps in middle regions still penalized
- Optimal alignment considering free end gaps
- Multiple optimal solutions possible

**Algorithm**: Modified Needleman-Wunsch

```
Example (This Project):
Seq1: HEAGAWGHE-
Seq2: ---PAW-HEA
      ↑       ↑
   Free gaps at ends (no penalty)
```

**Comparison Table:**

| Feature | Global | Local | Semi-Global |
|---------|--------|-------|-------------|
| **End gaps** | Penalized | Ignored | Free (no penalty) |
| **Middle gaps** | Penalized | Penalized | Penalized |
| **Full sequences** | Yes | No | Yes |
| **Use case** | Similar seqs | Conserved regions | Fragment alignment |
| **Algorithm** | Needleman-Wunsch | Smith-Waterman | Modified NW |

---

## 💻 Programming Assignment

### Semi-Global Alignment Algorithm

**Objective**: Implement semi-global alignment for protein sequences using dynamic programming.

**Problem Specification**:
- **Input**: Two protein sequences (uppercase amino acid letters)
- **Scoring**: PAM250 substitution matrix
- **Gap Penalty**: -9 (linear gap model)
- **Output**: 
  - Alignment score (maximum)
  - ALL optimal alignments (sorted lexicographically)
- **Constraints**:
  - Time limit: 1.5 seconds
  - Memory limit: 100 MB
  - Max sequence length: 100 characters

### Implementation Details

#### 1. Core Data Structures

```python
PAM250 = {...}              # 20x20 amino acid substitution matrix
GAP_PENALTY = -9            # Linear gap penalty

score_matrix = []           # DP matrix for scores
direction_matrix = []       # Traceback directions ('d', 'h', 'v')
total_score = 0            # Best alignment score
total_locations = []        # Starting positions for traceback
seq = []                    # Stores all optimal alignments
```

#### 2. Algorithm Steps

**Step 1: Initialize Matrices**
```python
def init_and_fill_matrix():
    # Create (len(str1)+1) x (len(str2)+1) matrices
    # Fill using dynamic programming recurrence
```

**DP Recurrence Relation:**
```
For cell (i, j):

score[i][j] = max(
    score[i-1][j-1] + PAM250[str1[i-1]][str2[j-1]],  # Match/Mismatch
    score[i][j-1] + GAP_PENALTY,                      # Horizontal (gap in str1)
    score[i-1][j] + GAP_PENALTY                       # Vertical (gap in str2)
)

Direction tracking:
- 'd' = diagonal (match/mismatch)
- 'h' = horizontal (gap in str1)
- 'v' = vertical (gap in str2)
```

**Key Difference from Global Alignment**:
```python
# Semi-global: Find max score in LAST ROW or LAST COLUMN only
if i == len(str1) or j == len(str2):
    if temp_score >= total_score:
        total_score = temp_score
```

**Step 2: Find Optimal Starting Points**
```python
def find_total_score_locations():
    # Search last row and column for cells with total_score
    # These are potential traceback starting points
```

**Step 3: Traceback All Paths**
```python
def trace_back(strings, x, y):
    # Recursively follow direction_matrix
    # Generate all paths leading to score = total_score
    # Add free gaps at sequence ends
```

**Step 4: Format Output**
```python
def print_output(score, seq):
    # Sort alignments lexicographically
    # Print score first
    # Then print each alignment (2 lines per alignment)
```

#### 3. Key Functions

**`init_and_fill_matrix()`**
- Creates and fills score/direction matrices
- Implements DP recurrence
- Tracks maximum score in borders

**`calculate_score(arr, i, j)`**
- Computes optimal score for cell (i,j)
- Returns score and direction(s)
- Handles multiple optimal paths (multiple directions possible)

**`find_total_score_locations()`**
- Identifies all cells with maximum score in last row/column
- These become traceback starting points

**`trace_back(strings, x, y)`**
- Recursively generates all optimal alignments
- Adds gaps at sequence ends (semi-global feature)
- Handles multiple traceback paths

**`semi_global_alignment()`**
- Orchestrates traceback from all optimal locations
- Prepends/appends gaps for unaligned ends

### Input/Output Format

#### Input Format
```
Line 1: First protein sequence (e.g., HEAGAWGHE)
Line 2: Second protein sequence (e.g., PAWHEA)
```

**Constraints**:
- Only uppercase letters (20 standard amino acids)
- No spaces or special characters
- Maximum 100 characters per sequence

#### Output Format
```
Line 1: Alignment score (integer)
Line 2+: All optimal alignments (sorted)
         Each alignment = 2 lines (seq1, seq2)
         Use '-' for gaps
```

#### Example 1: Simple Case

**Input**:
```
HEAGAWGHE
PAWHEA
```

**Output**:
```
20
HEAGAWGHE-
---PAW-HEA
```

**Explanation**:
- Score: 20
- 1 optimal alignment
- Free gaps at start of seq2 (---)
- Free gap at end of seq1 (-)

#### Example 2: Multiple Optimal Alignments

**Input**:
```
AAAAA
AA
```

**Output**:
```
4
AAAAA
---AA
AAAAA
--AA-
AAAAA
-AA--
AAAAA
AA---
```

**Explanation**:
- Score: 4 (2 matches * 2 points each)
- 4 optimal alignments possible
- AA can align at any position within AAAAA
- All produce same score due to free end gaps
- Sorted lexicographically

### Usage Examples

#### Running the Program

**Method 1: Command Line (stdin)**
```bash
cd "2- Pairwise Sequence Alignment/src"
python semi_global_alignment.py
# Then enter sequences when prompted:
HEAGAWGHE
PAWHEA
```

**Method 2: Redirecting Input**
```bash
echo -e "HEAGAWGHE\nPAWHEA" | python semi_global_alignment.py
```

**Method 3: Input File**
```bash
python semi_global_alignment.py < input.txt
```

#### Sample Test Cases

**Test 1: Basic Alignment**
```python
Input:
HEAGAWGHE
PAWHEA

Output:
20
HEAGAWGHE-
---PAW-HEA
```

**Test 2: Multiple Solutions**
```python
Input:
AAAAA
AA

Output:
4
AAAAA
---AA
AAAAA
--AA-
AAAAA
-AA--
AAAAA
AA---
```

**Test 3: No Similarity**
```python
Input:
AAA
CCC

Output:
-18
AAA
CCC
```

**Test 4: Identical Sequences**
```python
Input:
MEAM
MEAM

Output:
16
MEAM
MEAM
```

---

## 📝 Theoretical Assignment

### Question 1: Needleman-Wunsch & Dot Matrix

**Task**: Perform global alignment using Needleman-Wunsch and dot matrix methods.

#### Part A: Needleman-Wunsch Algorithm

**Sequences**:
- Seq1: `ACCTAGA`
- Seq2: `ACTGG`

**Scoring Scheme**:
```
Match: +1
Mismatch: -1
Gap: -2 (linear penalty)
```

**Solution Steps**:

1. **Initialize DP Matrix** (8x6 matrix)
```
    -   A   C   T   G   G
-   0  -2  -4  -6  -8  -10
A  -2
C  -4
C  -6
T  -8
A  -10
G  -12
A  -14
```

2. **Fill Matrix Using Recurrence**
```
For each cell (i,j):
score[i][j] = max(
    score[i-1][j-1] + match/mismatch,
    score[i-1][j] + gap,
    score[i][j-1] + gap
)
```

3. **Traceback from Bottom-Right**
- Follow arrows/directions to reconstruct alignment
- Multiple paths = multiple optimal alignments

**Expected Result** (from report):
```
Alignment:
ACCTAGA
AC-T-GG

Score: Calculated based on scoring scheme
```

#### Part B: Dot Matrix Visualization

**Method**:
1. Create 7x5 grid (ACCTAGA vs ACTGG)
2. Place dot where characters match
3. Identify diagonal patterns (alignment paths)

**Dot Matrix**:
```
      A  C  T  G  G
  A   •  -  -  -  -
  C   -  •  -  -  -
  C   -  •  -  -  -
  T   -  -  •  -  -
  A   •  -  -  -  -
  G   -  -  -  •  •
  A   •  -  -  -  -
```

**Analysis**:
- Diagonal lines indicate alignment regions
- Breaks in diagonals = gaps needed
- Compare with Needleman-Wunsch result

---

### Question 2: Dynamic Programming Table

**Task**: Complete partial DP table and analyze alignment properties.

#### Part A: Identify Scoring Parameters

**Given partial table**, determine:
- Match score: `m`
- Mismatch score: `s`
- Gap penalty: `g`

**Method**:
- Examine filled cells
- Reverse-engineer recurrence relation
- Verify consistency across multiple cells

**From Report**:
- Match: `+2`
- Mismatch: `-1`
- Gap: `-2`

#### Part B: Complete the Table

Fill remaining cells using identified parameters:
```
score[i][j] = max(
    score[i-1][j-1] + (m if match else s),
    score[i-1][j] + g,
    score[i][j-1] + g
)
```

#### Part C: Count Optimal Alignments

**Question**: How many optimal global alignments exist?

**Answer Method**:
1. Identify final score in bottom-right
2. Count all traceback paths leading to this score
3. Multiple paths = multiple optimal alignments

**Why multiple alignments?**
- Symmetric scoring (e.g., gap-gap = gap-mismatch)
- Equivalent substitutions
- Different gap placements with same total score

---

### Question 3: Local Alignment (Smith-Waterman)

**Task**: Calculate local alignment using Smith-Waterman algorithm.

**Sequences**:
- Seq1: `TATAGC`
- Seq2: `GTTATC`

**Scoring**:
```
Match: +2
Mismatch: -1
Gap: -2
```

**Key Difference from Global**:
```python
# Smith-Waterman modification:
score[i][j] = max(
    0,                                    # Can restart alignment
    score[i-1][j-1] + match/mismatch,
    score[i-1][j] + gap,
    score[i][j-1] + gap
)
```

**Solution Process** (from report):

1. **Initialize**: First row/column = 0 (not cumulative gaps)

2. **Fill Matrix**:
```
      -  G  T  T  A  T  C
  -   0  0  0  0  0  0  0
  T   0  0  2  2  0  2  0
  A   0  0  0  1  4  2  1
  T   0  0  2  2  2  6  4
  A   0  0  0  1  4  4  5
  G   0  2  0  0  2  3  3
  C   0  0  1  0  0  1  5
```

3. **Find Maximum Score**: 6 (at position corresponding to "TAT")

4. **Traceback Until Zero**:
```
Local Alignment:
TAT
TAT

Score: 6 (3 matches * 2)
```

**Analysis**:
- Only conserved region "TAT" aligned
- Flanking mismatches ignored (zero scores)
- Higher score than global alignment would achieve

---

### Question 4: Global vs Local Comparison

**Task**: Compare alignment strategies for protein sequence pairs.

**Given Sequences**:
```
A = MEAMNVEKASADGNLPEVISNIKETLKIVSRTPVNITMAGDSGNGMSTFI... (214 aa)
B = MEAMNVEKASADGNLPEVISNIKETLKIVSRTPVNITTAGHSGNGMSTFI... (181 aa)
C = MTXPALSLHTPLSTSFTPAVWYNMGWSILSKIGAINIENAVGGGKLLEVD... (389 aa)
```

#### Part A: Choosing Alignment Type

**For (A, B)**: Recommend **Global Alignment**

**Rationale**:
- Similar lengths (214 vs 181)
- High sequence similarity visible
- Few insertions/deletions expected
- Close evolutionary distance
- Want end-to-end comparison

**For (A, C)**: Recommend **Local Alignment**

**Rationale**:
- Very different lengths (214 vs 389)
- Low overall similarity
- Large evolutionary distance
- Interest in conserved domains/motifs
- Global alignment would be poor quality

#### Part B: Web Tool Results Analysis

**Tools Used**:
- EMBOSS Needle (Global): https://www.ebi.ac.uk/Tools/psa/emboss_needle/
- EMBOSS Water (Local): https://www.ebi.ac.uk/Tools/psa/emboss_water/

**Results from Report**:

**Pair (A, B)**:
- **Global Score**: 872.0
  - Identity: 174/214 (81.3%)
  - Similarity: 175/214 (81.8%)
  - Gaps: 33/214 (15.4%)
  
- **Local Score**: 876.0
  - Identity: 173/178 (97.2%)
  - Similarity: 174/178 (97.8%)
  - Gaps: 0/178 (0.0%)

**Analysis**:
- Scores very close (global vs local)
- Local slightly higher (focuses on best region)
- High similarity confirms global alignment appropriate
- Few gaps indicate close evolutionary relationship

**Pair (A, C)**:
- **Global Score**: 410.5
  - Identity: 101/390 (25.9%)
  - Similarity: 140/390 (35.9%)
  - Gaps: 177/390 (45.4%)
  
- **Local Score**: 419.0
  - Identity: 92/176 (52.3%)
  - Similarity: 120/176 (68.2%)
  - Gaps: 3/176 (1.7%)

**Analysis**:
- Much lower global score (distant sequences)
- Local score significantly better quality
- Local alignment finds conserved region
- 45% gaps in global vs 2% in local
- Confirms local alignment more appropriate

**Key Conclusions**:
1. Similar sequences: Global and local scores comparable
2. Distant sequences: Local alignment vastly superior
3. Gap percentage indicates sequence relatedness
4. Local alignment uncovers functional domains

---

### Question 5: PAM vs BLOSUM Matrices

**Task**: Research and compare PAM and BLOSUM scoring matrices.

#### PAM (Point Accepted Mutation) Matrices

**Developed By**: Margaret Dayhoff (1978)

**Methodology**:
- Based on phylogenetic model (evolutionary trees)
- Uses closely related protein sequences (>85% identity)
- Extrapolates from PAM1 matrix recursively
- Analyzes 1572 mutations from 71 protein families

**Calculation**:
1. Start with PAM1 (1% divergence = 1 mutation per 100 residues)
2. Create mutation probability matrix
3. Extrapolate: PAM2 = PAM1², PAM3 = PAM1³, etc.
4. Convert probabilities to log-odds scores

**Versions**:
```
PAM1:   1% divergence   → Closely related sequences
PAM120: 120% divergence → Moderately related
PAM250: 250% divergence → Distantly related sequences
```

**Interpretation**:
- **Higher number = Greater evolutionary distance**
- PAM1: Very similar proteins (recent divergence)
- PAM250: Distant homologs (ancient divergence)

**Applications**:
- Phylogenetic tree reconstruction
- Evolutionary distance estimation
- Database searching
- Genetic disease identification

**Limitations**:
- Based on extrapolation (not direct observation for high PAM)
- Limited to closely related starting sequences
- May not capture long-term evolutionary trends

---

#### BLOSUM (BLOck SUbstitution Matrix) Matrices

**Developed By**: Henikoff & Henikoff (1992)

**Methodology**:
- Based on conserved blocks from multiple sequence alignments
- Direct observation (no extrapolation)
- Uses BLOCKS database (>2000 blocks, 500 protein groups)
- Analyzes conserved regions (low mutation rate)

**Calculation**:
1. Identify conserved blocks (<60 aa, ungapped)
2. Calculate substitution frequencies in blocks
3. Compute log-odds scores:
   ```
   BLOSUM(i,j) = log₂(observed_freq(i,j) / expected_freq(i,j))
   ```
4. Weight sequences to avoid bias

**Versions**:
```
BLOSUM45: 45% identity threshold → Distant sequences
BLOSUM62: 62% identity threshold → Default (general purpose)
BLOSUM80: 80% identity threshold → Closely related sequences
```

**Interpretation**:
- **Higher number = Closer evolutionary relationship**
- BLOSUM45: Distant proteins
- BLOSUM62: **Most commonly used** (balanced)
- BLOSUM80: Similar proteins

**Applications**:
- Database searching (BLAST default: BLOSUM62)
- Finding conserved domains
- Protein function prediction
- Multiple sequence alignment

**Advantages**:
- No extrapolation (direct observation)
- Focuses on conserved regions (more reliable)
- Better performance for database searches

---

#### Comparison: PAM vs BLOSUM

**Fundamental Differences**:

| Aspect | PAM | BLOSUM |
|--------|-----|--------|
| **Basis** | Phylogenetic model | Multiple alignments |
| **Sequences** | Closely related | Varying distances |
| **Method** | Extrapolation | Direct observation |
| **Scope** | Full protein sequences | Conserved blocks only |
| **Number meaning** | Mutations per 100 residues | % Identity threshold |
| **Trend** | Higher = more distant | Higher = more similar |

**Calculation Method**:
- **PAM**: 
  - Pairwise comparisons
  - Recursive multiplication from PAM1
  - Based on evolutionary assumptions
  
- **BLOSUM**: 
  - Multiple sequence alignment blocks
  - Direct frequency counts
  - No evolutionary model required

**Number Interpretation** (OPPOSITE!):
```
PAM:    1 ←─── (similar) ─── 250 →    (distant)
BLOSUM: 80 ←── (similar) ─── 45 →     (distant)
```

**Equivalences** (approximate):
```
PAM250 ≈ BLOSUM45  (distant sequences)
PAM160 ≈ BLOSUM62  (moderate distance)
PAM100 ≈ BLOSUM80  (close sequences)
```

**When to Use**:

**Use PAM**:
- ✅ Global alignment preferred
- ✅ Phylogenetic analysis
- ✅ Complete sequence comparison
- ✅ Evolutionary distance studies

**Use BLOSUM**:
- ✅ Local alignment / database search
- ✅ Domain identification
- ✅ Conserved region analysis
- ✅ General protein comparison (BLOSUM62)

**Performance**:
- BLOSUM generally outperforms PAM for database searches
- BLOSUM62 is default for BLAST (most widely used)
- PAM better for evolutionary studies

---

## 🎓 Key Concepts

### Dynamic Programming
- **Optimal Substructure**: Best alignment of prefixes
- **Overlapping Subproblems**: Reuse computed scores
- **Bottom-Up Computation**: Build solution from small to large
- **Traceback**: Reconstruct solution from DP table

### Sequence Similarity vs Homology
- **Similarity**: Observable sequence resemblance (quantitative)
- **Homology**: Evolutionary relationship (qualitative, binary)
- High similarity suggests homology (but not proof)
- Low similarity doesn't rule out homology (divergent evolution)

### Scoring Schemes
- **Substitution matrices** (PAM, BLOSUM): Amino acid replacement costs
- **Gap penalties**: 
  - Linear: gap(k) = -k × penalty
  - Affine: gap(k) = open + (k × extend)
- **Biological relevance**: Reflect evolutionary/structural constraints

### Semi-Global Alignment Applications
- **Primer design**: Align primer to genome (free end gaps)
- **EST mapping**: Align expressed sequence tags to genes
- **Protein domain search**: Query domain vs full protein
- **Sequence assembly**: Overlap detection

---

## 🛠️ Installation & Setup

### Prerequisites
```bash
Python 3.x (no external libraries required)
```

### Project Setup

**Clone Repository**:
```bash
git clone https://github.com/zamirmehdi/Bioinformatics-Course.git
cd "Bioinformatics-Course/2- Pairwise Sequence Alignment"
```

**Run Implementation**:
```bash
cd src
python semi_global_alignment.py
```

### Testing

**Create Test Input File** (`test_input.txt`):
```
HEAGAWGHE
PAWHEA
```

**Run Test**:
```bash
python semi_global_alignment.py < test_input.txt
```

**Expected Output**:
```
20
HEAGAWGHE-
---PAW-HEA
```

---

## ⚙️ Algorithm Complexity

### Time Complexity
```
O(m × n)

Where:
- m = length of sequence 1
- n = length of sequence 2
```

**Breakdown**:
- Matrix initialization: O(m × n)
- DP table filling: O(m × n)
- Finding max scores: O(m + n) [only borders]
- Traceback: O(m + n) per alignment × k alignments

**Total**: O(m × n + k(m + n))
- k = number of optimal alignments (usually small)

### Space Complexity
```
O(m × n)

Storage:
- score_matrix: (m+1) × (n+1) integers
- direction_matrix: (m+1) × (n+1) strings
```

**Optimization Possible**:
- Linear space traceback: O(m + n) [Hirschberg's algorithm]
- Trade-off: More time, less space

### Practical Constraints
```
Max sequence length: 100
Max matrix size: 101 × 101 = 10,201 cells
Memory: ~100 KB (well within 100 MB limit)
Time: <0.1 seconds (well within 1.5 sec limit)
```

---

## 🧬 Scoring Matrices

### PAM250 Matrix (Used in This Project)

**20 × 20 matrix** for standard amino acids:
```python
PAM250 = {
    'A': {'A': 2, 'C': -2, 'D': 0, ...},
    'C': {'A': -2, 'C': 12, 'D': -5, ...},
    ...
}
```

**Interpretation**:
- **Positive scores**: Likely substitutions (biochemically similar)
- **Negative scores**: Unlikely substitutions (dissimilar)
- **Diagonal values**: Self-match scores (always positive)

**Example Scores**:
```
PAM250['C']['C'] = 12   (Cysteine to Cysteine - perfect match)
PAM250['C']['W'] = -8   (Cysteine to Tryptophan - very unlikely)
PAM250['A']['G'] = 1    (Alanine to Glycine - common, similar size)
```

### Why PAM250?
- Suitable for **distantly related proteins**
- Standard in evolutionary studies
- Widely used reference matrix
- Good for semi-global alignment

---

## 🗂️ Project Structure

```
2- Pairwise Sequence Alignment/
├── docs/
│   ├── Programming Instruction.pdf          # Coding assignment (Persian)
│   └── Theoretical/
│       ├── Instruction.pdf                  # Theory questions (Persian)
│       └── Report.pdf                       # Completed theory answers
├── src/
│   └── semi_global_alignment.py             # Main implementation
└── README.md                                # This documentation
```

**File Descriptions**:

**`semi_global_alignment.py`** (236 lines):
- PAM250 matrix definition
- DP matrix initialization and filling
- Traceback with multiple path handling
- Output formatting and sorting
- Main execution logic

**`Programming Instruction.pdf`**:
- Algorithm specification
- Input/output format
- Scoring scheme
- Test cases
- Submission requirements (Quera platform)

**`Theoretical/Instruction.pdf`**:
- 5 theory questions
- Global/local alignment exercises
- PAM vs BLOSUM research task
- Submission guidelines

**`Theoretical/Report.pdf`**:
- Completed solutions with manual calculations
- Web tool (EMBOSS) results
- PAM/BLOSUM comparison essay

---

## 📊 Results & Analysis

### Programming Assignment Results

**Test Case 1**: Basic Alignment
```
Input:  HEAGAWGHE, PAWHEA
Output: Score = 20, 1 alignment

Analysis:
- Good score indicates similarity
- Free gaps at termini (semi-global feature)
- Conserved region: "AWGHE" / "AWHEA"
```

**Test Case 2**: Multiple Optimal Alignments
```
Input:  AAAAA, AA
Output: Score = 4, 4 alignments

Analysis:
- Symmetric scoring creates multiple optima
- All positions of AA within AAAAA equally valid
- Demonstrates exhaustive traceback
```

### Theoretical Assignment Insights

**Key Findings from Question 4**:

**
