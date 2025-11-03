# Project 2 - Pairwise Sequence Alignment

A comprehensive implementation of **Semi-Global Alignment** algorithm for protein sequences using dynamic programming, PAM250 substitution matrix, and complete traceback for finding all optimal alignments.

[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![Algorithm](https://img.shields.io/badge/Algorithm-Dynamic%20Programming-green.svg)](#)
[![Course](https://img.shields.io/badge/Course-Bioinformatics-blue.svg)](#)
[![University](https://img.shields.io/badge/University-AUT-red.svg)](https://aut.ac.ir/en)

<details> <summary><h2>📚 Table of Contents</h2></summary>

- [Overview](#-overview)
- [Theoretical Background](#-theoretical-background)
  - [Alignment Types](#alignment-types)
  - [Semi-Global Alignment](#semi-global-alignment)
  - [Scoring Matrices](#scoring-matrices)
- [Implementation](#-implementation)
  - [Algorithm Details](#algorithm-details)
  - [Key Features](#key-features)
  - [Code Structure](#code-structure)
- [Installation](#-installation)
- [Usage](#-usage)
  - [Input Format](#input-format)
  - [Output Format](#output-format)
  - [Examples](#examples)
- [Theoretical Assignment](#-theoretical-assignment)
- [Performance Analysis](#-performance-analysis)
- [Project Structure](#%EF%B8%8F-project-structure)
- [Key Concepts Covered](#-key-concepts-covered)
- [Learning Outcomes](#-learning-outcomes)
- [Results Summary](#-results-summary)
- [Project Information](#ℹ%EF%B8%8F-project-information)
- [Contact](#-contact)

</details>

---

## 📋 Overview

This project implements **Semi-Global Alignment** (also known as **End-Gap Free Alignment**) for comparing protein sequences. Unlike global alignment, semi-global alignment allows gaps at the beginning and end of sequences without penalty, making it ideal for comparing sequences of different lengths or finding subsequence similarities.

**Project Components:**
- ✅ **Programming**: Python implementation of semi-global alignment
- ✅ **Theoretical**: Manual calculations using Needleman-Wunsch, dot matrix, local alignment
- ✅ **Analysis**: Comparison of PAM and BLOSUM scoring matrices

**Key Features:**
- Dynamic programming with O(mn) time complexity
- PAM250 substitution matrix for protein scoring
- Linear gap penalty model (gap = -9)
- Complete traceback to find **ALL** optimal alignments
- Lexicographically sorted output

---

## 📚 Theoretical Background

### Alignment Types

**1. Global Alignment (Needleman-Wunsch)**
- Aligns entire length of both sequences
- Penalizes gaps at sequence ends
- Best for: Similar-length, closely related sequences
- Use case: Comparing orthologs with similar evolutionary distance

**2. Local Alignment (Smith-Waterman)**
- Finds best matching subsequences
- No penalty for unaligned regions
- Best for: Sequences with conserved domains but different lengths
- Use case: Finding functional motifs or domains

**3. Semi-Global Alignment (This Project)**
- Aligns sequences end-to-end without penalizing end gaps
- Gaps at start/end are free, internal gaps penalized
- Best for: Finding subsequence within longer sequence
- Use case: Query sequence matching, sequence assembly overlaps

**Comparison:**
```
Sequence A: HEAGAWGHE
Sequence B: PAWHEA

Global:     HEAGAWGHE-     (penalizes end gap)
            ---PAW-HEA

Semi-Global: HEAGAWGHE-    (free end gaps)
             ---PAW-HEA

Local:       AGAWGH        (only best matching region)
             --AW-H
```

### Semi-Global Alignment

**Definition:**
Semi-global alignment is a variant of global alignment where gaps at the sequence termini (beginning and end) incur no penalty.

**Mathematical Formulation:**

**Initialization:**
```
F(0,0) = 0
F(i,0) = 0  for all i  (free gaps at start of seq1)
F(0,j) = 0  for all j  (free gaps at start of seq2)
```

**Recurrence Relation:**
```
F(i,j) = max {
    F(i-1, j-1) + S(xi, yj)    # Match/Mismatch
    F(i-1, j) + gap_penalty     # Deletion
    F(i, j-1) + gap_penalty     # Insertion
}

where:
  S(xi, yj) = PAM250[xi][yj]  (substitution score)
  gap_penalty = -9
```

**Termination:**
```
Optimal score = max { F(m,j) for all j, F(i,n) for all i }
(Maximum score in last row or last column)
```

**Why Semi-Global?**
- Useful when one sequence might be a subsequence of another
- Common in:
  - Database searching (query vs database sequence)
  - Sequence assembly (finding overlaps)
  - Gene finding (exon vs genome)
  - Primer design (primer vs template)

### Scoring Matrices

#### PAM (Point Accepted Mutation)

**Origin:**
- Developed by Margaret Dayhoff (1978)
- Based on evolutionary model
- Uses closely related protein families (>85% identity)

**Construction:**
```
1. Align closely related sequences (phylogenetically close)
2. Count observed amino acid substitutions
3. Calculate mutation probability matrix (PAM1)
4. Extrapolate to higher PAM numbers (PAM2, PAM3, ...)
```

**PAM Number Interpretation:**
- **PAM1**: 1 accepted point mutation per 100 amino acids
- **PAM250**: 250 accepted mutations per 100 amino acids
- Higher PAM → More evolutionary distance → Lower similarity

**Usage:**
- PAM1 - PAM100: Closely related sequences
- PAM120 - PAM160: Medium distance
- **PAM250**: Distantly related sequences (most divergent)

**Applications:**
- Phylogenetic tree construction
- Evolutionary analysis
- DNA-related studies

#### BLOSUM (BLOck SUbstitution Matrix)

**Origin:**
- Developed by Henikoff & Henikoff (1992)
- Based on conserved blocks in multiple alignments
- Uses sequences from BLOCKS database

**Construction:**
```
1. Identify conserved blocks in protein families
2. Count amino acid substitutions within blocks
3. Calculate substitution frequencies
4. Compute log-odds scores
```

**BLOSUM Number Interpretation:**
- **BLOSUM45**: Sequences with 45% identity (more divergent)
- **BLOSUM62**: Sequences with 62% identity (default, general use)
- **BLOSUM80**: Sequences with 80% identity (closely related)
- Higher BLOSUM → Higher similarity → Less evolutionary distance

**Usage:**
- BLOSUM45: Distantly related sequences
- **BLOSUM62**: General purpose (most common)
- BLOSUM80: Closely related sequences

**Applications:**
- Database searching (BLAST uses BLOSUM62)
- Finding conserved regions
- Protein domain identification

#### PAM vs BLOSUM Comparison

| Feature | PAM | BLOSUM |
|---------|-----|--------|
| **Base Data** | Closely related sequences | Divergent sequences |
| **Method** | Evolutionary extrapolation | Direct observation |
| **Calculation** | Pairwise comparisons | Multiple sequence alignment blocks |
| **Number Meaning** | Mutations per 100 residues | % identity threshold |
| **Number Direction** | Higher = more divergent | Higher = less divergent |
| **Best for Global** | ✅ Yes (whole sequence) | ❌ Less suitable |
| **Best for Local** | ❌ Less suitable | ✅ Yes (conserved blocks) |
| **Common Versions** | PAM250 (distant) | BLOSUM62 (general) |

**Rule of Thumb:**
```
Closely Related Sequences:
  PAM1 - PAM100  ≈  BLOSUM80 - BLOSUM90

Medium Distance:
  PAM120 - PAM160  ≈  BLOSUM62

Distantly Related:
  PAM200 - PAM250  ≈  BLOSUM45 - BLOSUM50
```

**This Project Uses:**
- **PAM250**: Suitable for distantly related protein sequences
- Gap penalty: -9 (linear model)

---

## 💻 Implementation

### Algorithm Details

**Dynamic Programming Approach:**

```python
# 1. INITIALIZATION
score_matrix[0][0] = 0
# Free gaps at boundaries (semi-global)
for i in range(1, len(seq1)+1):
    score_matrix[i][0] = 0
for j in range(1, len(seq2)+1):
    score_matrix[0][j] = 0

# 2. FILL MATRIX
for i in range(1, len(seq1)+1):
    for j in range(1, len(seq2)+1):
        diagonal = score_matrix[i-1][j-1] + PAM250[seq1[i-1]][seq2[j-1]]
        horizontal = score_matrix[i][j-1] + GAP_PENALTY
        vertical = score_matrix[i-1][j] + GAP_PENALTY
        
        score_matrix[i][j] = max(diagonal, horizontal, vertical)
        
        # Track all optimal paths
        if score_matrix[i][j] == diagonal:
            direction_matrix[i][j] += 'd'  # diagonal
        if score_matrix[i][j] == horizontal:
            direction_matrix[i][j] += 'h'  # horizontal
        if score_matrix[i][j] == vertical:
            direction_matrix[i][j] += 'v'  # vertical

# 3. FIND OPTIMAL SCORE
# Semi-global: check last row and last column
optimal_score = max(
    max(score_matrix[i][len(seq2)] for i in range(len(seq1)+1)),
    max(score_matrix[len(seq1)][j] for j in range(len(seq2)+1))
)

# 4. TRACEBACK (ALL optimal paths)
def traceback(i, j, aligned_seq1, aligned_seq2):
    if i == 0 or j == 0:
        # Add remaining gaps for semi-global
        while i > 0:
            aligned_seq1 += seq1[i-1]
            aligned_seq2 += '-'
            i -= 1
        while j > 0:
            aligned_seq1 += '-'
            aligned_seq2 += seq2[j-1]
            j -= 1
        return
    
    for direction in direction_matrix[i][j]:
        if direction == 'd':
            traceback(i-1, j-1, seq1[i-1] + aligned_seq1, seq2[j-1] + aligned_seq2)
        elif direction == 'h':
            traceback(i, j-1, '-' + aligned_seq1, seq2[j-1] + aligned_seq2)
        elif direction == 'v':
            traceback(i-1, j, seq1[i-1] + aligned_seq1, '-' + aligned_seq2)
```

### Key Features

**1. PAM250 Substitution Matrix**
- 20×20 matrix for all amino acid pairs
- Values range from -8 (W↔C, highly unfavorable) to 17 (W↔W, identical)
- Reflects evolutionary likelihood of substitutions

**2. Complete Traceback**
- Finds **ALL** optimal alignments, not just one
- Uses direction matrix to track multiple paths
- Essential when multiple equally good alignments exist

**3. Lexicographic Sorting**
- Outputs alignments in sorted order for reproducibility
- Consistent with judge/grading system requirements

**4. Free End Gaps**
- No penalty for gaps at sequence termini
- Implemented via initialization (all zeros)
- Optimal score searched in last row/column

### Code Structure

```
semi_global_alignment.py
├── PAM250 = {...}                    # 20×20 substitution matrix
├── GAP_PENALTY = -9                  # Linear gap penalty
│
├── init_and_fill_matrix()            # Initialize & fill DP matrices
│   ├── Creates score_matrix
│   ├── Creates direction_matrix
│   └── Returns filled matrices
│
├── calculate_score(arr, i, j)       # Calculate cell score
│   ├── Computes diagonal, horizontal, vertical
│   ├── Tracks direction(s) for optimal path
│   └── Returns (score, direction)
│
├── find_total_score_locations()     # Find optimal endpoints
│   ├── Searches last row & column
│   └── Returns list of (i, j) coordinates
│
├── trace_back(strings, x, y)        # Recursive traceback
│   ├── Follows direction matrix
│   ├── Handles multiple paths
│   ├── Adds end gaps when needed
│   └── Builds all optimal alignments
│
├── semi_global_alignment()          # Main alignment function
│   ├── Starts from optimal endpoints
│   ├── Adds trailing end gaps
│   └── Calls trace_back for each endpoint
│
└── print_output(score, seq)         # Format and display results
    ├── Prints optimal score
    ├── Sorts alignments lexicographically
    └── Outputs sequence pairs
```

---

## 📦 Installation

### Prerequisites
```bash
Python 3.x (no external libraries required)
```

### Setup

**Clone Repository:**
```bash
git clone https://github.com/zamirmehdi/Bioinformatics-Course.git
cd "Bioinformatics-Course/2- Pairwise Sequence Alignment/src"
```

**No Dependencies:**
The implementation uses only Python standard library - no pip installations needed!

---

## 🚀 Usage

### Input Format

**Two protein sequences** (uppercase letters, amino acid codes):
```
SEQUENCE1
SEQUENCE2
```

**Valid Amino Acids:**
```
A C D E F G H I K L M N P Q R S T V W Y
(20 standard amino acids)
```

**Constraints:**
- Maximum sequence length: 100 characters
- Time limit: 1.5 seconds
- Memory limit: 100 MB

### Output Format

```
SCORE
ALIGNMENT1_SEQ1
ALIGNMENT1_SEQ2
ALIGNMENT2_SEQ1
ALIGNMENT2_SEQ2
...
```

**Output Rules:**
1. First line: Optimal alignment score
2. Following lines: All optimal alignments (sorted)
3. Use `-` for gaps
4. Alignments sorted lexicographically

### Examples

#### Example 1: Simple Alignment

**Input:**
```
HEAGAWGHE
PAWHEA
```

**Output:**
```
20
HEAGAWGHE-
---PAW-HEA
```

**Explanation:**
- Score: 20 (from PAM250 matches minus gap penalties)
- One optimal alignment found
- Leading gaps (---) and trailing gap (-) are free (semi-global)

#### Example 2: Multiple Optimal Alignments

**Input:**
```
AAAAA
AA
```

**Output:**
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

**Explanation:**
- Score: 4 (2 × PAM250[A][A] = 2 × 2 = 4)
- Four equally optimal alignments
- AA can align with any two A's in AAAAA
- All sorted lexicographically

**Alignment Breakdown:**
```
Position:  1 2 3 4 5
Sequence:  A A A A A

Option 1:  - - - A A  (align with positions 4,5)
Option 2:  - - A A -  (align with positions 3,4)
Option 3:  - A A - -  (align with positions 2,3)
Option 4:  A A - - -  (align with positions 1,2)
```

### Running the Program

**Method 1: Standard Input**
```bash
python semi_global_alignment.py
# Enter sequences when prompted
HEAGAWGHE
PAWHEA
```

**Method 2: Input Redirection**
```bash
echo -e "HEAGAWGHE\nPAWHEA" | python semi_global_alignment.py
```

**Method 3: File Input**
```bash
python semi_global_alignment.py < input.txt
```

---

## 📝 Theoretical Assignment

The theoretical component covers manual calculations and conceptual understanding:

### Question 1: Global Alignment Manual Calculation

**Task:** 
a) Compute global alignment for `ACCTAGA` and `ACTGG` using Needleman-Wunsch
b) Verify using dot matrix method
c) Compare both results

**Scoring:**
- Match: +1
- Mismatch: -1  
- Gap penalty: -2

**Learning Objectives:**
- Understand DP matrix filling
- Practice traceback
- Compare algorithmic approaches

### Question 2: Scoring Parameter Identification

**Given:** Partially filled alignment matrix

**Tasks:**
a) Determine gap penalty, match score, mismatch score from matrix
b) Complete the matrix
c) Find optimal alignment and score
d) Identify all alternative optimal alignments

**Skills:**
- Reverse-engineering alignment parameters
- Pattern recognition in DP matrices
- Multiple traceback paths

### Question 3: Local Alignment (Smith-Waterman)

**Sequences:**
```
TATAGC
GTTATC
```

**Scoring:**
- Match: +2
- Mismatch: -1
- Gap penalty: -2

**Calculate:**
- Local alignment matrix
- Optimal local alignment
- Score

**Key Difference from Global:**
- Can start/end anywhere
- Minimum cell value = 0 (no negative scores)
- Traceback stops at 0

### Question 4: Real Protein Comparison

**Given 3 protein sequences:**
- Sequence A: 214 amino acids
- Sequence B: 181 amino acids (similar to A)
- Sequence C: 389 amino acids (divergent from A)

**Tasks:**

**a) Choose Appropriate Alignment:**
- A vs B: Which alignment type? Why?
- A vs C: Which alignment type? Why?

**Answer Summary:**
```
A vs B: Global Alignment
  Reason: Similar lengths, high similarity (81.3% identity)
          Close evolutionary distance
          
A vs C: Local Alignment  
  Reason: Very different lengths, low similarity (25.9% identity)
          Large evolutionary distance
          Only partial regions similar
```

**b) Use Online Tools:**
- EMBOSS Needle (global): https://www.ebi.ac.uk/Tools/psa/emboss_needle/
- EMBOSS Water (local): https://www.ebi.ac.uk/Tools/psa/emboss_water/

**Results Analysis:**
```
Global (A,B):  Score = 872.0   Identity = 81.3%   Gaps = 15.4%
Local  (A,B):  Score = 876.0   Identity = 97.2%   Gaps = 0.0%

Global (A,C):  Score = 410.5   Identity = 25.9%   Gaps = 45.4%
Local  (A,C):  Score = 419.0   Identity = 52.3%   Gaps = 1.7%
```

**Conclusions:**
- Local alignment scores higher (focuses on conserved regions)
- A-B: Scores similar (sequences are uniformly similar)
- A-C: Local much better (only partial similarity exists)
- Gap percentage indicates evolutionary distance

### Question 5: PAM vs BLOSUM Research

**Task:** Compare scoring matrix methodologies

**Required Topics:**
- Construction methods
- Versions and numbering systems
- Use cases and applications
- Advantages/disadvantages

**See:** [Scoring Matrices](#scoring-matrices) section above for complete comparison

---

## 📊 Performance Analysis

### Time Complexity

**Initialization:** O(m + n)
```python
for i in range(m+1):      # O(m)
    score_matrix[i][0] = 0
for j in range(n+1):      # O(n)
    score_matrix[0][j] = 0
```

**Matrix Filling:** O(m × n)
```python
for i in range(1, m+1):       # O(m)
    for j in range(1, n+1):   # O(n)
        calculate_score()      # O(1)
```

**Finding Optimal Endpoints:** O(m + n)
```python
for i in last_row:     # O(n)
    check_score()
for j in last_col:     # O(m)
    check_score()
```

**Traceback:** O(k × (m + n))
```
k = number of optimal alignments
Each traceback: O(m + n) path length
```

**Total:** O(m × n + k × (m + n))
- Dominated by matrix filling: **O(m × n)**
- Practical: Fast for sequences < 1000 residues

### Space Complexity

**Score Matrix:** O(m × n)
```python
score_matrix = [[0] * (n+1) for _ in range(m+1)]
```

**Direction Matrix:** O(m × n)
```python
direction_matrix = [[''] * (n+1) for _ in range(m+1)]
```

**Alignment Storage:** O(k × (m + n))
```
k alignments, each length ≤ m+n
```

**Total:** **O(m × n)** space
- Can be reduced to O(n) for score-only computation
- Full traceback requires O(m × n) for direction matrix

### Optimization Opportunities

**1. Space Optimization (Score Only):**
```python
# If only score needed, not alignments:
# Use only 2 rows instead of full matrix
current_row = [0] * (n+1)
previous_row = [0] * (n+1)
```

**2. Early Termination:**
```python
# If sequences very dissimilar, can detect early
# Not implemented in current version
```

**3. Banded Alignment:**
```python
# If sequences known to be similar
# Only compute diagonal band ±k
# Reduces to O(k × min(m,n))
```

**4. Parallel Processing:**
```python
# Anti-diagonal cells are independent
# Can compute in parallel
# Speedup up to min(m,n) processors
```

---

## 🗂️ Project Structure

```
2- Pairwise Sequence Alignment/
├── src/
│   └── semi_global_alignment.py    # Main implementation
├── docs/
│   ├── Programming Instruction.pdf  # Programming assignment (Persian)
│   └── Theoretical/
│       ├── Instruction.pdf          # Theoretical questions (Persian)
│       └── Report.pdf               # Completed answers (Persian)
└── README.md                        # This documentation
```

---

## 🎓 Key Concepts Covered

### Algorithms
- ✅ Dynamic Programming fundamentals
- ✅ Needleman-Wunsch algorithm (global)
- ✅ Smith-Waterman algorithm (local)
- ✅ Semi-global alignment variant
- ✅ Complete traceback for multiple solutions

### Bioinformatics
- ✅ Sequence alignment significance
- ✅ Substitution matrices (PAM, BLOSUM)
- ✅ Gap penalty models (linear vs affine)
- ✅ Evolutionary distance estimation
- ✅ Protein sequence comparison

### Data Structures
- ✅ 2D arrays for DP matrices
- ✅ Backtracking with direction pointers
- ✅ String manipulation and formatting
- ✅ Recursive solution enumeration

### Computational Biology
- ✅ Sequence homology detection
- ✅ Conserved region identification
- ✅ Database searching foundations
- ✅ Phylogenetic relationships

---

## 🎯 Learning Outcomes

After completing this project, students can:

### Programming Skills
✅ Implement dynamic programming algorithms efficiently  
✅ Handle multi-dimensional arrays and matrices  
✅ Perform recursive backtracking with multiple paths  
✅ Optimize code for time and space constraints  
✅ Debug complex algorithmic problems

### Bioinformatics Knowledge
✅ Understand alignment types and their applications  
✅ Choose appropriate scoring matrices (PAM vs BLOSUM)  
✅ Interpret alignment scores and identity percentages  
✅ Recognize biological significance of alignments  
✅ Apply computational methods to real proteins

### Analytical Thinking
✅ Analyze algorithm complexity (time/space)  
✅ Compare different algorithmic approaches  
✅ Validate implementation against known results  
✅ Optimize solutions for specific constraints  
✅ Interpret biological data computationally

---

## 📈 Results Summary

*(From submitted theoretical assignment)*

### Manual Calculations

**Question 1:** Global alignment of `ACCTAGA` vs `ACTGG`
- Successfully completed using Needleman-Wunsch
- Verified with dot matrix visualization
- Both methods yielded consistent results

**Question 2:** Parameter identification and matrix completion
- Correctly identified: match=+1, mismatch=-1, gap=-2
- Completed full DP matrix
- Found multiple optimal alignments

**Question 3:** Local alignment of `TATAGC` vs `GTTATC`
- Applied Smith-Waterman algorithm
- Identified optimal local alignment region
- Calculated score and verified

### Protein Comparison Analysis

**Sequences A & B (Similar):**
- Length: 214 vs 181 amino acids
- Global alignment recommended ✓
- Results:
  - Identity: 81.3%
  - Score: 872.0 (global) vs 876.0 (local)
  - Interpretation: Closely related, likely orthologs

**Sequences A & C (Divergent):**
- Length: 214 vs 389 amino acids  
- Local alignment recommended ✓
- Results:
  - Identity: 25.9% (global) vs 52.3% (local)
  - Score: 410.5 (global) vs 419.0 (local)
  - Interpretation: Distantly related, partial homology

**Key Insights:**
- Local alignment improves score for divergent sequences
- A-B pair: Uniform similarity (global ≈ local)
- A-C pair: Conserved domains only (local >> global)
- Gap percentage correlates with evolutionary distance

### PAM vs BLOSUM Understanding

Comprehensive comparison provided in theoretical report:
- Construction methodologies
- Numbering systems (opposite meanings!)
- Application guidelines
- Appropriate use cases

---

## ℹ️ Project Information

**Assignment:** Pairwise Sequence Alignment  
**Author:** Amirmehdi Zarrinnezhad  
**Course:** Bioinformatics  
**University:** Amirkabir University of Technology (Tehran Polytechnic) - Fall 2022  
**Language:** Python 3.x (Programming), Persian (Theoretical)  
**Submission Platform:** Quera  
**GitHub Link:** [2- Pairwise Sequence Alignment](https://github.com/zamirmehdi/Bioinformatics-Course/tree/main/2-%20Pairwise%20Sequence%20Alignment)

<div align="center">

**Part of Bioinformatics Course Projects**

[1: Basic Biology](../1-%20Basic%20biology) | [2: Sequence Alignment](.) | [3: MSA & DB Search](../3-%20Multiple%20Sequence%20Alignment%20-%20DB%20Search) | [4: Profile HMM](../4-%20Profile%20-%20Hidden%20Markov%20model) | [5: Phylogenetic Trees](../5-%20Phylogenetic%20Trees) | [Final: Virus Classification](../Virus%20Classification%20(Final%20Project))

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
