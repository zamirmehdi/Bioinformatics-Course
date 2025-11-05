# Project 5 - Phylogenetic Trees

A comprehensive theoretical exploration of **phylogenetic tree construction** methods, including distance-based algorithms (UPGMA, Neighbor-Joining), parsimony-based approaches, and maximum likelihood methods for inferring evolutionary relationships.

[![Theory](https://img.shields.io/badge/Type-Theoretical-purple.svg)](#)
[![Course](https://img.shields.io/badge/Course-Bioinformatics-blue.svg)](#)
[![University](https://img.shields.io/badge/University-AUT-red.svg)](https://aut.ac.ir/en)

<details> <summary><h2>📚 Table of Contents</h2></summary>

- [Overview](#-overview)
- [Phylogenetic Trees Fundamentals](#-phylogenetic-trees-fundamentals)
- [Assignment Questions](#-assignment-questions)
  - [Question 1: Regular Expressions for Motifs](#question-1-regular-expressions-for-motifs)
  - [Question 2: UPGMA vs Neighbor-Joining](#question-2-upgma-vs-neighbor-joining)
  - [Question 3: Parsimony Analysis](#question-3-parsimony-analysis)
  - [Question 4: Exhaustive Search vs Branch-and-Bound](#question-4-exhaustive-search-vs-branch-and-bound)
  - [Question 5: Maximum Likelihood Trees](#question-5-maximum-likelihood-trees)
- [Key Concepts Covered](#-key-concepts-covered)
- [Learning Outcomes](#-learning-outcomes)
- [Project Structure](#%EF%B8%8F-project-structure)
- [Project Information](#ℹ%EF%B8%8F-project-information)
- [Contact](#-contact)

</details>

---

## 📋 Overview

This theoretical project explores **phylogenetic tree reconstruction**, a fundamental problem in computational biology for inferring evolutionary relationships among species, genes, or sequences. Phylogenetic trees reveal:

- **Evolutionary History**: How organisms diverged from common ancestors
- **Species Classification**: Taxonomic relationships and biodiversity
- **Gene Function**: Protein family evolution and functional conservation
- **Disease Evolution**: Tracking pathogen spread and mutations
- **Conservation Biology**: Identifying endangered lineages

**Key Topics Covered:**
- Regular expressions for protein motif description
- Distance-based methods (UPGMA, Neighbor-Joining)
- Parsimony-based tree construction
- Exhaustive search vs branch-and-bound algorithms
- Maximum likelihood phylogenetic inference

**Assignment Type:** Theoretical Analysis (No Programming)  
**Submission:** Written answers with manual calculations and tree diagrams

---

## 🌳 Phylogenetic Trees Fundamentals

### What is a Phylogenetic Tree?

**Definition**: A branching diagram showing evolutionary relationships among biological entities based on their characteristics.

**Tree Components:**
```
         Root (common ancestor)
           |
      Internal Node (ancestor)
        /     \
    Branch   Branch
      /         \
   Leaf A      Leaf B
(taxa/species) (taxa/species)
```

**Tree Types:**

**1. Rooted Tree**
```
    Root
    /  \
   A    \
       / \
      B   C

Interpretation: Shows direction of evolution
Time flows from root → leaves
```

**2. Unrooted Tree**
```
    A
    |
    •---B
    |
    C

Interpretation: Shows relationships, not direction
No explicit common ancestor position
```

### Tree Construction Philosophies

**Three Main Approaches:**

**1. Distance-Based Methods**
- Input: Pairwise distances between sequences
- Algorithm: Cluster based on similarity
- Examples: UPGMA, Neighbor-Joining (NJ)
- Speed: Fast (polynomial time)
- Accuracy: Good for closely related species

**2. Character-Based Methods (Parsimony)**
- Input: Aligned sequences
- Principle: Minimize evolutionary changes (mutations)
- Algorithm: Find tree with fewest substitutions
- Speed: Slow (NP-hard, exponential)
- Accuracy: Good for slowly evolving sequences

**3. Probabilistic Methods (Maximum Likelihood)**
- Input: Aligned sequences + evolutionary model
- Principle: Find tree with highest probability
- Algorithm: Evaluate likelihood of all trees
- Speed: Very slow (computationally intensive)
- Accuracy: Best, especially with good model

**Comparison Table:**

| Method | Input | Principle | Speed | Accuracy | Use Case |
|--------|-------|-----------|-------|----------|----------|
| **UPGMA** | Distances | Average linkage | Fast | Moderate | Molecular clock data |
| **NJ** | Distances | Minimum evolution | Fast | Good | General purpose |
| **Parsimony** | Sequences | Minimize changes | Slow | Good | Conserved sequences |
| **ML** | Sequences | Maximize probability | Very slow | Best | Publication quality |

---

## 📝 Assignment Questions

### Question 1: Regular Expressions for Motifs

**Regular expressions** describe conserved patterns in biological sequences using concise syntax.

#### Part A: Creating Regular Expression

**Task**: Write a regex for the protein motif:
```
AYGTTSKK
AYPTTSIK
AVHTTSIK
AYMTTSIK
AVZTTSIK
```

**Analysis**:

**Position-by-position comparison:**
```
Position:  1 2 3 4 5 6 7 8
Seq1:      A Y G T T S K K
Seq2:      A Y P T T S I K
Seq3:      A V H T T S I K
Seq4:      A Y M T T S I K
Seq5:      A V Z T T S I K

Pattern:   A ? ? T T S ? K
```

**Conservation Analysis:**
- **Position 1**: Always A (conserved)
- **Position 2**: Y or V (two options)
- **Position 3**: G, P, H, M, or Z (any amino acid = X)
- **Position 4-5**: Always TT (conserved)
- **Position 6**: Always S (conserved)
- **Position 7**: I or K (two options)
- **Position 8**: Always K (conserved)

**Regular Expression:**
```
A-[YV]-X-T(2)-S-[IK]-K
```

**Syntax Explanation:**
- `A`: Exact match for amino acid A
- `-`: Separator (not part of sequence)
- `[YV]`: Either Y or V
- `X`: Any amino acid (wildcard)
- `T(2)` or `TT`: Two consecutive T's
- `S`: Exact match for S
- `[IK]`: Either I or K
- `K`: Exact match for K

---

#### Part B: Regular Expression Matching

**Given Regex**: `M-[TG]-X-{M}-A(2)-P-[YPC]`

**Syntax Meaning:**
- `M`: Must be M
- `[TG]`: Must be T or G
- `X`: Any character
- `{M}`: Anything EXCEPT M
- `A(2)`: Two consecutive A's (AA)
- `P`: Must be P
- `[YPC]`: Must be Y, P, or C

**Test Sequences:**

**1. MMTGAAPP**
```
Sequence:  M  M  T  G  A  A  P  P
Regex:     M [TG] X {M} A  A  P [YPC]
           ✓  ✗

Result: NO MATCH
Reason: Position 2 is M, but regex requires T or G
```

**2. MTTTAAPC**
```
Sequence:  M  T  T  T  A  A  P  C
Regex:     M [TG] X {M} A  A  P [YPC]
           ✓  ✓  ✓  ✓  ✓  ✓  ✓  ✓

Result: EXACT MATCH ✓
Explanation:
- M matches M
- T matches [TG]
- T matches X (anything)
- T matches {M} (not M)
- AA matches A(2)
- P matches P
- C matches [YPC]
```

**3. MGTMAAPP**
```
Sequence:  M  G  T  M  A  A  P  P
Regex:     M [TG] X {M} A  A  P [YPC]
           ✓  ✓  ✓  ✗

Result: NO MATCH
Reason: Position 4 is M, but {M} means "anything except M"
```

**4. MTGAAPPY**
```
Sequence:  M  T  G  A  A  P  P  Y
Regex:     M [TG] X {M} A  A  P [YPC]
           ✓  ✓  ✓  ✗

Alternative alignment:
Sequence:  M  T  G  A  A  P  P  Y
Regex:     M [TG] X {M} A(2) P [YPC]
           ✓  ✓  ✓  ✓  ?

Result: NO MATCH
Reason: Position 6 is P, but regex expects A (second A of A(2))
```

**Summary**: Only **MTTTAAPC** matches the regular expression.

---

### Question 2: UPGMA vs Neighbor-Joining

**Given Tree with Branch Lengths:**
```
        ┌─13─ A
     ┌──┤
     │  └──4─ B
  ───┤
     │  ┌──4─ C
     └──┤
        └─10─ D
```

**Task**: Extract distance matrix, build trees using UPGMA and NJ, compare results.

---

#### Part A: Distance Matrix Extraction

**Calculate pairwise distances** (sum branch lengths):

```
d(A,B) = 13 + 4 = 17
d(A,C) = 13 + 2 + 2 + 4 = 21
d(A,D) = 13 + 2 + 2 + 10 = 27
d(B,C) = 4 + 2 + 2 + 4 = 12
d(B,D) = 4 + 2 + 2 + 10 = 18
d(C,D) = 4 + 10 = 14
```

**Distance Matrix:**
```
     A   B   C   D
A    -  17  21  27
B   17   -  12  18
C   21  12   -  14
D   27  18  14   -
```

---

#### Part B: UPGMA Tree Construction

**UPGMA** (Unweighted Pair Group Method with Arithmetic Mean):
- Assumes **molecular clock** (constant evolution rate)
- Creates **ultrametric tree** (all leaves at same distance from root)
- Uses **average linkage clustering**

**Algorithm:**

**Round 1: Find minimum distance**
```
Minimum: d(B,C) = 12

Join B and C into cluster BC
Branch lengths: B─6─•─6─C (each gets half the distance)
```

**Update distances (average of B and C)**:
```
d(BC, A) = [d(B,A) + d(C,A)] / 2 = (17 + 21) / 2 = 19
d(BC, D) = [d(B,D) + d(C,D)] / 2 = (18 + 14) / 2 = 16
```

**New Matrix:**
```
      BC   A   D
BC     -  19  16
A     19   -  27
D     16  27   -
```

**Round 2: Find minimum distance**
```
Minimum: d(BC, D) = 16

Join BC and D into cluster BCD
Branch length from internal node: 16/2 = 8
```

**Update distances**:
```
d(BCD, A) = [d(BC,A) + d(D,A)] / 2 = (19 + 27) / 2 = 23
```

**Round 3: Join BCD and A**
```
Final tree with A joining at root
```

**UPGMA Tree:**
```
         ┌─── A
      ───┤
         │  ┌─ B
         └──┤
            │  C
            └─ D

All leaves at same height (ultrametric)
```

**Characteristics**:
- All taxa at same level from root
- Assumes equal evolution rates
- Branch lengths NOT shown explicitly

---

#### Part C: Neighbor-Joining Tree Construction

**Neighbor-Joining (NJ)**:
- Does NOT assume molecular clock
- Creates **additive tree** (unequal branch lengths)
- Minimizes total tree length

**Algorithm:**

**Step 1: Calculate r values** (sum of distances to all other taxa)
```
r_A = d(A,B) + d(A,C) + d(A,D) = 17 + 21 + 27 = 65
r_B = d(B,A) + d(B,C) + d(B,D) = 17 + 12 + 18 = 47
r_C = d(C,A) + d(C,B) + d(C,D) = 21 + 12 + 14 = 47
r_D = d(D,A) + d(D,B) + d(D,C) = 27 + 18 + 14 = 59

Normalized: r'_i = r_i / (n-2) = r_i / 2
r'_A = 32.5, r'_B = 23.5, r'_C = 23.5, r'_D = 29.5
```

**Step 2: Calculate d' (corrected distances)**
```
d'_ij = d_ij - (r'_i + r'_j)

d'_AB = 17 - (32.5 + 23.5) = -39  ← Minimum
d'_AC = 21 - (32.5 + 23.5) = -35
d'_AD = 27 - (32.5 + 29.5) = -35
d'_BC = 12 - (23.5 + 23.5) = -35
d'_BD = 18 - (23.5 + 29.5) = -35
d'_CD = 14 - (23.5 + 29.5) = -39  ← Also minimum
```

**Step 3: Join pair with minimum d'** (choose A-B)
```
Calculate branch lengths to new node U:
d_AU = [d_AB + (r'_A - r'_B)] / 2 = [17 + (32.5 - 23.5)] / 2 = 13
d_BU = d_AB - d_AU = 17 - 13 = 4
```

**Step 4: Update distances to new cluster U(AB)**
```
d_CU = [d_AC - d_AU + d_BC - d_BU] / 2 = [(21-13) + (12-4)] / 2 = 8
d_DU = [d_AD - d_AU + d_BD - d_BU] / 2 = [(27-13) + (18-4)] / 2 = 14
```

**New Matrix:**
```
      U(AB)  C   D
U(AB)   -    8  14
C       8    -  14
D      14   14   -
```

**Step 5: Repeat** (join U(AB) and C)
```
d'_U(AB)C = 8 - (11 + 11) = -14  ← Minimum

d_CU' = [d_CU + (r'_C - r'_U)] / 2 = [8 + 0] / 2 = 4
d_UU' = 4
d_DU' = [(14-4) + (14-4)] / 2 = 10
```

**Final NJ Tree:**
```
    A─13─┐
         ├─4─┐
    B─4──┘   ├─4─┐
             │   │
    C─4──────┘   └─10─ D

Branch lengths explicitly shown
```

---

#### Part D: Comparison with Original Tree

**Original Tree Pattern:**
```
A─13─┐
     ├─2─┐
B─4──┘   ├─2─┐
         │   │
C─4──────┘   └─10─ D
```

**Comparison Table:**

| Feature | Original | UPGMA | NJ |
|---------|----------|-------|-----|
| **First join** | A-B or C-D | B-C | A-B |
| **Second join** | - | BC-D | AB-C |
| **Final join** | - | BCD-A | ABC-D |
| **Branch lengths** | Explicit | Hidden | Explicit |
| **Ultrametric** | No | Yes | No |
| **Similarity to original** | - | Low | Higher |

**Key Differences:**

**1. Clustering Order**
- **UPGMA**: B-C → BC-D → BCD-A
- **NJ**: A-B → AB-C → ABC-D
- **Original**: Has AB pairing like NJ

**2. Molecular Clock Assumption**
- **UPGMA**: Forces equal distances from root (ultrametric)
  - All leaves at same level
  - Cannot represent varying evolution rates
- **NJ**: Allows unequal branch lengths
  - Matches original tree property
  - Reflects real evolution (some lineages faster)

**3. Distance Calculation**
- **UPGMA**: Simple averaging of distances
- **NJ**: Uses corrected distances considering all taxa
  - More comprehensive
  - Better handles rate variation

**4. Topological Accuracy**
- **NJ closer to original**: Both have AB grouping
- **UPGMA different**: Groups BC first
- NJ better preserves evolutionary relationships

**Why NJ Performs Better:**
- Original tree has unequal branch lengths (not ultrametric)
- UPGMA's molecular clock assumption violated
- NJ's flexibility handles rate heterogeneity
- NJ considers global tree structure (all pairwise distances)

---

### Question 3: Parsimony Analysis

**Parsimony Principle**: Evolution follows the path of fewest changes

**Task**: Find minimum number of mutations in rooted tree

**Given Tree:**
```
      Root (red dot)
      /    \
    ...    ...
   /  \    /  \
  T   T   A    C
  
... and other leaves
```

**Parsimony Algorithm: Fitch's Algorithm**

#### Step 1: Post-Order Traversal (Leaves → Root)

**At each internal node**: Calculate possible nucleotide sets

**Rules:**
- If children share nucleotides: **Intersection** (no mutation)
- If children don't overlap: **Union** (mutation required)

**Example Calculation:**

**Leaves**: T, T, A, G, C

**Internal Node 1** (subtree with T, T):
```
Left child: {T}
Right child: {T}
Intersection: {T} ∩ {T} = {T}
Mutations: 0 (no change needed)
```

**Internal Node 2** (subtree with A, G):
```
Left child: {A}
Right child: {G}
Intersection: {A} ∩ {G} = ∅ (empty)
Union: {A, G}
Mutations: 1 (one change required)
```

**Continue upward...**

#### Step 2: Root → Leaves (Traceback)

**Select nucleotide at each node** minimizing changes

**From Report:**
```
Minimum mutations: 5
Locations marked on tree branches
```

**Root Independence:**
- Changing root position doesn't change minimum mutations
- Why? Tree topology and relative relationships unchanged
- Only starting point of traversal changes, not total changes needed
- Parsimony score is an intrinsic property of the tree topology

**Example:**
```
Original root:     New root:
    A                  B
   / \                / \
  B   C              A   C
  
Same tree, just reoriented
Same minimum mutations (3 in this example)
```

---

### Question 4: Exhaustive Search vs Branch-and-Bound

**Problem**: Finding optimal phylogenetic tree is NP-hard

**Why Hard?**
- Number of possible trees grows super-exponentially
- For n taxa: (2n-3)!! unrooted trees
- Examples:
  - 5 taxa: 105 trees
  - 10 taxa: 34,459,425 trees
  - 20 taxa: ~10²⁰ trees (infeasible!)

---

#### Exhaustive Search

**Method**: Evaluate ALL possible trees

**Algorithm:**
```
1. Start with 3-taxon tree (any 3 species)
2. Add 4th taxon in all possible positions
   - Score each resulting tree
3. Add 5th taxon to each tree from step 2
   - Score all new trees
4. Continue until all taxa added
5. Return tree with best score
```

**Advantages:**
✅ **Guaranteed optimal**: Finds best tree (global optimum)  
✅ **Complete search**: No tree missed  
✅ **Simple to understand**: Straightforward logic  

**Disadvantages:**
❌ **Extremely slow**: O((2n-3)!!) complexity  
❌ **Limited scalability**: ~10 taxa maximum  
❌ **Wasteful**: Evaluates many poor trees  

**Time Complexity:**
```
n=5:  105 trees           → seconds
n=10: 34,459,425 trees    → hours
n=15: 2×10¹³ trees        → years
n=20: 8×10²⁰ trees        → universe lifetime
```

---

#### Branch-and-Bound

**Method**: Exhaustive search with **pruning** of suboptimal branches

**Key Idea**: If partial tree already worse than best known complete tree, don't continue expanding it

**Algorithm:**
```
1. Build initial tree (e.g., using UPGMA/NJ)
   - Calculate its score → Initial upper bound
   
2. Start building trees systematically (like exhaustive)
   
3. While building:
   - Calculate lower bound for current partial tree
   - If lower_bound > upper_bound:
     → PRUNE (don't expand this branch)
     → Skip all descendant trees
   - Else:
     → Continue expanding
     
4. Update upper bound when better tree found
   
5. Return best tree after exploring all unpruned branches
```

**Example:**
```
Best known score: 12 (upper bound)

Partial tree A: minimum possible score = 8
  → Continue (might improve to <12)

Partial tree B: minimum possible score = 15
  → PRUNE (cannot be better than 12)
  → Skip millions of descendant trees!
```

**Advantages:**
✅ **Guaranteed optimal**: Still finds best tree  
✅ **Faster than exhaustive**: Prunes many branches  
✅ **Scales better**: ~20 taxa feasible  

**Disadvantages:**
❌ **Still exponential**: Just better constant factor  
❌ **Complex implementation**: Requires good bounds  
❌ **Bound-dependent**: Performance varies with data  

**Comparison:**

| Feature | Exhaustive | Branch-and-Bound |
|---------|-----------|------------------|
| **Optimality** | Guaranteed | Guaranteed |
| **Completeness** | Full search | Full search (with pruning) |
| **Speed** | Slowest | Faster (via pruning) |
| **Scalability** | ~10 taxa | ~20 taxa |
| **Implementation** | Simple | Complex |
| **Pruning** | None | Aggressive |

**Pruning Efficiency Example:**
```
15 taxa problem:
Exhaustive: Evaluate 2×10¹² trees
Branch-and-bound: Evaluate ~10⁸ trees (99.995% pruned!)
```

**Both methods find same optimal tree**, but branch-and-bound gets there faster by avoiding provably suboptimal paths.

---

### Question 5: Maximum Likelihood Trees

**Maximum Likelihood (ML)**: Find tree that maximizes probability of observed data

**Given:**
- Equal prior probabilities: P(A) = P(C) = P(G) = P(T) = 0.25
- Substitution matrix with transition probabilities

**Substitution Matrix** (example values):
```
     A    C    G    T
A   0.55 0.20 0.15 0.10
C   0.05 0.70 0.10 0.15
G   0.15 0.10 0.60 0.15
T   0.25 0.05 0.10 0.60
```

**Task**: Calculate probability of three different tree topologies

---

#### Tree (A): T → T, T, A, T, C, G

**Tree Structure:**
```
        T (root)
       /|\
      / | \
     T  •  •
       /|  |\
      T A  T •
             |\
             C G
```

**Probability Calculation:**

**Method 1: With Prior Probabilities**
```
P(Tree) = P(root=T) × 
          P(T→T) × P(T→T) × P(T→A) × 
          P(T→T) × P(T→C) × P(T→G)

= 0.25 × 0.6 × 0.6 × 0.25 × 0.6 × 0.05 × 0.1
= 0.0000675
```

**Method 2: Without Prior** (for comparison only)
```
P(Tree) = P(T→T) × P(T→T) × P(T→A) × 
          P(T→T) × P(T→C) × P(T→G)

= 0.6 × 0.6 × 0.25 × 0.6 × 0.05 × 0.1
= 0.00027
```

---

#### Tree (B): G → A, A, A, A, T, C, G

**Tree Structure:**
```
        G (root)
       /|\
      / | \
     •  •  G
    /|  |\
   A A  A •
         |\
         T C
```

**Probability Calculation:**
```
P(Tree) = P(root=G) × 
          P(G→A) × P(A→A) × P(A→A) × 
          P(A→T) × P(A→C) × P(G→G)

= 0.25 × 0.15 × 0.55 × 0.55 × 0.1 × 0.2 × 0.6
= 0.00013612
```

**Without Prior:**
```
= 0.15 × 0.55 × 0.55 × 0.1 × 0.2 × 0.6
= 0.0005445  ← Highest probability!
```

---

#### Tree (C): C → C, A, A, A, T, C, G

**Tree Structure:**
```
        C (root)
       /|\
      / | \
     C  •  •
       /|  |\
      A A  • C
          /|
         A T
```

**Probability Calculation:**
```
P(Tree) = P(root=C) × 
          P(C→C) × P(C→A) × P(A→A) × 
          P(A→T) × P(C→C) × P(C→G)

= 0.25 × 0.7 × 0.05 × 0.55 × 0.1 × 0.7 × 0.15
= 0.00005053
```

**Without Prior:**
```
= 0.7 × 0.05 × 0.55 × 0.1 × 0.7 × 0.15
= 0.0002
```

---

#### Maximum Likelihood Tree Selection

**With Priors:**
```
Tree A: 0.0000675
Tree B: 0.00013612  ← Maximum!
Tree C: 0.00005053

ML Tree: Tree B (highest probability)
```

**Without Priors:**
```
Tree A: 0.00027
Tree B: 0.0005445  ← Maximum!
Tree C: 0.0002

ML Tree: Tree B (still highest)
```

**Interpretation:**
- **Tree B has highest likelihood** under given model
- Would be selected as best phylogenetic hypothesis
- Prior probabilities affect absolute values but not ranking (if equal)
- More likely tree = better fit to evolutionary model

**Why Tree B Wins:**
- Fewer low-probability transitions
- More self-transitions (A→A, G→G are likely)
- Better match to substitution matrix patterns

---

## 🎓 Key Concepts Covered

### Sequence Patterns
- ✅ Regular expression syntax for biological motifs
- ✅ Pattern matching and validation
- ✅ Conserved position identification
- ✅ Wildcard and exclusion operators

### Distance-Based Methods
- ✅ UPGMA algorithm and molecular clock
- ✅ Neighbor-Joining and additive trees
- ✅ Ultrametric vs non-ultrametric trees
- ✅ Distance matrix manipulation
- ✅ Comparison of tree construction methods

### Parsimony Analysis
- ✅ Minimum evolution principle
- ✅ Fitch's algorithm for parsimony scoring
- ✅ Rooted tree analysis
- ✅ Mutation counting and optimization
- ✅ Root position independence

### Search Algorithms
- ✅ Exhaustive search guarantees
- ✅ Branch-and-bound pruning
- ✅ Computational complexity (NP-hardness)
- ✅ Scalability limitations
- ✅ Optimality trade-offs

### Maximum Likelihood
- ✅ Probabilistic phylogenetic inference
- ✅ Substitution models and matrices
- ✅ Tree probability calculation
- ✅ Model-based evolution
- ✅ Likelihood comparison

---

## 🎯 Learning Outcomes

After completing this project, students can:

### Theoretical Understanding
✅ Explain different phylogenetic tree construction philosophies  
✅ Compare UPGMA, Neighbor-Joining, Parsimony, and ML methods  
✅ Understand molecular clock assumption and its violations  
✅ Recognize when each method is appropriate  
✅ Interpret phylogenetic trees correctly  

### Algorithmic Skills
✅ Manually construct trees using UPGMA algorithm  
✅ Apply Neighbor-Joining with corrected distances  
✅ Perform parsimony analysis with Fitch's algorithm  
✅ Calculate maximum likelihood for given trees  
✅ Understand branch-and-bound pruning strategy  

### Analytical Thinking
✅ Compare tree topologies and identify differences  
✅ Evaluate computational complexity of algorithms  
✅ Assess trade-offs between speed and accuracy  
✅ Recognize limitations of different methods  
✅ Select appropriate method for given data  

### Practical Skills
✅ Write regular expressions for sequence motifs  
✅ Extract and manipulate distance matrices  
✅ Draw phylogenetic trees with correct topology  
✅ Calculate evolutionary probabilities  
✅ Interpret biological meaning of trees  

---

## 🗂️ Project Structure

```
5- Phylogenetic Trees/
├── Instruction.pdf          # Assignment questions (Persian)
├── Report.pdf              # Completed solutions (Persian)
└── README.md               # This documentation
```

**File Descriptions:**

**`Instruction.pdf`**:
- 5 comprehensive theory questions
- Regular expression exercises
