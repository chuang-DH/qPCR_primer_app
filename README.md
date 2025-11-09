🧬 qPCR Primer Designer GUI
Overview

This project is a qPCR Primer Design Tool built with Python and Tkinter.
It allows users to input a DNA sequence (up to 15,000 bp) and automatically generates several primer pairs suitable for qPCR experiments.

The program ensures that selected primers meet common qPCR conditions, including:

GC content around 50% (acceptable range: 40–60%)

Melting temperature (Tm) close to 60 °C

Amplicon length suitable for qPCR detection

It also provides a clean and intuitive Graphical User Interface (GUI) for quick analysis and result export.

🚀 Features

End-to-end automation
From raw DNA input to optimal primer pair selection, the algorithm automatically performs sequence cleaning, GC and Tm evaluation, structural filtering (self-dimer, 3′ GC-rich, hairpin), and amplicon pairing—no manual parameter tuning required.

Heuristic + Thermodynamic Scoring System
Each primer and primer pair is ranked through a multi-parameter scoring system that considers GC deviation, Tm deviation, structural penalties, and amplicon size uniformity. This hybrid strategy ensures that selected primers are both stable and experimentally realistic.

Adaptive Melting Temperature (Tm) Calculation
The system dynamically selects the most appropriate Tm formula according to primer length and composition—using simplified nearest-neighbor approximations for short oligos and thermodynamic corrections for longer sequences—thereby preventing distortion in temperature estimation and improving thermodynamic accuracy.

Structural Awareness
Unlike basic GC/Tm filters, this algorithm performs secondary structure checks (self-complementarity, 3′-GC richness, and hairpin loops), greatly reducing false-positive primer designs that would fail in wet-lab qPCR validation.

High Efficiency for Long Sequences
Capable of handling up to 20,000 bp efficiently through pruning and position-indexed pairing, allowing fast primer design even on consumer-grade hardware.

Customizable Design Parameters and Cleavage Information
Users can freely adjust target primer length ranges, melting temperature (Tm), GC content percentage, Tm tolerance, and amplicon length. The system automatically calculates and displays the cleavage positions and resulting nucleotide lengths for each candidate primer, enabling precise control over experimental design details.

Interactive Web Interface (Streamlit)
The Streamlit GUI allows users to adjust constraints (e.g., primer length, GC%, Tm, amplicon range) and visualize results instantly, with real-time feedback and ranked tables of candidate primer pairs. This makes the tool accessible to both bioinformatics researchers and experimental biologists.

Lightweight and Deployable Anywhere
Written entirely in Python with no external dependencies beyond Streamlit, it can be deployed locally or on a web server (e.g., Streamlit Cloud, Hugging Face Spaces) for immediate use by lab members.

🧩 Technologies Used

The core algorithm for primer design integrates multi-parameter optimization to identify reliable primer candidates suitable for qPCR applications. It operates in two main stages: primer candidate generation and primer pair evaluation.

In the first stage, the function find_primer_candidates() scans the input DNA sequence across all possible windows within the specified primer length range (min_len to max_len). For each subsequence, it computes the GC content and melting temperature (Tm) using established thermodynamic rules. Candidates are filtered based on user-defined constraints for GC percentage (gc_min, gc_max) and temperature tolerance (tm_target, tm_tol). Each surviving primer undergoes additional quality checks, including self-complementarity, 3′-end GC richness, and potential hairpin formation, which are flagged as warnings. The resulting list of candidates is annotated with their positions, lengths, GC%, Tm, and any detected structural issues.

In the second stage, design_qpcr_primers() orchestrates the design of primer pairs. Both forward and reverse primer candidates are generated independently and pruned based on a composite score function that prioritizes thermodynamic stability and uniformity. This scoring penalizes deviations from the ideal GC content (50%) and target Tm, while incorporating additional penalties for structural warnings such as self-dimers and hairpins. The algorithm then iteratively pairs each forward primer with downstream reverse primers, calculating the corresponding amplicon length and reverse-complement sequence. Only pairs producing amplicons within the acceptable range (amp_min–amp_max) are retained. Each pair is further ranked by an aggregate performance score that considers GC balance, Tm deviation, structural warnings, and amplicon size uniformity.

By combining heuristic filtering with parametric scoring, this implementation efficiently searches the sequence space for the most stable and specific primer pairs. The final ranked list of primer pairs (top_n) represents optimal candidates for qPCR assays, balancing accuracy, thermodynamic compatibility, and structural robustness.

⚡technical workflow

┌───────────────────────────┐
│ 1. User Input DNA Sequence │
│    (Streamlit GUI)        │
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│ 2. Sequence Sanitization  │
│ sanitize_seq(seq)         │
│ - Remove invalid bases    │
│ - Convert to uppercase    │
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│ 3. Thermodynamic Parameter│
│    Calculation & Filtering│
│ - GC content              │
│   (calculate_gc_content)  │
│ - Melting temperature (Tm)│
│   (calculate_tm, adaptive │
│    formula based on length│
│    and composition)       │
│ - Structural checks:      │
│   - Self-complementarity  │
│   - 3′ GC-rich tail       │
│   - Hairpin formation     │
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│ 4. Primer Candidate       │
│    Generation (Forward &  │
│    Reverse)               │
│ find_primer_candidates()  │
│ - Record position, length,│
│   GC%, Tm, warnings       │
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│ 5. Primer Pair Evaluation │
│ design_qpcr_primers()    │
│ - Compute amplicon length │
│   (Amp_len)               │
│ - Compute composite score │
│   (GC, Tm, structural     │
│    penalties, Amp_len)    │
│ - Sort & select Top N     │
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│ 6. Results Display & Output│
│    (Streamlit GUI)        │
│ - Show Top N candidate    │
│   primer pairs            │
│ - Include F/R sequence,   │
│   position, Tm, GC%,      │
│   Amp_len, warnings       │
└───────────────────────────┘

🧠 How It Works

Paste or upload your DNA sequence.

Click "Design Primers".

The program will scan the sequence and output several optimal primer pairs.

You can view GC%, Tm, and other details directly in the interface.

🌐 Online Demo

Try it on Streamlit:
👉 https://qpcrprimerapp-smypfaljr4wtha97gtyub5.streamlit.app/

🧑‍💻 Author

Developed by Chuang, Dong-Hua
mail: benalu85853@gmail.com
If you find this useful, please ⭐ star the repository and share feedback!
