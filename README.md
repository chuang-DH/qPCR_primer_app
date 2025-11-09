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

Accepts long input sequences (≤ 15,000 bp)

Automatically finds primer pairs matching qPCR standards

Displays primer details (Tm, GC%, product size, etc.)

Simple GUI built with Tkinter

Easy to use for both beginners and researchers

🧩 Technologies Used

The core algorithm for primer design integrates multi-parameter optimization to identify reliable primer candidates suitable for qPCR applications. It operates in two main stages: primer candidate generation and primer pair evaluation.

In the first stage, the function find_primer_candidates() scans the input DNA sequence across all possible windows within the specified primer length range (min_len to max_len). For each subsequence, it computes the GC content and melting temperature (Tm) using established thermodynamic rules. Candidates are filtered based on user-defined constraints for GC percentage (gc_min, gc_max) and temperature tolerance (tm_target, tm_tol). Each surviving primer undergoes additional quality checks, including self-complementarity, 3′-end GC richness, and potential hairpin formation, which are flagged as warnings. The resulting list of candidates is annotated with their positions, lengths, GC%, Tm, and any detected structural issues.

In the second stage, design_qpcr_primers() orchestrates the design of primer pairs. Both forward and reverse primer candidates are generated independently and pruned based on a composite score function that prioritizes thermodynamic stability and uniformity. This scoring penalizes deviations from the ideal GC content (50%) and target Tm, while incorporating additional penalties for structural warnings such as self-dimers and hairpins. The algorithm then iteratively pairs each forward primer with downstream reverse primers, calculating the corresponding amplicon length and reverse-complement sequence. Only pairs producing amplicons within the acceptable range (amp_min–amp_max) are retained. Each pair is further ranked by an aggregate performance score that considers GC balance, Tm deviation, structural warnings, and amplicon size uniformity.

By combining heuristic filtering with parametric scoring, this implementation efficiently searches the sequence space for the most stable and specific primer pairs. The final ranked list of primer pairs (top_n) represents optimal candidates for qPCR assays, balancing accuracy, thermodynamic compatibility, and structural robustness.

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
