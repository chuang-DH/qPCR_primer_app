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

Python 3.x

Tkinter (for GUI)

Biopython / primer3-py (for sequence and primer analysis)

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
