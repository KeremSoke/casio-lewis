# Lewis Structure + VSEPR for CASIO Calculators

This is a MicroPython program (`lewis.py`) for the CASIO fx-9750GIII and fx-9860GIII graphing calculators. It is designed to help chemistry students by automatically generating plausible Lewis structures, predicting molecular geometry with VSEPR theory, and identifying resonance forms for a given chemical formula.

The program is entirely self-contained and has been specifically tailored to run on the limited, non-standard MicroPython environment of these calculators.

## Features

-   **Lewis Structure Generation:** Automatically calculates the total number of valence electrons and determines the arrangement of bonds and lone pairs.
-   **Formal Charge Calculation:** Computes the formal charge for each atom in the structure to help determine stability.
-   **Resonance Structures:** Identifies and displays multiple valid resonance forms when they exist, based on minimizing formal charges.
-   **VSEPR Theory Prediction:** Based on the most optimal structure, it predicts:
    -   AXE Notation (e.g., AX₂E₂)
    -   Molecular Shape (e.g., Bent, Tetrahedral)
    -   Ideal Bond Angles (e.g., 109.5°, <120°)
    -   Central Atom Hybridization (e.g., sp³, sp²)
-   **Supports:** Neutral molecules, anions (negative ions), and cations (positive ions).

## Compatibility

This script is specifically designed for the MicroPython 1.9.4 environment found on the following calculators:
-   **CASIO fx-9750GIII**
-   **CASIO fx-9860GIII**

Due to the unique and limited nature of this MicroPython version (e.g., missing `enumerate`, `sys`, and having non-standard method behavior), this script is **not guaranteed** to work on other Python or MicroPython platforms without modification.

## Installation

1.  **Connect Your Calculator:** Use a USB cable to connect your CASIO calculator to your computer.
2.  **Enter USB Mode:** On the calculator, navigate to the main menu, select the "LINK" or "USB" icon, and choose "USB Mass Storage". The calculator will appear as a removable drive on your computer.
3.  **Create the File:** On your computer, create a new plain text file and name it `lewis.py`.
4.  **Copy the Code:** Copy the entire Python code for the program into this `lewis.py` file.
5.  **Transfer the File:** Drag and drop the `lewis.py` file from your computer into the root directory of the calculator's drive.
6.  **Eject Safely:** Safely eject the calculator's drive from your computer.

## How to Use

1.  On your calculator, go to the main menu and select the **Python** icon.
2.  You will see a list of `.py` files. Find and select `lewis.py`.
3.  Press the **[F1]** key to **Run** the script.
4.  The screen will display a prompt: `Formula:`.
5.  Type the chemical formula you want to analyze and press the **[EXE]** key.
6.  The program will output the most optimal Lewis structure, any resonance forms, and a complete VSEPR prediction.

## Input Format Rules

To ensure the formula is parsed correctly, follow these rules:

-   **Case Sensitive:** Element symbols must be capitalized correctly (e.g., `CO2`, not `co2`).
-   **Atom Counts:** Write numbers immediately after the element symbol (e.g., `H2O`, `CH4`).
-   **Ionic Charge:**
    -   For charges of `+1` or `-1`, simply add `+` or `-` at the end (e.g., `NH4+`, `OH-`).
    -   For charges greater than 1, add the sign and then the number (e.g., `SO4-2`, `CO3-2`).
    -   **Do not use spaces** in the formula.

---

## Examples

### Example 1: Water (`H2O`)

A simple, neutral molecule.

**Input:**
```
Formula: H2O
```

**Output:**
```
=========================
--- Most Optimal Structure ---

[Structure]
  O0 - H0
  O0 - H1

[Lone Pairs]
  O0: 4e (2 pairs)

[Formal Charges]
  All charges are zero.
=========================

--- VSEPR Prediction ---
Central atom: O
Bonded atoms (X): 2
Lone pairs (E): 2
Notation: AX2E2
Shape: Bent
Bond angle(s): <109.5
Hybridization: sp3
```

### Example 2: Carbonate Ion (`CO3-2`)

An ion with double bonds and resonance.

**Input:**
```
Formula: CO3-2
```
**Output:**
```
=========================
--- Most Optimal Structure ---

[Structure]
  C0 = O0
  C0 - O1
  C0 - O2

[Lone Pairs]
  O0: 4e (2 pairs)
  O1: 6e (3 pairs)
  O2: 6e (3 pairs)

[Formal Charges]
  O1: -1
  O2: -1
=========================

--- VSEPR Prediction ---
Central atom: C
Bonded atoms (X): 3
Lone pairs (E): 0
Notation: AX3
Shape: Trigonal Planar
Bond angle(s): 120
Hybridization: sp2

This molecule exhibits resonance.

=========================
--- Resonance Form 2 ---

[Structure]
  C0 - O0
  C0 = O1
  C0 - O2
...
(Additional resonance forms will be printed)
```

### Example 3: Sulfur Hexafluoride (`SF6`)

A molecule with an expanded octet.

**Input:**
```
Formula: SF6
```
**Output:**
```
=========================
--- Most Optimal Structure ---

[Structure]
  S0 - F0
  S0 - F1
  S0 - F2
  S0 - F3
  S0 - F4
  S0 - F5

[Lone Pairs]
  F0: 6e (3 pairs)
  F1: 6e (3 pairs)
  ... (and so on for all F atoms)

[Formal Charges]
  All charges are zero.
=========================

--- VSEPR Prediction ---
Central atom: S
Bonded atoms (X): 6
Lone pairs (E): 0
Notation: AX6
Shape: Octahedral
Bond angle(s): 90
Hybridization: sp3d2
```

## Limitations

-   **Text-Based:** The program provides a textual description of the structure, not a visual drawing.
-   **Idealized Angles:** The VSEPR bond angles are ideal (e.g., "109.5" or "<109.5") and do not account for distortions from lone pairs or different terminal atoms.
-   **Algorithm Heuristics:** The "optimal" structure is determined by minimizing formal charges, which is a powerful guideline but may not perfectly represent every complex molecule.
-   **No Radicals:** The program is designed for molecules with an even number of valence electrons and does not handle radicals (odd number of electrons).

## License

This project is licensed under the terms of the GNU General Public License v3.0.

| [![GNU Head][gnu-logo]][gpl-v3-link] | This work is free software; you can redistribute it and/or modify it under the terms of the **[GNU General Public License v3.0][gpl-v3-link]** as published by the Free Software Foundation. <br><br> This work is distributed in the hope that it will be useful, but **WITHOUT ANY WARRANTY**; without even the implied warranty of **MERCHANTABILITY** or **FITNESS FOR A PARTICULAR PURPOSE**. See the [GNU General Public License][gpl-v3-link] for more details. |
| :--- | :--- |

## Disclaimer of Liability

The GPLv3 license explicitly states that this software is provided "AS IS" with **ABSOLUTELY NO WARRANTY**. The user assumes all risks associated with its use. The author is not liable for any damages, including but not limited to inaccuracies in calculations, data loss, or any other issues arising from the use of this program. Please verify all calculations before entering them into official logbooks or records.

<!-- Reference links used in the table above for cleaner Markdown -->
[gnu-logo]: https://upload.wikimedia.org/wikipedia/commons/thumb/2/22/Heckert_GNU_white.svg/100px-Heckert_GNU_white.svg.png
[gpl-v3-link]: https://www.gnu.org/licenses/gpl-3.0.html

