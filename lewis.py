# =============================================================================
# lewis.py
#
# A Lewis Structure and VSEPR Geometry calculator for the CASIO fx-9750GIII
# and compatible calculators running MicroPython 1.9.4.
#
# This program determines plausible Lewis structures for a given chemical
# formula, including ions, resonance structures, and molecules with
# expanded octets. It then uses the most stable structure to predict
# molecular geometry based on VSEPR theory.
#
# This script is specifically tailored to the limitations of the CASIO
# MicroPython environment (e.g., no 're' or 'sys' modules, broken 'raise'
# command, and missing 'enumerate' function).
#
# =============================================================================

# --- VSEPR Data Table ---
# A dictionary that maps the number of bonded atoms (X) and lone pairs (E)
# on a central atom to its VSEPR geometry.
# Format: (X, E): ("AXE Notation", "Shape", "Ideal Angle", "Hybridization")
VSEPR = {
    (1,0): ("AX","Linear","180","s"),
    (2,0): ("AX2","Linear","180","sp"),
    (2,1): ("AX2E","Bent","<120","sp2"),
    (3,0): ("AX3","Trigonal Planar","120","sp2"),
    (3,1): ("AX3E","Trigonal Pyramidal","<109.5","sp3"),
    (2,2): ("AX2E2","Bent","<109.5","sp3"),
    (4,0): ("AX4","Tetrahedral","109.5","sp3"),
    (4,1): ("AX4E","Seesaw","<120 & <90","sp3d"),
    (3,2): ("AX3E2","T-shaped","<90","sp3d"),
    (2,3): ("AX2E3","Linear","180","sp3d"),
    (5,0): ("AX5","Trigonal Bipyramidal","120 & 90","sp3d2"),
    (5,1): ("AX5E","Square Pyramidal","<90","sp3d2"),
    (4,2): ("AX4E2","Square Planar","90","sp3d2"),
    (6,0): ("AX6","Octahedral","90","sp3d2"),
}

def show_vsepr_info(structure, central_atom):
    """
    Calculates and prints the VSEPR geometry prediction based on the
    number of bonded atoms and lone pairs from the final Lewis structure.
    """
    X = len(structure['bonds'])
    E = structure['lone_pairs'].get(central_atom, 0) // 2
    
    print("\n--VSEPR Prediction--")
    print("Central atom: " + str(central_atom[:-1]))
    print("Bonded atoms (X): " + str(X))
    print("Lone pairs (E): " + str(E))
    
    info = VSEPR.get((X, E))
    if info:
        code, shape, angle, hyb = info
        print("Notation: " + str(code))
        print("Shape: " + str(shape))
        print("Bond angle(s): " + str(angle))
        print("Hybridization: " + str(hyb))
    else:
        print("No VSEPR entry for X=" + str(X) + ", E=" + str(E))

def combinations(iterable, r):
    """
    A custom implementation of the 'combinations' function from Python's
    itertools library, which is not available in this MicroPython environment.
    """
    pool = tuple(iterable)
    n = len(pool)
    if r > n: return
    indices = list(range(r))
    yield tuple(pool[i] for i in indices)
    while True:
        # Loop backwards through the indices to find the next combination
        for i in range(r - 1, -1, -1):
            if indices[i] != i + n - r: break
        else: return
        indices[i] += 1
        for j in range(i + 1, r):
            indices[j] = indices[j - 1] + 1
        yield tuple(pool[i] for i in indices)

# =============================================================================
# Main Lewis Structure Class
# =============================================================================
class LewisStructureCreator:
    # --- Chemical Data ---
    VALENCE_ELECTRONS = {'H':1,'B':3,'C':4,'N':5,'O':6,'F':7,'Be':2,'Li':1,'Na':1,'Mg':2,'Al':3,'Si':4,'P':5,'S':6,'Cl':7,'K':1,'Ca':2,'Ga':3,'Ge':4,'As':5,'Se':6,'Br':7,'I':7,'Xe':8,'Kr':8}
    ELECTRONEGATIVITY = {'H':2.2,'B':2.04,'C':2.55,'N':3.04,'O':3.44,'F':3.98,'Be':1.57,'Li':0.98,'Na':0.93,'Mg':1.31,'Al':1.61,'Si':1.9,'P':2.19,'S':2.58,'Cl':3.16,'K':0.82,'Ca':1.0,'Ga':1.81,'Ge':2.01,'As':2.18,'Se':2.55,'Br':2.96,'I':2.66,'Xe':2.6,'Kr':3.0}

    def __init__(self, formula_string):
        """
        Initializes the creator, parses the formula, and sets up
        the necessary variables for structure generation.
        """
        self.is_valid = False
        # The `is_valid` flag is used for error handling, as the `raise`
        # command is broken on the target MicroPython platform.
        self.atom_counts, self.charge = self._parse_formula(formula_string)
        if self.atom_counts is None: return # Stop if parsing failed
        
        self.total_valence_electrons = self._calculate_total_valence_electrons()
        self.all_valid_structures = []
        self.central_atom = ""
        self.terminal_atoms = []
        self.is_valid = True

    def _parse_formula(self, formula):
        """
        Parses a chemical formula string (e.g., "SO4-2") into a dictionary
        of atom counts and an integer charge, without using the 're' module.
        """
        atom_counts = {}
        charge = 0
        temp = formula
        
        # Manually parse the charge from the end of the string.
        # This handles formats like "NO3-", "NH4+", and "SO4-2".
        # The check `temp.endswith('-') or temp.endswith('+')` is a critical
        # bug fix, as this version of MicroPython does not support tuples
        # in the .endswith() method.
        if len(temp) >= 1 and (temp.endswith('-') or temp.endswith('+')):
            charge = -1 if temp.endswith('-') else 1
            temp = temp[:-1]
        elif len(temp) > 1 and temp[-1].isdigit() and (temp[-2] == '+' or temp[-2] == '-'):
            charge_val = int(temp[-1])
            charge = charge_val * (-1 if temp[-2] == '-' else 1)
            temp = temp[:-2]
        
        # Manually parse the atoms and their counts.
        i = 0
        while i < len(temp):
            s = temp[i]; i += 1
            # Handle two-letter element symbols (e.g., "Cl", "Mg").
            if i < len(temp) and temp[i].islower(): s += temp[i]; i += 1
            
            if s not in self.VALENCE_ELECTRONS:
                print("Error: Element '" + str(s) + "' not supported.")
                return None, None
            
            num = ""
            while i < len(temp) and temp[i].isdigit(): num += temp[i]; i += 1
            
            count = int(num) if num else 1
            atom_counts[s] = atom_counts.get(s, 0) + count
            
        if not atom_counts:
            print("Error: No elements found.")
            return None, None
            
        return atom_counts, charge

    def _calculate_total_valence_electrons(self):
        """Calculates the total number of valence electrons for the molecule."""
        total = 0
        for s, c in self.atom_counts.items(): total += self.VALENCE_ELECTRONS[s] * c
        return total - self.charge

    def _find_central_atom(self):
        """
        Determines the central atom based on a hierarchy of rules:
        1. If there's only one type of atom, it's central.
        2. Hydrogen is never central.
        3. The atom that appears only once is usually central.
        4. Otherwise, the least electronegative atom is central.
        """
        keys = list(self.atom_counts.keys())
        if len(keys) == 1: return keys[0]
        
        candidates = {s: c for s, c in self.atom_counts.items() if s != 'H'}
        if not candidates: return min(keys, key=lambda s: self.ELECTRONEGATIVITY[s])
        
        singles = [s for s, c in candidates.items() if c == 1]
        if len(singles) == 1: return singles[0]
        
        return min(list(candidates.keys()), key=lambda s: self.ELECTRONEGATIVITY[s])

    def _generate_structures(self):
        """
        The main engine for generating the initial, single-bonded skeleton
        of the molecule and distributing electrons.
        """
        if self.total_valence_electrons % 2 != 0: return # Handle radicals
        
        c_symbol = self._find_central_atom()
        
        # Create a list of all terminal atoms.
        terms = []
        for s, cnt in self.atom_counts.items():
            num = cnt - 1 if s == c_symbol else cnt
            for _ in range(num): terms.append(s)
        
        # Name all atoms uniquely to prevent bugs (e.g., in O3).
        # Central atom is always 'Symbol0'.
        self.central_atom = c_symbol + "0"
        self.terminal_atoms = []
        temp_counts = {}
        for s in terms:
            # Terminal atoms are 'Symbol1', 'Symbol2', etc. This prevents a
            # name collision with the central atom.
            idx = temp_counts.get(s, 1) 
            self.terminal_atoms.append(s + str(idx))
            temp_counts[s] = idx + 1
        
        # Build the single-bonded skeleton.
        nb = len(self.terminal_atoms)
        e_in_b = nb * 2
        if e_in_b > self.total_valence_electrons: return # Not enough electrons
        
        bonds = {a: 1 for a in self.terminal_atoms}
        lp_electrons = self.total_valence_electrons - e_in_b
        
        # Distribute lone pairs, filling terminal atoms first.
        lone_pairs = {a: 0 for a in [self.central_atom] + self.terminal_atoms}
        for t in self.terminal_atoms:
            sym = t[:-1] if t[-1].isdigit() else t
            need = 0 if sym == 'H' else 6 # H needs 0, others need 6 to complete octet.
            take = min(lp_electrons, need)
            lone_pairs[t] += take; lp_electrons -= take
        
        # Any remaining electrons go to the central atom.
        lone_pairs[self.central_atom] += lp_electrons
        
        # Start the recursive process to find all plausible structures.
        self._satisfy_central_octet(bonds, lone_pairs)

    def _satisfy_central_octet(self, bonds, lone_pairs):
        """
        A recursive function that explores all plausible bonding patterns by
        promoting lone pairs from terminal atoms to form multiple bonds.
        It includes intelligent pruning to prevent freezes on the calculator.
        """
        # --- Pre-calculation and Pruning ---
        fc = {}
        all_atoms = [self.central_atom] + self.terminal_atoms
        for atom in all_atoms:
            sym = atom[:-1] if atom[-1].isdigit() else atom
            val = self.VALENCE_ELECTRONS[sym]
            nb = bonds.get(atom, 0) if atom != self.central_atom else sum(bonds.values())
            lp = lone_pairs.get(atom, 0)
            fc[atom] = val - lp - nb

        # Pruning Rule 1: Stop if a terminal halogen has a positive charge.
        # This is a chemically impossible situation and prevents freezes on molecules like SCl4.
        for atom in self.terminal_atoms:
            if (atom[:-1] if atom[-1].isdigit() else atom) in ['F', 'Cl', 'Br', 'I']:
                if fc.get(atom, 0) > 0:
                    return

        # If the structure is plausible, store it.
        self._store_if_valid(bonds, lone_pairs)

        # --- Intelligent "Stop" Rule ---
        central_symbol = self.central_atom[:-1]
        bonds_electrons = sum(bonds.values()) * 2
        central_electrons = bonds_electrons + lone_pairs.get(self.central_atom, 0)
        
        # Check if all terminals are neutral.
        all_terminals_neutral = True
        for atom in self.terminal_atoms:
            if fc.get(atom, 0) != 0:
                all_terminals_neutral = False
                break
        
        # If the central atom has a stable octet and good formal charges, there's no
        # reason to expand the octet. Pruning this branch stops freezes on PF6-, etc.
        if central_electrons == 8 and fc.get(self.central_atom, 0) <= 0 and all_terminals_neutral:
            return

        # --- Recursive Step ---
        # If the central atom is from Period 2, it cannot exceed an octet.
        is_period_2 = central_symbol in ['B', 'C', 'N', 'O', 'F']
        if is_period_2 and central_electrons >= 8:
            return

        # For each terminal atom that can donate a lone pair, create a new
        # structure with a multiple bond and analyze it recursively.
        for terminal_atom in self.terminal_atoms:
            if lone_pairs.get(terminal_atom, 0) >= 2 and bonds.get(terminal_atom, 1) < 3:
                new_bonds = bonds.copy()
                new_lone_pairs = lone_pairs.copy()
                new_bonds[terminal_atom] += 1
                new_lone_pairs[terminal_atom] -= 2
                
                self._satisfy_central_octet(new_bonds, new_lone_pairs)

    def _store_if_valid(self, bonds, lone_pairs):
        """Calculates formal charges and stores the complete structure if it's unique."""
        structure = {'bonds': bonds, 'lone_pairs': lone_pairs}
        fc = {}
        all_atoms = [self.central_atom] + self.terminal_atoms
        for atom in all_atoms:
            sym = atom[:-1] if atom[-1].isdigit() else atom
            val = self.VALENCE_ELECTRONS[sym]
            nb = bonds.get(atom, 0) if atom != self.central_atom else sum(bonds.values())
            lp = lone_pairs.get(atom, 0)
            charge = val - lp - nb
            fc[atom] = charge
        structure['formal_charges'] = fc
        
        # Add the structure only if it's not already in the list.
        if structure not in self.all_valid_structures:
            self.all_valid_structures.append(structure)

    def get_optimal_structures(self):
        """
        Finds the best structure(s) from all generated possibilities based on
        a sophisticated scoring system that minimizes formal charges and places
        them on the most appropriate atoms.
        """
        self._generate_structures()
        if not self.all_valid_structures: return None, []

        best_score = (100, 100) # (primary_score, penalty_score)
        optimal_structures = []

        for s in self.all_valid_structures:
            # 1. Primary Score: The sum of absolute formal charges (lower is better).
            primary_score = sum([abs(c) for c in s['formal_charges'].values()])

            # 2. Penalty Score (Tie-breaker): Penalize poor charge placement.
            penalty_score = 0
            for atom, charge in s['formal_charges'].items():
                symbol = atom[:-1] if atom[-1].isdigit() else atom
                en = self.ELECTRONEGATIVITY[symbol]
                
                # Penalize positive charges on electronegative atoms.
                if charge > 0:
                    penalty_score += charge * en 
                # Penalize negative charges on less electronegative atoms.
                elif charge < 0:
                    penalty_score += abs(charge) * (4.0 - en)
            
            current_score = (primary_score, penalty_score)

            # Compare to the best score found so far.
            if current_score[0] < best_score[0]:
                best_score = current_score
                optimal_structures = [s] # Found a new best, start the list over.
            elif current_score[0] == best_score[0]:
                if current_score[1] < best_score[1]:
                    best_score = current_score
                    optimal_structures = [s] # Same primary score, but better tie-breaker.
                elif current_score[1] == best_score[1]:
                    optimal_structures.append(s) # Identical scores, this is a resonance form.

        if not optimal_structures:
            return None, []

        return optimal_structures[0], optimal_structures[1:]

    def format_structure(self, s, title):
        """Formats a single structure's data into a readable string."""
        out = ["--- " + str(title) + " ---", "\n[Structure]"]
        bond_symbols = {1:'-', 2:'=', 3:'~'} # Using ~ for triple bond
        for atom, order in s['bonds'].items():
            out.append("  " + self.central_atom + " " + bond_symbols.get(order, '?') + " " + atom)
        
        out.append("\n[Lone Pairs]")
        # Sorting keys before printing is necessary for consistent output.
        keys_lp = sorted(list(s['lone_pairs'].keys()))
        for k in keys_lp:
            lp = s['lone_pairs'][k]
            if lp > 0: out.append("  " + k + ": " + str(lp) + "e (" + str(lp//2) + " pairs)")
        
        out.append("\n[Formal Charges]")
        anyc = False
        keys_fc = sorted(list(s['formal_charges'].keys()))
        for k in keys_fc:
            ch = s['formal_charges'][k]
            if ch != 0:
                anyc = True
                sign = "+" if ch > 0 else ""
                out.append("  " + k + ": " + sign + str(ch))
        if not anyc: out.append("  All charges are zero.")
        
        return "\n".join(out)

# =============================================================================
# Main Execution Block
# =============================================================================
def main():
    """Main function to handle user input and run the program."""
    print("Lewis Structure + VSEPR")
    try:
        formula = input("Formula: ")
        creator = LewisStructureCreator(formula)
        if not creator.is_valid: return

        opt, resonances = creator.get_optimal_structures()
        
        if opt:
            print("\n" + "="*21)
            print(creator.format_structure(opt, "Most Optimal Structure"))
            print("="*21)
            show_vsepr_info(opt, creator.central_atom)

            if resonances:
                print("\nThis molecule exhibits resonance.")
                i = 1
                for res_structure in resonances:
                    print("\n" + "="*25)
                    print(creator.format_structure(res_structure, "Resonance Form " + str(i + 1)))
                    print("="*25)
                    i += 1
        else:
            print("\nCould not find a valid structure.")
    except Exception as e:
        # A general error handler for any unexpected issues during development.
        print("\nAn unexpected error occurred: " + str(e))

main()
