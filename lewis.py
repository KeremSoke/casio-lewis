# MicroPython 1.9.4 - Full-Featured Version for CASIO fx-9750GIII
# Includes fixes for expanded octets and unique atom naming.

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
    pool = tuple(iterable)
    n = len(pool)
    if r > n: return
    indices = list(range(r))
    yield tuple(pool[i] for i in indices)
    while True:
        for i in range(r - 1, -1, -1):
            if indices[i] != i + n - r: break
        else: return
        indices[i] += 1
        for j in range(i + 1, r):
            indices[j] = indices[j - 1] + 1
        yield tuple(pool[i] for i in indices)

class LewisStructureCreator:
    VALENCE_ELECTRONS = {'H':1,'B':3,'C':4,'N':5,'O':6,'F':7,'Be':2,'Li':1,'Na':1,'Mg':2,'Al':3,'Si':4,'P':5,'S':6,'Cl':7,'K':1,'Ca':2,'Ga':3,'Ge':4,'As':5,'Se':6,'Br':7,'I':7,'Xe':8,'Kr':8}
    ELECTRONEGATIVITY = {'H':2.2,'B':2.04,'C':2.55,'N':3.04,'O':3.44,'F':3.98,'Be':1.57,'Li':0.98,'Na':0.93,'Mg':1.31,'Al':1.61,'Si':1.9,'P':2.19,'S':2.58,'Cl':3.16,'K':0.82,'Ca':1.0,'Ga':1.81,'Ge':2.01,'As':2.18,'Se':2.55,'Br':2.96,'I':2.66,'Xe':2.6,'Kr':3.0}

    def __init__(self, formula_string):
        self.is_valid = False
        self.atom_counts, self.charge = self._parse_formula(formula_string)
        if self.atom_counts is None: return
        self.total_valence_electrons = self._calculate_total_valence_electrons()
        self.all_valid_structures = []
        self.central_atom = ""
        self.terminal_atoms = []
        self.is_valid = True

    def _parse_formula(self, formula):
        atom_counts = {}
        charge = 0
        temp = formula
        if len(temp) >= 1 and (temp.endswith('-') or temp.endswith('+')):
            charge = -1 if temp.endswith('-') else 1
            temp = temp[:-1]
        elif len(temp) > 1 and temp[-1].isdigit() and (temp[-2] == '+' or temp[-2] == '-'):
            charge_val = int(temp[-1])
            charge = charge_val * (-1 if temp[-2] == '-' else 1)
            temp = temp[:-2]
        i = 0
        while i < len(temp):
            s = temp[i]; i += 1
            if i < len(temp) and temp[i].islower(): s += temp[i]; i += 1
            if s not in self.VALENCE_ELECTRONS:
                print("Error: Element '" + str(s) + "' not supported.")
                return None, None
            num = ""
            while i < len(temp) and temp[i].isdigit(): num += temp[i]; i += 1
            atom_counts[s] = atom_counts.get(s, 0) + (int(num) if num else 1)
        if not atom_counts:
            print("Error: No elements found.")
            return None, None
        return atom_counts, charge

    def _calculate_total_valence_electrons(self):
        total = 0
        for s, c in self.atom_counts.items(): total += self.VALENCE_ELECTRONS[s] * c
        return total - self.charge

    def _find_central_atom(self):
        keys = list(self.atom_counts.keys())
        if len(keys) == 1: return keys[0]
        candidates = {s: c for s, c in self.atom_counts.items() if s != 'H'}
        if not candidates: return min(keys, key=lambda s: self.ELECTRONEGATIVITY[s])
        singles = [s for s, c in candidates.items() if c == 1]
        if len(singles) == 1: return singles[0]
        return min(list(candidates.keys()), key=lambda s: self.ELECTRONEGATIVITY[s])

    def _generate_structures(self):
        if self.total_valence_electrons % 2 != 0: return
        c_symbol = self._find_central_atom()
        terms = []
        for s, cnt in self.atom_counts.items():
            num = cnt - 1 if s == c_symbol else cnt
            for _ in range(num): terms.append(s)
        
        self.central_atom = c_symbol + "0"
        self.terminal_atoms = []
        temp_counts = {}
        for s in terms:
            idx = temp_counts.get(s, 1) 
            self.terminal_atoms.append(s + str(idx))
            temp_counts[s] = idx + 1
        
        nb = len(self.terminal_atoms)
        e_in_b = nb * 2
        if e_in_b > self.total_valence_electrons: return
        bonds = {a: 1 for a in self.terminal_atoms}
        lp_electrons = self.total_valence_electrons - e_in_b
        lone_pairs = {a: 0 for a in [self.central_atom] + self.terminal_atoms}
        for t in self.terminal_atoms:
            sym = t[:-1] if t[-1].isdigit() else t
            need = 0 if sym == 'H' else 6
            take = min(lp_electrons, need)
            lone_pairs[t] += take; lp_electrons -= take
        lone_pairs[self.central_atom] += lp_electrons
        self._satisfy_central_octet(bonds, lone_pairs)

    def _satisfy_central_octet(self, bonds, lone_pairs):
        self._store_if_valid(bonds, lone_pairs)

        central_symbol = self.central_atom[:-1]
        bonds_electrons = sum(bonds.values()) * 2
        central_electrons = bonds_electrons + lone_pairs.get(self.central_atom, 0)
        
        is_period_2 = central_symbol in ['B', 'C', 'N', 'O', 'F']
        if is_period_2 and central_electrons >= 8:
            return

        for terminal_atom in self.terminal_atoms:
            if lone_pairs.get(terminal_atom, 0) >= 2 and bonds.get(terminal_atom, 1) < 3:
                
                new_bonds = bonds.copy()
                new_lone_pairs = lone_pairs.copy()

                new_bonds[terminal_atom] += 1
                new_lone_pairs[terminal_atom] -= 2
                
                self._satisfy_central_octet(new_bonds, new_lone_pairs)

    def _store_if_valid(self, bonds, lone_pairs):
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
        if structure not in self.all_valid_structures:
            self.all_valid_structures.append(structure)

    def get_optimal_structures(self):
        self._generate_structures()
        if not self.all_valid_structures: return None, []

        min_score = 100 
        for s in self.all_valid_structures:
            score = sum([abs(c) for c in s['formal_charges'].values()])
            if score < min_score:
                min_score = score

        best_score_structures = []
        for s in self.all_valid_structures:
            score = sum([abs(c) for c in s['formal_charges'].values()])
            if score == min_score:
                best_score_structures.append(s)

        if len(best_score_structures) <= 1:
            return best_score_structures[0], []

        final_optimal_structures = []
        max_electronegativity_for_neg_charge = -1.0 

        for s in best_score_structures:
            min_en_for_this_structure = 10.0 
            has_negative_charge = False
            for atom, charge in s['formal_charges'].items():
                if charge < 0:
                    has_negative_charge = True
                    symbol = atom[:-1] if atom[-1].isdigit() else atom
                    en = self.ELECTRONEGATIVITY[symbol]
                    if en < min_en_for_this_structure:
                        min_en_for_this_structure = en
            
            if not has_negative_charge:
                 min_en_for_this_structure = 10.0 

            if min_en_for_this_structure > max_electronegativity_for_neg_charge:
                max_electronegativity_for_neg_charge = min_en_for_this_structure
                final_optimal_structures = [s]
            elif min_en_for_this_structure == max_electronegativity_for_neg_charge:
                final_optimal_structures.append(s)

        if not final_optimal_structures:
            return best_score_structures[0], best_score_structures[1:]

        return final_optimal_structures[0], final_optimal_structures[1:]

    def format_structure(self, s, title):
        out = ["--- " + str(title) + " ---", "\n[Structure]"]
        bond_symbols = {1:'-', 2:'=', 3:'~'}
        for atom, order in s['bonds'].items():
            out.append("  " + self.central_atom + " " + bond_symbols.get(order, '?') + " " + atom)
        out.append("\n[Lone Pairs]")
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

def main():
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
        print("\nAn unexpected error occurred: " + str(e))

main()
print("\n=-=-=-=-=-=-=-=-=-=-=\n")
    
