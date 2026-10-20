"""Resolve exact occupancy ties against fixed local peptide geometry.

Strict occupancy winners and unrelated alternate labels remain independent.
A tied choice changes only when its original geometry conflicts with a fixed
adjacent peptide and exactly one complete alternative satisfies those anchors.
Every final evaluable peptide edge touching a tie is checked again; unresolved
incompatible tied chains hold. These compatibility checks do not replace
downstream geometry validation.
"""
import hashlib
import math

STANDARD = frozenset({'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                      'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'})
POLICY = 'residue_max_occupancy_fixed_peptide_final_edge_guard_v3'


def resolve_occupancy_ties(rows, groups, selected, decisions, peptide_context=None):
    """Never alter a strict maximum, an atom coordinate or a nonlocal label.

    The 2 A C--N bound is existing SkinPipe fragment evidence. The 1.644 A
    O...N conflict is Meeko 0.8's 1.2 * (O 0.66 + N 0.71) bond perception.
    Both are refusal/compatibility criteria, not physical validation.
    """
    evidence = {'policy': POLICY, 'source_atom_records_sha256': hashlib.sha256('\n'.join(rows).encode()).hexdigest(),
                'context_available': peptide_context is not None, 'events': [], 'coordinate_changes': 0,
                'strictly_higher_occupancy_changes': 0, 'global_altloc_label_matching': False}
    if peptide_context is None:
        evidence['reason'] = 'No full TER/model context: retain occupancy/lexical choices'
        return selected, evidence
    evidence['context_sha256'] = hashlib.sha256(peptide_context.encode()).hexdigest()
    indexed = {(row['chain'], row['residue'], row['name']): row for row in decisions}
    by_key = {key: indexed.get((key[0], key[1].strip(), key[2].strip())) for key in groups}
    domains, choices = {}, {}
    for key in groups:
        decision = by_key.get(key)
        if decision is None:
            domains[key], choices[key] = [' '], ' '
        else:
            scores = decision['mean_heavy_occupancy']
            maximum = max(scores.values())
            domains[key] = sorted((label for label, score in scores.items() if score == maximum), key=lambda label: (label != 'A', label))
            choices[key] = decision['selected_altloc']

    def candidate(key, label):
        atoms = {}
        for index in groups[key]:
            line = rows[index]
            if line[16] not in {' ', label}:
                continue
            name = line[12:16].strip()
            if name in atoms:
                previous = rows[atoms[name]]
                if previous[16] == line[16]:
                    raise ValueError('duplicate atom in tied peptide alternative; selection held')
                winner = min((atoms[name], index), key=lambda item: (-float(rows[item][54:60]), rows[item][16] != ' '))
                atoms[name] = winner
            else:
                atoms[name] = index
        heavy = {name: index for name, index in atoms.items() if rows[index][76:78].strip().upper() not in {'H', 'D'}}
        xyz = {name: tuple(float(rows[index][start:start + 8]) for start in (30, 38, 46)) for name, index in heavy.items()}
        if any(not all(math.isfinite(value) for value in point) for point in xyz.values()):
            raise ValueError('nonfinite tied peptide coordinates; selection held')
        return set(atoms.values()), xyz

    # Full context supplies adjacency. Never infer across TER, chain changes,
    # missing/ambiguous numbering, nonstandard/hetero residues or extra models.
    segments, locations, segment = {}, {}, 0
    previous_chain = None
    for line in peptide_context.splitlines():
        if line.startswith('ENDMDL'):
            break
        if line.startswith('TER'):
            segment += 1
            previous_chain = None
            continue
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        if line[21] != previous_chain:
            segment += 1
            previous_chain = line[21]
        key = (line[21], line[22:27], line[17:20])
        if not line.startswith('ATOM  ') or key[2].strip() not in STANDARD:
            segment += 1
            previous_chain = None
            continue
        if key not in groups:
            continue
        order = segments.setdefault(segment, [])
        if not order or order[-1] != key:
            order.append(key)
        locations.setdefault(key, set()).add(segment)
    neighbors = {key: [] for key in groups}
    for order in segments.values():
        for left, right in zip(order, order[1:]):
            a, b = left[1].strip(), right[1].strip()
            if (len(locations[left]) != 1 or len(locations[right]) != 1 or not a.lstrip('-').isdigit()
                    or not b.lstrip('-').isdigit() or int(b) != int(a) + 1):
                continue
            neighbors[left].append((left, right))
            neighbors[right].append((left, right))

    updated = set(selected)
    for key, labels in domains.items():
        if len(labels) < 2:
            continue
        event = {'chain': key[0], 'residue': key[1].strip(), 'name': key[2].strip(),
                 'original_selected_altloc': choices[key], 'equal_maximum_labels': labels, 'candidate_edges': {}}
        candidates = {label: candidate(key, label) for label in labels}
        anchors = [(left, right) for left, right in neighbors[key] if len(domains[right if left == key else left]) == 1]
        for label, (_indices, coordinates) in candidates.items():
            edges = []
            for left, right in anchors:
                neighbor = right if left == key else left
                _fixed_indices, fixed = candidate(neighbor, choices[neighbor])
                first, second = (coordinates, fixed) if left == key else (fixed, coordinates)
                if not {'C', 'O'} <= first.keys() or 'N' not in second:
                    continue
                cn, on = math.dist(first['C'], second['N']), math.dist(first['O'], second['N'])
                edges.append({'left': [left[0], left[1].strip()], 'right': [right[0], right[1].strip()],
                              'C_N_A': cn, 'O_N_A': on, 'compatible': 0 < cn <= 2.0 and on >= 1.2 * (0.66 + 0.71)})
            event['candidate_edges'][label] = edges
        baseline = event['candidate_edges'][choices[key]]
        if not baseline or all(edge['compatible'] for edge in baseline):
            event['decision'] = 'retain_occupancy_lexical_choice_no_fixed_conflict'
            evidence['events'].append(event)
            continue
        same_coverage = len({frozenset(coords) for _ids, coords in candidates.values()}) == 1
        winners = [label for label in labels if len(event['candidate_edges'][label]) == len(baseline)
                   and all(edge['compatible'] for edge in event['candidate_edges'][label])]
        if not same_coverage or len(winners) != 1:
            raise ValueError('tied peptide alternatives have no unique full-coverage compatible choice; selection held')
        winner = winners[0]
        updated.difference_update(groups[key])
        updated.update(candidates[winner][0])
        decision = by_key[key]
        decision['original_occupancy_lexical_altloc'] = choices[key]
        decision['selected_altloc'] = winner
        decision['discarded_serials'] = [rows[index][6:11].strip() for index in groups[key] if index not in candidates[winner][0]]
        event.update(decision='unique_equal_occupancy_fixed_peptide_compatible_choice', selected_altloc=winner)
        evidence['events'].append(event)
    # Re-evaluate every adjacent edge touching a tie with the final choices.
    # Local unique choices can conflict with another tied residue; do not
    # silently publish that greedy mixture as a resolved peptide chain.
    final_edges = []
    seen_edges = set()
    for key, labels in domains.items():
        if len(labels) < 2:
            continue
        for left, right in neighbors[key]:
            if (left, right) in seen_edges:
                continue
            seen_edges.add((left, right))
            final_left = by_key[left]['selected_altloc'] if by_key[left] else choices[left]
            final_right = by_key[right]['selected_altloc'] if by_key[right] else choices[right]
            _, first = candidate(left, final_left)
            _, second = candidate(right, final_right)
            if not {'C', 'O'} <= first.keys() or 'N' not in second:
                continue
            cn, on = math.dist(first['C'], second['N']), math.dist(first['O'], second['N'])
            edge = {'left': [left[0], left[1].strip()], 'right': [right[0], right[1].strip()],
                    'C_N_A': cn, 'O_N_A': on,
                    'compatible': 0 < cn <= 2.0 and on >= 1.2 * (0.66 + 0.71)}
            final_edges.append(edge)
            if not edge['compatible']:
                raise ValueError('final tied peptide choices have incompatible adjacent geometry; selection held')
    evidence['final_adjacent_edges'] = final_edges
    evidence['final_adjacent_edge_scope'] = 'all evaluable peptide edges touching equal-occupancy alternatives; compatibility only, not geometry validation'
    return updated, evidence
