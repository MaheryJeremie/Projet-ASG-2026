"""
Module Lot 2 : Alignement par Programmation Dynamique (variante LCS)
Implémente la Plus Longue Sous-Séquence Commune (LCS) pour valider
les chevauchements entre reads.
"""
import numpy as np
from typing import Tuple, List


def lcs_matrix(seq1: str, seq2: str) -> np.ndarray:
    """
    Construit la matrice LCS (Longest Common Subsequence).
    Complexité : O(n*m) temps et espace.
    """
    n, m = len(seq1), len(seq2)
    dp = np.zeros((n + 1, m + 1), dtype=np.int32)
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1] + 1
            else:
                dp[i][j] = max(dp[i-1][j], dp[i][j-1])
    
    return dp


def traceback_lcs(dp: np.ndarray, seq1: str, seq2: str) -> str:
    """Remonte la matrice pour retrouver la LCS."""
    i, j = len(seq1), len(seq2)
    lcs_chars = []
    
    while i > 0 and j > 0:
        if seq1[i-1] == seq2[j-1]:
            lcs_chars.append(seq1[i-1])
            i -= 1
            j -= 1
        elif dp[i-1][j] > dp[i][j-1]:
            i -= 1
        else:
            j -= 1
    
    return ''.join(reversed(lcs_chars))


def find_overlap(seq1: str, seq2: str) -> Tuple[int, int, int]:
    """
    Trouve la position et la longueur du meilleur chevauchement entre seq1 et seq2.
    Teste si la fin de seq1 correspond au début de seq2 (suffix-prefix overlap).
    Retourne: (overlap_length, score, position_in_seq2)
    """
    best_overlap = 0
    best_score = 0
    
    min_overlap = 10  # Chevauchement minimum significatif
    
    for overlap_len in range(min_overlap, min(len(seq1), len(seq2)) + 1):
        suffix = seq1[-overlap_len:]
        prefix = seq2[:overlap_len]
        
        # Score : nombre de bases identiques
        matches = sum(a == b for a, b in zip(suffix, prefix))
        score = matches / overlap_len  # Identité normalisée
        
        if score > best_score:
            best_score = score
            best_overlap = overlap_len
    
    return best_overlap, best_score, 0


def align_reads(seq1: str, seq2: str) -> dict:
    """
    Alignement complet entre deux reads.
    Retourne toutes les informations d'alignement.
    """
    seq1 = seq1.upper().strip()
    seq2 = seq2.upper().strip()
    
    # Matrice LCS
    dp = lcs_matrix(seq1, seq2)
    lcs = traceback_lcs(dp, seq1, seq2)
    lcs_score = len(lcs)
    
    # Score d'identité normalisé
    max_len = max(len(seq1), len(seq2))
    identity = lcs_score / max_len if max_len > 0 else 0
    
    # Chevauchement suffix-prefix
    overlap_len, overlap_score, _ = find_overlap(seq1, seq2)
    
    # Alignement visuel
    visual = _build_visual_alignment(seq1, seq2, lcs)
    
    return {
        'seq1': seq1,
        'seq2': seq2,
        'lcs': lcs,
        'lcs_length': lcs_score,
        'identity': identity,
        'overlap_length': overlap_len,
        'overlap_identity': overlap_score,
        'matrix': dp,
        'visual': visual,
        'dp_dimensions': f"{len(seq1)+1} × {len(seq2)+1}",
    }


def _build_visual_alignment(seq1: str, seq2: str, lcs: str) -> dict:
    """
    Construit une représentation visuelle de l'alignement.
    Marque les matches, mismatches et gaps.
    """
    # Alignement séquence 1 vs LCS
    aligned1, aligned2, match_line = [], [], []
    
    i, j, k = 0, 0, 0
    lcs_set_positions1 = _get_lcs_positions(seq1, lcs)
    lcs_set_positions2 = _get_lcs_positions(seq2, lcs)
    
    # Alignement simplifié : suffix-prefix pour visualisation du chevauchement
    overlap = _compute_overlap_alignment(seq1, seq2)
    
    return overlap


def _get_lcs_positions(seq: str, lcs: str) -> List[int]:
    """Retourne les positions dans seq qui participent à la LCS."""
    positions = []
    j = 0
    for i, c in enumerate(seq):
        if j < len(lcs) and c == lcs[j]:
            positions.append(i)
            j += 1
    return positions


def _compute_overlap_alignment(seq1: str, seq2: str) -> dict:
    """
    Calcule et retourne l'alignement de chevauchement lisible.
    """
    best_score = -1
    best_offset = 0
    
    # Test tous les offsets possibles
    for offset in range(-len(seq2) + 1, len(seq1)):
        matches = 0
        total = 0
        for i in range(len(seq1)):
            j = i - offset
            if 0 <= j < len(seq2):
                total += 1
                if seq1[i] == seq2[j]:
                    matches += 1
        
        if total >= 10:
            score = matches / total
            if score > best_score:
                best_score = score
                best_offset = offset
    
    # Construire les lignes d'alignement
    offset = best_offset
    start = min(0, offset)
    end = max(len(seq1), len(seq2) + offset)
    
    line1, line2, line_match = [], [], []
    
    for pos in range(start, end):
        i = pos
        j = pos - offset
        
        c1 = seq1[i] if 0 <= i < len(seq1) else '-'
        c2 = seq2[j] if 0 <= j < len(seq2) else '-'
        
        line1.append(c1)
        line2.append(c2)
        
        if c1 == '-' or c2 == '-':
            line_match.append(' ')
        elif c1 == c2:
            line_match.append('|')
        else:
            line_match.append('.')
    
    return {
        'seq1_aligned': ''.join(line1),
        'seq2_aligned': ''.join(line2),
        'match_line': ''.join(line_match),
        'offset': offset,
        'identity': best_score,
        'overlap_region': _get_overlap_region(seq1, seq2, offset),
    }


def _get_overlap_region(seq1: str, seq2: str, offset: int) -> Tuple[int, int]:
    """Retourne les bornes de la région de chevauchement."""
    start = max(0, offset)
    end = min(len(seq1), len(seq2) + offset)
    return (start, end)