from collections import Counter
from typing import Dict, List, Tuple
import re


def parse_fastq(content: str) -> List[Tuple[str, str, str]]:
    """Parse un fichier FASTQ et retourne (id, sequence, qualite)."""
    reads = []
    lines = [l.strip() for l in content.strip().split('\n') if l.strip()]
    i = 0
    while i + 3 < len(lines):
        header = lines[i]
        seq = lines[i + 1]
        plus = lines[i + 2]
        qual = lines[i + 3]
        if header.startswith('@') and plus.startswith('+'):
            reads.append((header[1:], seq, qual))
        i += 4
    return reads


def parse_fasta(content: str) -> List[Tuple[str, str]]:
    """Parse un fichier FASTA et retourne (id, sequence)."""
    reads = []
    current_id = None
    current_seq = []
    for line in content.strip().split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_id is not None:
                reads.append((current_id, ''.join(current_seq)))
            current_id = line[1:]
            current_seq = []
        elif line:
            current_seq.append(line.upper())
    if current_id is not None:
        reads.append((current_id, ''.join(current_seq)))
    return reads


def fastq_to_fasta(reads: List[Tuple[str, str, str]]) -> str:
    """Convertit des reads FASTQ en FASTA (filtre qualité minimale)."""
    fasta_lines = []
    for rid, seq, qual in reads:
        # Score Phred moyen
        scores = [ord(c) - 33 for c in qual]
        avg = sum(scores) / len(scores) if scores else 0
        if avg >= 20:  # Seuil Q20
            fasta_lines.append(f">{rid}")
            fasta_lines.append(seq)
    return '\n'.join(fasta_lines)


def extract_kmers(sequence: str, k: int) -> List[str]:
    """Extrait tous les k-mers d'une séquence."""
    seq = sequence.upper()
    # Filtre les bases ambiguës
    valid = set('ACGT')
    return [
        seq[i:i+k]
        for i in range(len(seq) - k + 1)
        if all(c in valid for c in seq[i:i+k])
    ]


def count_kmers(sequences: List[str], k: int) -> Counter:
    """Compte tous les k-mers dans une liste de séquences."""
    counts = Counter()
    for seq in sequences:
        counts.update(extract_kmers(seq, k))
    return counts


def estimate_error_rate(kmer_counts: Counter) -> dict:
    """
    Estime le taux d'erreur à partir de l'histogramme de fréquence.
    Les k-mers de fréquence 1 (hapax) sont souvent des erreurs de séquençage.
    """
    total = sum(kmer_counts.values())
    hapax = sum(1 for v in kmer_counts.values() if v == 1)
    solid_threshold = _find_valley(kmer_counts)
    
    error_kmers = sum(v for v in kmer_counts.values() if v < solid_threshold)
    solid_kmers = sum(v for v in kmer_counts.values() if v >= solid_threshold)
    
    return {
        'total_kmers': total,
        'distinct_kmers': len(kmer_counts),
        'hapax_count': hapax,
        'solid_threshold': solid_threshold,
        'error_kmers': error_kmers,
        'solid_kmers': solid_kmers,
        'estimated_error_rate': error_kmers / total if total > 0 else 0,
    }


def _find_valley(counts: Counter, max_freq: int = 100) -> int:
    """
    Trouve la vallée entre le pic d'erreurs (fréq. basse) et le pic 
    de couverture (fréq. élevée) pour déterminer le seuil 'solide'.
    """
    hist = [0] * (max_freq + 1)
    for v in counts.values():
        if v <= max_freq:
            hist[v] += 1
    
    # Cherche le premier minimum local après freq=1
    for i in range(2, max_freq - 1):
        if hist[i] < hist[i-1] and hist[i] <= hist[i+1]:
            return i
    return 3  # Valeur par défaut raisonnable


def generate_toy_dataset(ref_sequence: str, num_reads: int, read_length: int,
                          error_rate: float = 0.01) -> str:
    """
    Génère un dataset FASTQ synthétique à partir d'une séquence de référence.
    Simule un taux d'erreur de 1% (substitutions aléatoires).
    """
    import random
    bases = list('ACGT')
    lines = []
    
    for i in range(num_reads):
        # Position de départ aléatoire (avec chevauchements)
        max_start = max(0, len(ref_sequence) - read_length)
        start = random.randint(0, max_start)
        read = list(ref_sequence[start:start + read_length])
        
        # Introduire des erreurs
        for j in range(len(read)):
            if random.random() < error_rate:
                read[j] = random.choice([b for b in bases if b != read[j]])
        
        seq = ''.join(read)
        # Score de qualité Phred ~30 (bon) avec variation
        qual = ''.join(chr(random.randint(55, 73)) for _ in seq)
        
        lines.append(f"@read_{i+1}_pos{start}")
        lines.append(seq)
        lines.append('+')
        lines.append(qual)
    
    return '\n'.join(lines)