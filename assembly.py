"""
Module Lot 3 : Moteur d'Assemblage Memory-Efficient (approche Minia 2)
- Filtre de Bloom pour stocker les k-mers "solides"
- Traversée "on-the-fly" du Graphe de de Bruijn implicite
"""
import mmh3
import math
from typing import List, Set, Tuple, Optional
from collections import Counter


class BloomFilter:
    """
    Filtre de Bloom pour test d'appartenance des k-mers.
    
    Paramètres :
    - m : nombre de bits du filtre
    - k_hash : nombre de fonctions de hachage
    
    Taux de faux positifs théorique : (1 - e^(-k_hash*n/m))^k_hash
    où n est le nombre d'éléments insérés.
    """
    
    def __init__(self, capacity: int, false_positive_rate: float = 0.01):
        """
        Initialise le filtre de Bloom avec la capacité souhaitée
        et le taux de faux positifs cible.
        """
        # Formule optimale pour m (bits)
        self.m = self._optimal_m(capacity, false_positive_rate)
        # Formule optimale pour k_hash
        self.k_hash = self._optimal_k(self.m, capacity)
        
        # Tableau de bits (stocké comme bytearray pour l'efficacité)
        self.bits = bytearray(math.ceil(self.m / 8))
        self.n_inserted = 0
        self.capacity = capacity
        self.target_fpr = false_positive_rate
    
    @staticmethod
    def _optimal_m(n: int, p: float) -> int:
        """m = -n * ln(p) / (ln 2)^2"""
        return max(1, int(-n * math.log(p) / (math.log(2) ** 2)))
    
    @staticmethod
    def _optimal_k(m: int, n: int) -> int:
        """k = (m/n) * ln(2)"""
        if n == 0:
            return 1
        return max(1, round((m / n) * math.log(2)))
    
    def _get_bit_positions(self, item: str) -> List[int]:
        """Génère k_hash positions de bits via des fonctions de hachage MurmurHash3."""
        positions = []
        for seed in range(self.k_hash):
            h = mmh3.hash(item, seed, signed=False)
            positions.append(h % self.m)
        return positions
    
    def add(self, item: str):
        """Insère un élément dans le filtre."""
        for pos in self._get_bit_positions(item):
            byte_idx = pos // 8
            bit_idx = pos % 8
            self.bits[byte_idx] |= (1 << bit_idx)
        self.n_inserted += 1
    
    def contains(self, item: str) -> bool:
        """
        Teste l'appartenance. 
        Peut retourner True pour des éléments absents (faux positifs).
        Ne retourne JAMAIS False pour des éléments présents.
        """
        for pos in self._get_bit_positions(item):
            byte_idx = pos // 8
            bit_idx = pos % 8
            if not (self.bits[byte_idx] & (1 << bit_idx)):
                return False
        return True
    
    def estimated_fpr(self) -> float:
        """Taux de faux positifs réel estimé selon les éléments insérés."""
        if self.m == 0:
            return 1.0
        return (1 - math.exp(-self.k_hash * self.n_inserted / self.m)) ** self.k_hash
    
    def memory_bytes(self) -> int:
        """Mémoire utilisée en octets."""
        return math.ceil(self.m / 8)
    
    def stats(self) -> dict:
        return {
            'm_bits': self.m,
            'k_hash_functions': self.k_hash,
            'n_inserted': self.n_inserted,
            'memory_bytes': self.memory_bytes(),
            'memory_kb': self.memory_bytes() / 1024,
            'target_fpr': self.target_fpr,
            'estimated_fpr': self.estimated_fpr(),
        }


def build_bloom_filter(sequences: List[str], k: int, 
                        min_count: int = 2,
                        fpr: float = 0.01) -> Tuple[BloomFilter, Counter]:
    """
    Construit un Filtre de Bloom avec les k-mers "solides".
    Un k-mer est solide s'il apparaît >= min_count fois (filtre les erreurs).
    """
    # Étape 1 : Compter tous les k-mers
    kmer_counts = Counter()
    valid_bases = set('ACGT')
    
    for seq in sequences:
        seq = seq.upper()
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if all(c in valid_bases for c in kmer):
                kmer_counts[kmer] += 1
    
    # Étape 2 : Filtrer les k-mers solides
    solid_kmers = {km for km, cnt in kmer_counts.items() if cnt >= min_count}
    
    # Étape 3 : Construire le Filtre de Bloom
    bf = BloomFilter(capacity=max(len(solid_kmers), 1), false_positive_rate=fpr)
    for kmer in solid_kmers:
        bf.add(kmer)
    
    return bf, kmer_counts


def get_extensions(kmer: str) -> List[str]:
    """
    Génère les 4 extensions possibles en 3' d'un k-mer.
    Ex: "ACGT" → ["CGTA", "CGTC", "CGTG", "CGTT"]
    """
    prefix = kmer[1:]  # Retire la première base
    return [prefix + base for base in 'ACGT']


def get_predecessors(kmer: str) -> List[str]:
    """
    Génère les 4 prédécesseurs possibles en 5' d'un k-mer.
    """
    suffix = kmer[:-1]  # Retire la dernière base
    return [base + suffix for base in 'ACGT']


def traverse_graph(seed: str, bloom: BloomFilter, 
                    max_length: int = 10000,
                    max_branches: int = 3) -> List[str]:
    """
    Traversée "on-the-fly" du Graphe de de Bruijn implicite.
    
    Le graphe n'est JAMAIS construit explicitement.
    Le Filtre de Bloom sert d'oracle pour tester la connectivité.
    
    Retourne une liste de contigs assemblés.
    """
    contigs = []
    visited_seeds = set()
    
    def extend_contig(start_kmer: str) -> str:
        """Étend un contig depuis un k-mer de départ."""
        contig = list(start_kmer)
        current = start_kmer
        
        for _ in range(max_length):
            extensions = get_extensions(current)
            valid_ext = [e for e in extensions if bloom.contains(e)]
            
            if len(valid_ext) == 0:
                break  # Fin du contig
            elif len(valid_ext) == 1:
                # Chemin unique : continuer
                next_kmer = valid_ext[0]
                contig.append(next_kmer[-1])  # Ajoute seulement la nouvelle base
                current = next_kmer
            else:
                # Bifurcation : marquer et s'arrêter
                break
        
        return ''.join(contig)
    
    # Démarrer depuis le seed
    if bloom.contains(seed):
        contig = extend_contig(seed)
        if len(contig) >= len(seed):
            contigs.append(contig)
            visited_seeds.add(seed)
    
    return contigs


def assemble(sequences: List[str], k: int, 
             min_count: int = 2, fpr: float = 0.01,
             max_contigs: int = 20) -> dict:
    """
    Pipeline d'assemblage complet.
    """
    # Construire le Filtre de Bloom
    bloom, kmer_counts = build_bloom_filter(sequences, k, min_count, fpr)
    
    # Collecter les k-mers solides comme seeds potentiels
    valid_bases = set('ACGT')
    solid_seeds = [
        km for km, cnt in kmer_counts.items() 
        if cnt >= min_count and all(c in valid_bases for c in km)
    ]
    
    # Traversée on-the-fly depuis les meilleurs seeds (plus fréquents)
    solid_seeds.sort(key=lambda x: -kmer_counts[x])
    
    contigs = []
    seen_sequences = set()
    
    for seed in solid_seeds[:500]:  # Limiter les seeds explorés
        if len(contigs) >= max_contigs:
            break
        
        result = traverse_graph(seed, bloom, max_length=5000)
        
        for contig in result:
            # Dédupliquer les contigs qui se chevauchent
            is_duplicate = False
            for existing in seen_sequences:
                if contig in existing or existing in contig:
                    is_duplicate = True
                    break
            
            if not is_duplicate and len(contig) >= k:
                contigs.append(contig)
                seen_sequences.add(contig)
    
    # Trier par longueur décroissante
    contigs.sort(key=len, reverse=True)
    
    # Statistiques d'assemblage
    total_bases = sum(len(c) for c in contigs)
    n50 = _compute_n50(contigs)
    
    # Comparaison mémoire Bloom vs dictionnaire
    dict_memory = sum(len(km.encode()) + 28 for km in kmer_counts)  # Estimation dict Python
    
    return {
        'contigs': contigs,
        'bloom_stats': bloom.stats(),
        'kmer_counts': kmer_counts,
        'solid_kmers_count': len(solid_seeds),
        'total_bases_assembled': total_bases,
        'n50': n50,
        'num_contigs': len(contigs),
        'dict_memory_bytes': dict_memory,
        'bloom_memory_bytes': bloom.memory_bytes(),
        'memory_ratio': dict_memory / bloom.memory_bytes() if bloom.memory_bytes() > 0 else 0,
    }


def _compute_n50(contigs: List[str]) -> int:
    """Calcule le N50 : longueur au-delà de laquelle 50% des bases sont couvertes."""
    if not contigs:
        return 0
    lengths = sorted([len(c) for c in contigs], reverse=True)
    total = sum(lengths)
    cumsum = 0
    for l in lengths:
        cumsum += l
        if cumsum >= total / 2:
            return l
    return lengths[-1]


def compare_memory(kmer_counts: Counter, bloom: BloomFilter) -> dict:
    """
    Analyse comparative de la consommation mémoire :
    Filtre de Bloom vs dictionnaire Python standard.
    """
    n = len(kmer_counts)
    
    # Mémoire dictionnaire Python (estimation)
    # Chaque entrée : ~56 bytes overhead + longueur de la clé + valeur int (~28 bytes)
    avg_key_len = sum(len(k) for k in list(kmer_counts.keys())[:100]) / min(100, n) if n > 0 else 0
    dict_mem = n * (56 + avg_key_len + 28)
    
    bloom_mem = bloom.memory_bytes()
    
    return {
        'n_distinct_kmers': n,
        'dict_memory_bytes': int(dict_mem),
        'dict_memory_mb': dict_mem / (1024 * 1024),
        'bloom_memory_bytes': bloom_mem,
        'bloom_memory_mb': bloom_mem / (1024 * 1024),
        'compression_ratio': dict_mem / bloom_mem if bloom_mem > 0 else 0,
    }