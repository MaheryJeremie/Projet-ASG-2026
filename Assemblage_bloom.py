"""
Lot 3 : Moteur d'Assemblage De Novo — Filtre de Bloom + Graphe de de Bruijn Implicite
"""
import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from modules.assembly import (
    BloomFilter, build_bloom_filter, traverse_graph,
    assemble, compare_memory, get_extensions
)
from modules.kmer import generate_toy_dataset, parse_fastq, count_kmers

st.set_page_config(page_title="Lot 3 — Assemblage Bloom", page_icon="⚙️", layout="wide")

st.markdown("""
<style>
[data-testid="stAppViewContainer"] { background: #0d1117; color: #e6edf3; }
[data-testid="stSidebar"] { background: #161b22; border-right: 1px solid #30363d; }
h1, h2, h3 { color: #58a6ff; }
.contig-card {
    background: #161b22;
    border: 1px solid #30363d;
    border-left: 4px solid #3fb950;
    border-radius: 6px;
    padding: 0.8rem 1rem;
    margin-bottom: 0.6rem;
    font-family: monospace;
    font-size: 0.8rem;
}
.bloom-stat {
    background: #21262d;
    border-radius: 6px;
    padding: 0.6rem 1rem;
    margin: 0.3rem 0;
}
</style>
""", unsafe_allow_html=True)

st.markdown("## ⚙️ Lot 3 — Moteur d'Assemblage Memory-Efficient")
st.markdown(
    "Reconstruction de contigs via un **Filtre de Bloom** comme oracle d'appartenance "
    "et traversée **on-the-fly** du graphe de de Bruijn implicite (approche Minia 2)."
)
st.divider()

# Séquence référence
TOY_REF = (
    "ATGGAGTTCAACACCAATGTCACAGGCGAATCGTACAACCAGCACATCATTGTGGCTGTGATCAACATCATCAAC"
    "GGTTTGAATGTGGACTTCAAGGGCGAGTACTTGTACACGGCGAAGATCCCTGCGGTCGAGTTCCGCAAGCAGGTC"
    "GAGTTCAAGTACGCCAACAAGGTCACCGTCACGTACACCAACAGCGGCGTCGTCATCGCCGGCGAGATCATCAAG"
    "GCCATCAAGGAGACGCTGAACGCCGCGGTGAACAACATCGTCGACGGCAAGATCATCGCCAACGACTTCACCAAG"
    "TTCGAGCGCAAGGCGATCAAGGAGTTCGGCGAGAACAAGGTCGACATCGCCAAGCAGTTCGTCATCGACGAGAAC"
)

# ─── Sidebar : paramètres ────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("### ⚙️ Paramètres d'assemblage")
    k_assembly = st.slider("k (taille des k-mers)", 11, 31, 21, step=2,
                            help="k plus grand → plus spécifique, moins de faux positifs dans le graphe")
    min_count = st.slider("Seuil 'solide' (min count)", 1, 10, 2,
                           help="k-mers apparaissant < min_count sont considérés comme des erreurs")
    fpr_target = st.select_slider("Taux de faux positifs cible (Bloom)",
                                    options=[0.001, 0.005, 0.01, 0.05, 0.1],
                                    value=0.01,
                                    format_func=lambda x: f"{x*100:.1f}%")
    
    st.divider()
    st.markdown("### 📖 Principe")
    st.markdown("""
    **Graphe de de Bruijn implicite :**
    - Nœuds = k-mers
    - Arêtes = chevauchements de (k-1) bases
    
    **Sans jamais construire le graphe**, le Filtre de Bloom répond à la question :
    *"Ce k-mer est-il dans mes données ?"*
    
    **Traversée :**
    ```
    kmer → 4 extensions possibles
         → filtre Bloom
         → 1 valide → continuer
         → 0 valide → fin contig  
         → 2+ valides → bifurcation
    ```
    """)

# ─── Onglets ─────────────────────────────────────────────────────────────────
tab1, tab2, tab3, tab4 = st.tabs([
    "🔵 Filtre de Bloom",
    "⚙️ Assemblage",
    "📊 Analyse mémoire",
    "🔬 Analyse critique",
])

# ══════════════════════════════════════════════════════════════════════════════
with tab1:
    st.markdown("### 🔵 Construction du Filtre de Bloom")
    st.markdown(
        "Le Filtre de Bloom est une structure de données probabiliste qui stocke les "
        "k-mers **solides** de manière compacte pour tester l'appartenance en O(k_hash)."
    )
    
    # Données source
    data_source = st.radio("Source de données", ["🧪 Toy Dataset", "📋 Données du Lot 1"],
                            horizontal=True)
    
    sequences = []
    if "🧪 Toy" in data_source:
        n_reads = st.number_input("Nombre de reads", 100, 5000, 1000, 100)
        read_len = st.number_input("Longueur des reads", 50, 150, 75, 5)
        err_rate = st.slider("Taux d'erreur", 0.0, 0.05, 0.01, 0.001)
        
        if st.button("Générer & construire le filtre", type="primary"):
            with st.spinner("Génération et construction du filtre..."):
                fq = generate_toy_dataset(TOY_REF, n_reads, read_len, err_rate)
                reads = parse_fastq(fq)
                sequences = [r[1] for r in reads]
                st.session_state['assembly_seqs'] = sequences
    else:
        if 'reads_data' in st.session_state and st.session_state['reads_data']:
            sequences = [r[1] for r in st.session_state['reads_data']]
            st.info(f"✅ {len(sequences)} reads chargés depuis le Lot 1.")
            st.session_state['assembly_seqs'] = sequences
        else:
            st.warning("Aucune donnée dans le Lot 1. Utilisez le Toy Dataset.")
    
    if 'assembly_seqs' in st.session_state and st.session_state['assembly_seqs']:
        seqs = st.session_state['assembly_seqs']
        
        with st.spinner("Construction du Filtre de Bloom..."):
            bloom, kmer_counts = build_bloom_filter(seqs, k_assembly, min_count, fpr_target)
            st.session_state['bloom'] = bloom
            st.session_state['kmer_counts_asm'] = kmer_counts
            st.session_state['assembly_k'] = k_assembly
        
        bstats = bloom.stats()
        
        st.success("✅ Filtre de Bloom construit")
        
        # Métriques du filtre
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("m (bits)", f"{bstats['m_bits']:,}")
        c2.metric("k_hash (fonctions)", bstats['k_hash_functions'])
        c3.metric("Mémoire", f"{bstats['memory_kb']:.1f} KB")
        c4.metric("FPR estimé", f"{bstats['estimated_fpr']*100:.3f}%")
        
        col_left, col_right = st.columns(2)
        
        with col_left:
            st.markdown("#### Détail du filtre")
            st.markdown(f"""
            <div class="bloom-stat">📐 Capacité cible : <b>{bstats['n_inserted']:,}</b> k-mers solides</div>
            <div class="bloom-stat">🎯 FPR cible : <b>{bstats['target_fpr']*100:.1f}%</b></div>
            <div class="bloom-stat">✅ FPR réel estimé : <b>{bstats['estimated_fpr']*100:.3f}%</b></div>
            <div class="bloom-stat">💾 m = <b>{bstats['m_bits']:,}</b> bits = <b>{bstats['memory_kb']:.1f}</b> KB</div>
            <div class="bloom-stat">🔑 k_hash = <b>{bstats['k_hash_functions']}</b> fonctions MurmurHash3</div>
            """, unsafe_allow_html=True)
        
        with col_right:
            # Visualisation du taux de remplissage
            st.markdown("#### Taux de remplissage")
            total_bits = bstats['m_bits']
            # Estimation bits à 1
            bits_set = int(total_bits * (1 - (1 - 1/total_bits)**(bstats['k_hash_functions']*bstats['n_inserted'])))
            fill_rate = bits_set / total_bits if total_bits > 0 else 0
            
            fig, ax = plt.subplots(figsize=(5, 2), facecolor='#0d1117')
            ax.set_facecolor('#161b22')
            ax.barh([0], [fill_rate], color='#3fb950', height=0.5)
            ax.barh([0], [1 - fill_rate], left=[fill_rate], color='#21262d', height=0.5)
            ax.set_xlim(0, 1)
            ax.set_yticks([])
            ax.set_xlabel("Taux de remplissage du tableau de bits", color='#8b949e')
            ax.tick_params(colors='#8b949e')
            ax.spines[:].set_color('#30363d')
            ax.set_title(f"{fill_rate*100:.1f}% des bits à 1", color='#3fb950')
            plt.tight_layout()
            st.pyplot(fig)
            plt.close()
            
            st.caption("Plus le filtre est rempli, plus le FPR augmente.")
        
        # Démonstration de la vérification d'appartenance
        st.divider()
        st.markdown("#### 🔍 Test d'appartenance interactif")
        test_kmer = st.text_input(f"Saisir un {k_assembly}-mer à tester",
                                   placeholder=f"Ex: {'A'*k_assembly}")
        if test_kmer:
            test_kmer = test_kmer.upper().strip()
            if len(test_kmer) == k_assembly:
                result = bloom.contains(test_kmer)
                count = kmer_counts.get(test_kmer, 0)
                
                if result and count >= min_count:
                    st.success(f"✅ **PRÉSENT** dans le filtre (fréquence réelle : {count}×) — Vrai Positif")
                elif result and count < min_count:
                    st.error(f"⚠️ **PRÉSENT dans le filtre** MAIS fréquence réelle = {count} — **FAUX POSITIF !**")
                else:
                    st.info(f"❌ **ABSENT** du filtre (fréquence réelle : {count}) — Vrai Négatif")
            else:
                st.warning(f"Le k-mer doit avoir exactement {k_assembly} caractères.")

# ══════════════════════════════════════════════════════════════════════════════
with tab2:
    st.markdown("### ⚙️ Assemblage De Novo")
    
    if 'bloom' not in st.session_state:
        st.info("Construisez d'abord le Filtre de Bloom dans l'onglet précédent.")
    else:
        bloom = st.session_state['bloom']
        kmer_counts_asm = st.session_state['kmer_counts_asm']
        k_asm = st.session_state['assembly_k']
        seqs_asm = st.session_state['assembly_seqs']
        
        col_seed, col_max = st.columns(2)
        max_contigs = col_seed.number_input("Nombre max de contigs", 1, 50, 10)
        
        if st.button("🚀 Lancer l'assemblage", type="primary"):
            with st.spinner("Traversée on-the-fly du graphe de de Bruijn implicite..."):
                result = assemble(seqs_asm, k_asm, min_count, fpr_target, max_contigs)
                st.session_state['assembly_result'] = result
        
        if 'assembly_result' in st.session_state:
            res = st.session_state['assembly_result']
            contigs = res['contigs']
            
            # Métriques d'assemblage
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Contigs assemblés", res['num_contigs'])
            c2.metric("Bases totales", f"{res['total_bases_assembled']:,} bp")
            c3.metric("N50", f"{res['n50']:,} bp")
            c4.metric("Ratio compression", f"{res['memory_ratio']:.0f}×",
                      help="Bloom vs dictionnaire Python")
            
            st.divider()
            
            # Validation vs référence
            if st.session_state.get('ref_seq'):
                ref = st.session_state['ref_seq']
                st.markdown("#### 🎯 Validation vs séquence de référence")
                
                if contigs:
                    longest = contigs[0]
                    # Calcul d'identité simplifié : pourcentage de k-mers du contig dans la référence
                    from modules.kmer import extract_kmers
                    ref_kmers = set(extract_kmers(ref, k_asm))
                    contig_kmers = set(extract_kmers(longest, k_asm))
                    
                    if contig_kmers:
                        identity = len(contig_kmers & ref_kmers) / len(contig_kmers)
                        recall = len(contig_kmers & ref_kmers) / len(ref_kmers) if ref_kmers else 0
                        
                        col_v1, col_v2, col_v3 = st.columns(3)
                        col_v1.metric("Identité (precision k-mers)", f"{identity*100:.1f}%",
                                       delta="✅ Validé" if identity >= 0.98 else "⚠️ Sous le seuil")
                        col_v2.metric("Rappel (k-mers ref couverts)", f"{recall*100:.1f}%")
                        col_v3.metric("Seuil de validation", "98%")
                        
                        if identity >= 0.98:
                            st.success("🎉 **Critère de validation atteint** : identité ≥ 98%")
                        else:
                            st.warning(f"Identité {identity*100:.1f}% < 98%. Essayez d'ajuster k ou min_count.")
            
            # Liste des contigs
            st.markdown("#### 📋 Contigs assemblés")
            
            for i, contig in enumerate(contigs):
                gc = (contig.count('G') + contig.count('C')) / len(contig) * 100 if contig else 0
                
                st.markdown(f"""
                <div class="contig-card">
                <b style="color:#58a6ff">Contig {i+1}</b>
                &nbsp;|&nbsp; <span style="color:#3fb950">{len(contig)} bp</span>
                &nbsp;|&nbsp; GC: {gc:.1f}%
                <br/><span style="color:#8b949e; word-break:break-all;">{contig[:120]}{'...' if len(contig) > 120 else ''}</span>
                </div>
                """, unsafe_allow_html=True)
            
            # Graphique longueurs
            if contigs:
                st.markdown("#### 📊 Distribution des longueurs de contigs")
                fig, ax = plt.subplots(figsize=(8, 3), facecolor='#0d1117')
                ax.set_facecolor('#161b22')
                lengths = [len(c) for c in contigs]
                ax.bar(range(1, len(lengths)+1), lengths, color='#3fb950', alpha=0.8)
                ax.axhline(y=res['n50'], color='#e3b341', linestyle='--',
                           label=f'N50 = {res["n50"]} bp')
                ax.set_xlabel("Contig #", color='#8b949e')
                ax.set_ylabel("Longueur (bp)", color='#8b949e')
                ax.set_title("Longueurs des contigs assemblés", color='#58a6ff')
                ax.tick_params(colors='#8b949e')
                ax.spines[:].set_color('#30363d')
                ax.legend(facecolor='#21262d', edgecolor='#30363d', labelcolor='#e6edf3')
                plt.tight_layout()
                st.pyplot(fig)
                plt.close()
            
            # Export FASTA
            fasta_out = '\n'.join(
                f">contig_{i+1}_len{len(c)}\n{c}"
                for i, c in enumerate(contigs)
            )
            st.download_button("⬇️ Exporter les contigs (FASTA)", fasta_out,
                               "contigs.fasta", "text/plain")

# ══════════════════════════════════════════════════════════════════════════════
with tab3:
    st.markdown("### 📊 Analyse comparative de la mémoire")
    st.markdown(
        "Justification du choix du Filtre de Bloom par rapport à un dictionnaire Python standard."
    )
    
    if 'bloom' not in st.session_state:
        st.info("Construisez le Filtre de Bloom dans le premier onglet.")
    else:
        bloom = st.session_state['bloom']
        kmer_counts_asm = st.session_state['kmer_counts_asm']
        
        mem = compare_memory(kmer_counts_asm, bloom)
        
        c1, c2, c3 = st.columns(3)
        c1.metric("k-mers distincts", f"{mem['n_distinct_kmers']:,}")
        c2.metric("Dictionnaire Python", f"{mem['dict_memory_mb']:.2f} MB")
        c3.metric("Filtre de Bloom", f"{mem['bloom_memory_mb']:.2f} MB",
                  delta=f"÷ {mem['compression_ratio']:.0f}× plus léger")
        
        # Graphique comparatif
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), facecolor='#0d1117')
        
        # Bar chart
        ax = axes[0]
        ax.set_facecolor('#161b22')
        vals = [mem['dict_memory_mb'], mem['bloom_memory_mb']]
        labels = ['Dictionnaire\nPython', 'Filtre de\nBloom']
        colors = ['#f85149', '#3fb950']
        bars = ax.bar(labels, vals, color=colors, alpha=0.85, width=0.4, edgecolor='none')
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(vals)*0.01,
                    f'{val:.2f} MB', ha='center', va='bottom', color='#e6edf3', fontsize=10)
        ax.set_ylabel("Mémoire (MB)", color='#8b949e')
        ax.set_title("Comparaison mémoire", color='#58a6ff')
        ax.tick_params(colors='#8b949e')
        ax.spines[:].set_color('#30363d')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Scalabilité : projection pour différentes tailles
        ax2 = axes[1]
        ax2.set_facecolor('#161b22')
        n_range = np.logspace(3, 8, 50)
        k_len = st.session_state['assembly_k']
        
        # Mémoire dict : ~84 bytes par k-mer (overhead + chaîne + int)
        dict_mem_curve = n_range * (84 + k_len) / (1024**2)
        # Mémoire Bloom : -n*ln(p)/(ln2)^2 bits
        import math
        fpr = fpr_target
        bloom_mem_curve = (-n_range * math.log(fpr) / (math.log(2)**2)) / (8 * 1024**2)
        
        ax2.plot(n_range, dict_mem_curve, color='#f85149', label='Dictionnaire Python', linewidth=2)
        ax2.plot(n_range, bloom_mem_curve, color='#3fb950', label='Filtre de Bloom', linewidth=2)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_xlabel("Nombre de k-mers distincts", color='#8b949e')
        ax2.set_ylabel("Mémoire (MB)", color='#8b949e')
        ax2.set_title("Scalabilité (projection log-log)", color='#58a6ff')
        ax2.tick_params(colors='#8b949e')
        ax2.spines[:].set_color('#30363d')
        ax2.legend(facecolor='#21262d', edgecolor='#30363d', labelcolor='#e6edf3', fontsize=9)
        ax2.grid(color='#21262d', alpha=0.5)
        
        plt.tight_layout()
        st.pyplot(fig)
        plt.close()
        
        st.markdown(f"""
        **Interprétation** : Pour **{mem['n_distinct_kmers']:,}** k-mers de taille {k_assembly},
        le Filtre de Bloom est **{mem['compression_ratio']:.0f}× plus économe** en mémoire.
        En échelle logarithmique, l'écart croît avec le volume de données,
        ce qui justifie son utilisation pour des projets à grande échelle (millions de k-mers).
        
        **Complexité spatiale :**
        - Dictionnaire : $O(n \\cdot k)$ où n = k-mers distincts, k = longueur
        - Filtre de Bloom : $O(m)$ bits, indépendant de k, avec $m = -n \\cdot \\ln(p) / (\\ln 2)^2$
        """)

# ══════════════════════════════════════════════════════════════════════════════
with tab4:
    st.markdown("### 🔬 Analyse critique : Impact des faux positifs")
    
    st.markdown("""
    Le Filtre de Bloom peut retourner **vrai** pour un k-mer absent des données réelles.
    Ces **faux positifs** créent des **chemins fantômes** dans le graphe de de Bruijn implicite.
    """)
    
    st.markdown("#### Impact sur l'assemblage")
    
    col_l, col_r = st.columns(2)
    with col_l:
        st.markdown("""
        **Conséquences des faux positifs :**
        
        1. **Chemins fantômes** : une extension invalide est acceptée, créant une branche fictive dans le graphe
        2. **Chimères** : deux contigs distincts sont incorrectement fusionnés
        3. **Bifurcations artificielles** : une extension valide génère aussi une fausse extension, forçant l'arrêt prématuré du contig
        4. **Réduction du N50** : les contigs sont plus courts et plus nombreux que nécessaire
        """)
    
    with col_r:
        st.markdown("""
        **Paramètres clés et compromis :**
        
        | Paramètre | ↑ Effet | ↓ Effet |
        |-----------|---------|---------|
        | m (bits) | ↓ FPR, ↑ mémoire | ↑ FPR, ↓ mémoire |
        | k_hash | ↓ FPR (optimal) | — |
        | k (k-mer) | ↓ spécificité | ↑ sensibilité |
        | min_count | ↓ bruit | ↑ perte signal |
        """)
    
    # Simulation de l'impact du FPR sur la qualité d'assemblage
    st.markdown("#### 📈 Simulation : FPR vs qualité d'assemblage théorique")
    
    fpr_vals = np.logspace(-4, -0.5, 100)
    
    # Modèle simplifié : probabilité d'une bifurcation fantôme par extension
    # P(au moins un faux positif parmi 4 extensions) = 1 - (1-FPR)^4
    p_phantom_branch = 1 - (1 - fpr_vals)**4
    
    # Longueur moyenne d'un contig sans bifurcation fantôme
    # L_contig = k / P(bifurcation) (approximation géométrique)
    k_val = 21
    mean_contig_len = k_val / (p_phantom_branch + 0.001)  # Éviter div/0
    mean_contig_len = np.clip(mean_contig_len, 0, 5000)
    
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), facecolor='#0d1117')
    
    ax1 = axes[0]
    ax1.set_facecolor('#161b22')
    ax1.semilogx(fpr_vals * 100, p_phantom_branch * 100, color='#f85149', linewidth=2)
    ax1.axvline(x=1.0, color='#e3b341', linestyle='--', label='FPR = 1%')
    ax1.axvline(x=0.1, color='#3fb950', linestyle='--', label='FPR = 0.1%')
    ax1.set_xlabel("FPR (%)", color='#8b949e')
    ax1.set_ylabel("P(branche fantôme) %", color='#8b949e')
    ax1.set_title("Probabilité de chemin fantôme", color='#58a6ff')
    ax1.tick_params(colors='#8b949e')
    ax1.spines[:].set_color('#30363d')
    ax1.legend(facecolor='#21262d', edgecolor='#30363d', labelcolor='#e6edf3', fontsize=9)
    ax1.grid(color='#21262d', alpha=0.5)
    
    ax2 = axes[1]
    ax2.set_facecolor('#161b22')
    ax2.semilogx(fpr_vals * 100, mean_contig_len, color='#3fb950', linewidth=2)
    ax2.axvline(x=1.0, color='#e3b341', linestyle='--', label='FPR = 1%')
    ax2.set_xlabel("FPR (%)", color='#8b949e')
    ax2.set_ylabel("Longueur contig théorique (bp)", color='#8b949e')
    ax2.set_title("Impact sur la longueur des contigs (N50)", color='#58a6ff')
    ax2.tick_params(colors='#8b949e')
    ax2.spines[:].set_color('#30363d')
    ax2.legend(facecolor='#21262d', edgecolor='#30363d', labelcolor='#e6edf3', fontsize=9)
    ax2.grid(color='#21262d', alpha=0.5)
    
    plt.tight_layout()
    st.pyplot(fig)
    plt.close()
    
    st.markdown("#### 💡 Recommandations pratiques")
    st.info("""
    **Pour minimiser l'impact des faux positifs :**
    
    - **FPR cible ≤ 1%** est un bon compromis mémoire/qualité pour la plupart des projets
    - **Augmenter min_count** réduit les k-mers d'erreur et améliore le signal/bruit
    - **k plus grand** réduit la probabilité de collision, au prix d'une moindre sensibilité
    - **Post-traitement** : filtrer les contigs courts (<2k bases) après assemblage élimine les artefacts
    - L'algorithme Minia 2 original utilise en plus une **étape de correction** des branches courtes
    """)
    
    st.markdown("#### 📐 Formule du FPR")
    st.latex(r"\text{FPR} = \left(1 - e^{-k_h \cdot n / m}\right)^{k_h}")
    st.markdown("""
    Avec :
    - $k_h$ = nombre de fonctions de hachage  
    - $n$ = nombre d'éléments insérés  
    - $m$ = taille du filtre en bits
    
    La valeur optimale de $k_h$ est : $k_h^* = \\frac{m}{n} \\ln 2$
    """)