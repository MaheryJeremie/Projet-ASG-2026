import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from modules.kmer import (
    parse_fastq, parse_fasta,
    count_kmers, estimate_error_rate,
    generate_toy_dataset,
)
from aff.ui import inject_styles, page_header, sidebar_brand, sidebar_section, lot_chart, LOT_COLOR

st.set_page_config(page_title="Lot 1 — Ingestion & Qualité", layout="wide")
inject_styles(lot=1)
CHART = lot_chart(1)

TOY_REF = (
    "ATGGAGTTCAACACCAATGTCACAGGCGAATCGTACAACCAGCACATCATTGTGGCTGTGATCAACATCATCAAC"
    "GGTTTGAATGTGGACTTCAAGGGCGAGTACTTGTACACGGCGAAGATCCCTGCGGTCGAGTTCCGCAAGCAGGTC"
    "GAGTTCAAGTACGCCAACAAGGTCACCGTCACGTACACCAACAGCGGCGTCGTCATCGCCGGCGAGATCATCAAG"
    "GCCATCAAGGAGACGCTGAACGCCGCGGTGAACAACATCGTCGACGGCAAGATCATCGCCAACGACTTCACCAAG"
    "TTCGAGCGCAAGGCGATCAAGGAGTTCGGCGAGAACAAGGTCGACATCGCCAAGCAGTTCGTCATCGACGAGAAC"
    "TACAAGATCGCCAAGCAGTTCATCGCCGACGAGTTCGTCATCAACGAGAACTTCAAGGTCGCCAAGCAGTTCATC"
    "GCGGACGAGTTCGTCATCAACGAGAACTTCAAGGTCGCCAAGCAGTTCATCGCGGACGAGTTCGTCATCAACGAG"
)

page_header(
    "Lot 1 — Ingestion et contrôle qualité",
    "Lecture des reads, comptage des k-mers, filtrage Phred, export FASTA.",
)

with st.sidebar:
    sidebar_brand(lot=1)
    sidebar_section("Paramètres")
    k = st.slider("k (taille des k-mers)", min_value=7, max_value=31, value=21, step=2)
    min_count = st.slider("Seuil min_count (k-mer solide)", min_value=1, max_value=10, value=2)
    max_freq_display = st.slider("Fréquence max (histogramme)", 10, 100, 50)

    sidebar_section("Rappels")
    st.markdown("""
    **k-mer** : sous-séquence de longueur k.

    **Hapax** : k-mer observé une seule fois.

    **Phred** : score de qualité de séquençage.

    **bp** : paire de bases (base pair).
    """)

tab1, tab2, tab3 = st.tabs(["Données", "k-mers", "FASTA"])

with tab1:
    data_source = st.radio(
        "Source",
        ["Toy dataset", "Fichier FASTQ ou FASTA"],
        horizontal=True,
    )

    if "Toy" in data_source:
        col_a, col_b = st.columns(2)
        with col_a:
            num_reads = st.number_input("Nombre de reads", 100, 10000, 1000, 100)
            read_length = st.number_input("Longueur (bp)", 50, 150, 75, 5)
        with col_b:
            error_rate = st.slider("Taux d'erreur", 0.0, 0.05, 0.01, 0.001, format="%.3f")

        if st.button("Générer le toy dataset", type="primary"):
            with st.spinner("Génération..."):
                fastq_content = generate_toy_dataset(TOY_REF, num_reads, read_length, error_rate)
                reads_raw = parse_fastq(fastq_content)
                st.session_state['reads_data'] = reads_raw
                st.session_state['file_type'] = 'fastq'
                st.session_state['ref_seq'] = TOY_REF
            st.success(f"{len(reads_raw)} reads générés.")
    else:
        uploaded = st.file_uploader("Fichier FASTQ ou FASTA", type=["fastq", "fq", "fasta", "fa", "txt"])
        if uploaded:
            content = uploaded.read().decode('utf-8', errors='ignore')
            if uploaded.name.lower().endswith(('.fastq', '.fq')):
                reads_data = parse_fastq(content)
                file_type = 'fastq'
            else:
                reads_data = [(rid, seq, None) for rid, seq in parse_fasta(content)]
                file_type = 'fasta'
            st.session_state['reads_data'] = reads_data
            st.session_state['file_type'] = file_type
            st.session_state['ref_seq'] = None
            st.success(f"{len(reads_data)} reads chargés ({file_type.upper()}).")

    if st.session_state.get('reads_data'):
        rd = st.session_state['reads_data']
        seqs = [r[1] for r in rd]
        lengths = [len(s) for s in seqs]
        c1, c2, c3 = st.columns(3)
        c1.metric("Reads", len(rd))
        c2.metric("Longueur moyenne (bp)", f"{np.mean(lengths):.1f}")
        c3.metric("Bases totales (bp)", f"{sum(lengths):,}")

        if st.session_state.get('file_type') == 'fastq':
            quals = [r[2] for r in rd if r[2]]
            if quals:
                q_mean = np.mean([
                    np.mean([ord(c) - 33 for c in q]) for q in quals
                ])
                st.metric("Score Phred moyen", f"Q{q_mean:.1f}")

with tab2:
    if not st.session_state.get('reads_data'):
        st.info("Chargez des données dans l'onglet Données.")
    else:
        seqs = [r[1] for r in st.session_state['reads_data']]
        with st.spinner("Comptage des k-mers..."):
            kmer_counts = count_kmers(seqs, k)
            stats = estimate_error_rate(kmer_counts)

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("k-mers distincts", f"{stats['distinct_kmers']:,}")
        c2.metric("k-mers solides", f"{sum(1 for v in kmer_counts.values() if v >= min_count):,}")
        c3.metric("Hapax", f"{stats['hapax_count']:,}")
        c4.metric("Taux d'erreur estimé", f"{stats['estimated_error_rate'] * 100:.2f} %")

        threshold = stats['solid_threshold']
        hist_data = [0] * (max_freq_display + 1)
        for v in kmer_counts.values():
            if v <= max_freq_display:
                hist_data[v] += 1

        fig, ax = plt.subplots(figsize=(11, 4), facecolor=CHART['bg'])
        ax.set_facecolor(CHART['face'])
        freqs = list(range(1, max_freq_display + 1))
        counts_plot = hist_data[1:max_freq_display + 1]
        colors = [CHART["error"] if f < threshold else LOT_COLOR[1] for f in freqs]
        ax.bar(freqs, counts_plot, color=colors, alpha=0.88, width=0.8, edgecolor='none')
        ax.axvline(x=threshold, color=CHART['warning'], linewidth=1.5, linestyle='--',
                   label=f'Seuil f = {threshold}')
        patch_err = mpatches.Patch(color=CHART['error'], label="Erreurs (faible fréq.)")
        patch_sol = mpatches.Patch(color=LOT_COLOR[1], label="k-mers solides")
        ax.legend(handles=[patch_err, patch_sol], facecolor=CHART['legend_bg'],
                  edgecolor=CHART['legend_edge'], labelcolor=CHART['legend_text'], fontsize=9)
        ax.set_xlabel(f"Fréquence du {k}-mer", color=CHART['text'])
        ax.set_ylabel("Nombre de k-mers", color=CHART['text'])
        ax.set_title(f"Histogramme — seuil f = {threshold}", color=CHART['title'])
        ax.tick_params(colors=CHART['text'])
        for spine in ('bottom', 'left'):
            ax.spines[spine].set_color(CHART['spine'])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', color=CHART['grid'], linewidth=0.5, alpha=0.7)
        plt.tight_layout()
        st.pyplot(fig)
        plt.close()

        st.session_state['kmer_counts'] = kmer_counts
        st.session_state['k'] = k
        st.session_state['min_count'] = min_count

with tab3:
    if not st.session_state.get('reads_data'):
        st.info("Chargez des données FASTQ dans l'onglet Données.")
    elif st.session_state.get('file_type') != 'fastq':
        st.warning("La conversion nécessite un fichier FASTQ.")
    else:
        q_threshold = st.slider("Seuil Phred", 10, 40, 20)
        rd = st.session_state['reads_data']
        fasta_lines, kept, discarded = [], 0, 0
        for rid, seq, qual in rd:
            if qual:
                avg = sum(ord(c) - 33 for c in qual) / len(qual)
                if avg >= q_threshold:
                    fasta_lines.extend([f">{rid}", seq])
                    kept += 1
                else:
                    discarded += 1
            else:
                fasta_lines.extend([f">{rid}", seq])
                kept += 1

        c1, c2 = st.columns(2)
        c1.metric("Reads conservés", kept)
        c2.metric("Reads filtrés", discarded)

        st.download_button(
            "Télécharger FASTA",
            '\n'.join(fasta_lines),
            "reads_filtered.fasta",
            "text/plain",
            type="primary",
        )
