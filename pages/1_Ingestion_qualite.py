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
from aff.ui import inject_styles, page_header, CHART

st.set_page_config(page_title="Lot 1 — Ingestion & Qualité", layout="wide")
inject_styles()

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
    "Lectures brutes, comptage des k-mers, estimation du taux d'erreur.",
)
st.divider()

with st.sidebar:
    st.markdown("### Paramètres")
    k = st.slider("Taille des k-mers (k)", min_value=7, max_value=31, value=21, step=2,
                   help="k impair recommandé pour éviter les palindromes")
    min_count = st.slider("Seuil 'solide' (min count)", min_value=1, max_value=10, value=2)
    max_freq_display = st.slider("Fréquence max affichée", 10, 100, 50)

    st.divider()
    st.markdown("### Théorie")
    st.markdown("""
    **k-mer** : sous-séquence de longueur k.

    **k-mers solides** : fréquence ≥ seuil (signal).

    **Hapax** : fréquence = 1 (bruit/erreur probable).

    *L'histogramme montre typiquement 2 pics :*
    - **Pic gauche** (faible fréq.) → erreurs
    - **Pic droit** (haute fréq.) → vraies séquences
    """)

tab1, tab2, tab3 = st.tabs(["Chargement des données", "Analyse k-mers", "Conversion FASTA"])

with tab1:
    st.markdown("### Source des données")

    data_source = st.radio(
        "Source",
        ["Toy Dataset (généré automatiquement)", "Uploader un fichier FASTQ/FASTA"],
        horizontal=True,
    )

    reads_data = []
    file_type = None

    if "Toy Dataset" in data_source:
        col_toy1, col_toy2 = st.columns(2)
        with col_toy1:
            num_reads = st.number_input("Nombre de reads", min_value=100, max_value=10000,
                                         value=1000, step=100)
            read_length = st.number_input("Longueur des reads", min_value=50, max_value=150,
                                           value=75, step=5)
        with col_toy2:
            error_rate = st.slider("Taux d'erreur simulé", 0.0, 0.05, 0.01, step=0.001,
                                    format="%.3f")
            st.metric("Séquence référence", f"{len(TOY_REF)} bp")

        if st.button("Générer le Toy Dataset", type="primary"):
            with st.spinner("Génération des reads simulés..."):
                fastq_content = generate_toy_dataset(TOY_REF, num_reads, read_length, error_rate)
                reads_raw = parse_fastq(fastq_content)
                reads_data = reads_raw
                file_type = 'fastq'
                st.session_state['reads_data'] = reads_data
                st.session_state['file_type'] = file_type
                st.session_state['ref_seq'] = TOY_REF

            st.success(f"{len(reads_data)} reads générés avec {error_rate * 100:.1f} % d'erreur")

            with st.expander("Aperçu du fichier FASTQ généré (10 premiers reads)"):
                preview_lines = fastq_content.split('\n')[:40]
                st.code('\n'.join(preview_lines), language="text")

    else:
        uploaded = st.file_uploader(
            "Déposer un fichier FASTQ ou FASTA",
            type=["fastq", "fq", "fasta", "fa", "txt"],
        )

        if uploaded:
            content = uploaded.read().decode('utf-8', errors='ignore')
            fname = uploaded.name.lower()

            if fname.endswith(('.fastq', '.fq')):
                reads_raw = parse_fastq(content)
                reads_data = reads_raw
                file_type = 'fastq'
            else:
                reads_fasta = parse_fasta(content)
                reads_data = [(rid, seq, None) for rid, seq in reads_fasta]
                file_type = 'fasta'

            st.session_state['reads_data'] = reads_data
            st.session_state['file_type'] = file_type
            st.session_state['ref_seq'] = None
            st.success(f"{len(reads_data)} reads chargés ({file_type.upper()})")

    if 'reads_data' in st.session_state and st.session_state['reads_data']:
        rd = st.session_state['reads_data']
        seqs = [r[1] for r in rd]

        st.divider()
        st.markdown("#### Statistiques de base")
        c1, c2, c3, c4 = st.columns(4)
        lengths = [len(s) for s in seqs]
        c1.metric("Reads", len(rd))
        c2.metric("Longueur moy.", f"{np.mean(lengths):.1f} bp")
        c3.metric("Min / Max", f"{min(lengths)} / {max(lengths)} bp")
        c4.metric("Bases totales", f"{sum(lengths):,} bp")

        if st.session_state['file_type'] == 'fastq':
            quals = [r[2] for r in rd if r[2]]
            if quals:
                avg_scores = []
                for q in quals:
                    scores = [ord(c) - 33 for c in q]
                    avg_scores.append(np.mean(scores))
                q_mean = np.mean(avg_scores)
                st.metric(
                    "Score Phred moyen", f"Q{q_mean:.1f}",
                    delta="Bonne qualité" if q_mean >= 20 else "Qualité faible",
                )

with tab2:
    st.markdown("### Analyse des k-mers")

    if 'reads_data' not in st.session_state or not st.session_state['reads_data']:
        st.info("Chargez ou générez des données dans l'onglet précédent.")
    else:
        rd = st.session_state['reads_data']
        seqs = [r[1] for r in rd]

        with st.spinner(f"Comptage des {k}-mers sur {len(seqs)} reads..."):
            kmer_counts = count_kmers(seqs, k)
            stats = estimate_error_rate(kmer_counts)

        col1, col2, col3, col4 = st.columns(4)
        col1.metric("k-mers distincts", f"{stats['distinct_kmers']:,}")
        col2.metric(
            "k-mers solides",
            f"{sum(1 for v in kmer_counts.values() if v >= min_count):,}",
            help=f"Apparaissant ≥ {min_count} fois",
        )
        col3.metric("Hapax (erreurs probables)", f"{stats['hapax_count']:,}")
        col4.metric("Taux d'erreur estimé", f"{stats['estimated_error_rate'] * 100:.2f} %")

        st.divider()
        st.markdown("#### Histogramme de fréquence des k-mers")

        hist_data = [0] * (max_freq_display + 1)
        for v in kmer_counts.values():
            if v <= max_freq_display:
                hist_data[v] += 1

        threshold = stats['solid_threshold']

        fig, ax = plt.subplots(figsize=(12, 5), facecolor=CHART['bg'])
        ax.set_facecolor(CHART['face'])

        freqs = list(range(1, max_freq_display + 1))
        counts_plot = hist_data[1:max_freq_display + 1]

        colors = [CHART['error'] if f < threshold else CHART['success'] for f in freqs]
        ax.bar(freqs, counts_plot, color=colors, alpha=0.85, width=0.8, edgecolor='none')

        ax.axvline(
            x=threshold, color=CHART['warning'], linewidth=2, linestyle='--',
            label=f'Seuil solide (freq = {threshold})',
        )

        ax.annotate(
            f'Seuil\nf = {threshold}',
            xy=(threshold, max(counts_plot) * 0.7),
            color=CHART['warning'], fontsize=9, ha='left',
            xytext=(threshold + 0.5, max(counts_plot) * 0.75),
        )

        patch_err = mpatches.Patch(color=CHART['error'], alpha=0.85, label="k-mers d'erreur (bruit)")
        patch_sol = mpatches.Patch(color=CHART['success'], alpha=0.85, label='k-mers solides (signal)')
        ax.legend(
            handles=[
                patch_err, patch_sol,
                plt.Line2D([0], [0], color=CHART['warning'], linestyle='--',
                           label=f'Seuil (f={threshold})'),
            ],
            facecolor=CHART['legend_bg'], edgecolor=CHART['legend_edge'],
            labelcolor=CHART['legend_text'], fontsize=9,
        )

        ax.set_xlabel(f'Fréquence du {k}-mer', color=CHART['text'], fontsize=11)
        ax.set_ylabel("Nombre de k-mers distincts", color=CHART['text'], fontsize=11)
        ax.set_title(f"Histogramme des {k}-mers — {len(seqs)} reads", color=CHART['title'], fontsize=13)
        ax.tick_params(colors=CHART['text'])
        ax.spines['bottom'].set_color(CHART['spine'])
        ax.spines['left'].set_color(CHART['spine'])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', color=CHART['grid'], linewidth=0.5, alpha=0.7)

        plt.tight_layout()
        st.pyplot(fig)
        plt.close()

        st.caption(
            f"**Interprétation** : le pic à faible fréquence (rouge) correspond aux k-mers issus "
            f"d'erreurs de séquençage. Le pic à plus haute fréquence (vert) représente les k-mers "
            f"vrais de la séquence cible. Le seuil automatique détecté est f = {threshold}."
        )

        with st.expander("Top 20 k-mers les plus fréquents"):
            top20 = kmer_counts.most_common(20)
            col_a, col_b = st.columns(2)
            for i, (km, cnt) in enumerate(top20):
                if i < 10:
                    col_a.code(f"{km}  ×{cnt}", language="text")
                else:
                    col_b.code(f"{km}  ×{cnt}", language="text")

        st.session_state['kmer_counts'] = kmer_counts
        st.session_state['k'] = k
        st.session_state['min_count'] = min_count

with tab3:
    st.markdown("### Conversion FASTQ → FASTA")
    st.markdown(
        "Filtre les reads selon un seuil de qualité Phred (Q ≥ 20 par défaut) "
        "et exporte au format FASTA."
    )

    if 'reads_data' not in st.session_state or not st.session_state['reads_data']:
        st.info("Chargez des données FASTQ dans l'onglet 1.")
    elif st.session_state.get('file_type') != 'fastq':
        st.warning("La conversion FASTQ → FASTA nécessite un fichier d'entrée au format FASTQ.")
    else:
        rd = st.session_state['reads_data']
        q_threshold = st.slider("Seuil de qualité Phred", min_value=10, max_value=40, value=20)

        fasta_lines = []
        kept, discarded = 0, 0
        for rid, seq, qual in rd:
            if qual:
                scores = [ord(c) - 33 for c in qual]
                avg = sum(scores) / len(scores)
                if avg >= q_threshold:
                    fasta_lines.append(f">{rid}")
                    fasta_lines.append(seq)
                    kept += 1
                else:
                    discarded += 1
            else:
                fasta_lines.append(f">{rid}")
                fasta_lines.append(seq)
                kept += 1

        fasta_output = '\n'.join(fasta_lines)

        col1, col2 = st.columns(2)
        col1.metric("Reads conservés", kept, delta=f"-{discarded} filtrés")
        col2.metric(
            "Taux de conservation",
            f"{100 * kept / (kept + discarded):.1f} %" if (kept + discarded) > 0 else "–",
        )

        st.download_button(
            label="Télécharger le fichier FASTA filtré",
            data=fasta_output,
            file_name="reads_filtered.fasta",
            mime="text/plain",
            type="primary",
        )

        with st.expander("Aperçu FASTA (20 premières entrées)"):
            preview = '\n'.join(fasta_lines[:40])
            st.code(preview, language="text")
