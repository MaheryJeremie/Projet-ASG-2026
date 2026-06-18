"""
Lot 3 : Assemblage par filtre de Bloom et graphe de de Bruijn implicite
"""
import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from modules.assembly import build_bloom_filter, assemble, compare_memory
from modules.kmer import generate_toy_dataset, parse_fastq, extract_kmers
from aff.ui import inject_styles, page_header, sidebar_brand, sidebar_section, lot_chart, LOT_COLOR

st.set_page_config(page_title="Lot 3 — Assemblage Bloom", layout="wide")
inject_styles(lot=3)
CHART = lot_chart(3)

TOY_REF = (
    "ATGGAGTTCAACACCAATGTCACAGGCGAATCGTACAACCAGCACATCATTGTGGCTGTGATCAACATCATCAAC"
    "GGTTTGAATGTGGACTTCAAGGGCGAGTACTTGTACACGGCGAAGATCCCTGCGGTCGAGTTCCGCAAGCAGGTC"
    "GAGTTCAAGTACGCCAACAAGGTCACCGTCACGTACACCAACAGCGGCGTCGTCATCGCCGGCGAGATCATCAAG"
    "GCCATCAAGGAGACGCTGAACGCCGCGGTGAACAACATCGTCGACGGCAAGATCATCGCCAACGACTTCACCAAG"
    "TTCGAGCGCAAGGCGATCAAGGAGTTCGGCGAGAACAAGGTCGACATCGCCAAGCAGTTCGTCATCGACGAGAAC"
)

page_header(
    "Lot 3 — Assemblage",
    "Filtre de Bloom, graphe de de Bruijn implicite, production de contigs.",
)

with st.sidebar:
    sidebar_brand(lot=3)
    sidebar_section("Paramètres")
    k_assembly = st.slider("k (taille des k-mers)", 11, 31, 21, step=2)
    min_count = st.slider("min_count (seuil k-mer solide)", 1, 10, 2)
    fpr_target = st.select_slider(
        "FPR (taux de faux positifs)",
        options=[0.001, 0.005, 0.01, 0.05, 0.1],
        value=0.01,
        format_func=lambda x: f"{x * 100:.1f} %",
    )

    sidebar_section("Rappels")
    st.markdown("""
    **FPR** : Taux de Faux Positifs (filtre de Bloom).

    **N50** : longueur du contig médian.

    **bp** : paire de bases.

    Le graphe de de Bruijn n'est jamais stocké : on teste
    chaque extension via le filtre de Bloom.
    """)

tab1, tab2, tab3 = st.tabs(["Filtre de Bloom", "Assemblage", "Mémoire"])

with tab1:
    data_source = st.radio("Source", ["Toy dataset", "Données du lot 1"], horizontal=True)

    if "Toy" in data_source:
        n_reads = st.number_input("Nombre de reads", 100, 5000, 1000, 100)
        read_len = st.number_input("Longueur (bp)", 50, 150, 75, 5)
        err_rate = st.slider("Taux d'erreur", 0.0, 0.05, 0.01, 0.001)
        if st.button("Construire le filtre", type="primary"):
            with st.spinner("Construction..."):
                fq = generate_toy_dataset(TOY_REF, n_reads, read_len, err_rate)
                st.session_state['assembly_seqs'] = [r[1] for r in parse_fastq(fq)]
                st.session_state['ref_seq'] = TOY_REF
    elif st.session_state.get('reads_data'):
        st.session_state['assembly_seqs'] = [r[1] for r in st.session_state['reads_data']]
        st.info(f"{len(st.session_state['assembly_seqs'])} reads du lot 1.")
    else:
        st.warning("Aucune donnée. Utilisez le toy dataset.")

    if st.session_state.get('assembly_seqs'):
        seqs = st.session_state['assembly_seqs']
        bloom, kmer_counts = build_bloom_filter(seqs, k_assembly, min_count, fpr_target)
        st.session_state.update({
            'bloom': bloom, 'kmer_counts_asm': kmer_counts, 'assembly_k': k_assembly,
        })
        b = bloom.stats()
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("m (bits)", f"{b['m_bits']:,}")
        c2.metric("k_hash (fonctions de hachage)", b['k_hash_functions'])
        c3.metric("Mémoire (Ko)", f"{b['memory_kb']:.1f}")
        c4.metric("FPR estimé", f"{b['estimated_fpr'] * 100:.3f} %")

with tab2:
    if 'bloom' not in st.session_state:
        st.info("Construisez le filtre de Bloom d'abord.")
    else:
        max_contigs = st.number_input("Nombre max de contigs", 1, 50, 10)
        if st.button("Lancer l'assemblage", type="primary"):
            with st.spinner("Assemblage..."):
                res = assemble(
                    st.session_state['assembly_seqs'],
                    st.session_state['assembly_k'],
                    min_count, fpr_target, max_contigs,
                )
                st.session_state['assembly_result'] = res

        if st.session_state.get('assembly_result'):
            res = st.session_state['assembly_result']
            contigs = res['contigs']

            c1, c2, c3 = st.columns(3)
            c1.metric("Contigs", res['num_contigs'])
            c2.metric("N50 (bp)", f"{res['n50']:,}")
            c3.metric("Bases assemblées (bp)", f"{res['total_bases_assembled']:,}")

            if st.session_state.get('ref_seq') and contigs:
                ref = st.session_state['ref_seq']
                k_asm = st.session_state['assembly_k']
                ck = set(extract_kmers(contigs[0], k_asm))
                rk = set(extract_kmers(ref, k_asm))
                if ck:
                    identity = len(ck & rk) / len(ck)
                    st.metric("Identité vs référence", f"{identity * 100:.1f} %",
                              help="Critère du toy dataset : ≥ 98 %")

            for i, contig in enumerate(contigs):
                st.markdown(f"""
                <div class="contig-block">
                <b>Contig {i + 1}</b> — {len(contig)} bp<br/>
                {contig[:100]}{'…' if len(contig) > 100 else ''}
                </div>
                """, unsafe_allow_html=True)

            if contigs:
                fig, ax = plt.subplots(figsize=(7, 3), facecolor=CHART['bg'])
                ax.set_facecolor(CHART['face'])
                lengths = [len(c) for c in contigs]
                ax.bar(range(1, len(lengths) + 1), lengths, color=LOT_COLOR[3])
                ax.axhline(y=res['n50'], color=CHART['warning'], linestyle='--', label=f'N50 = {res["n50"]}')
                ax.set_xlabel("Contig", color=CHART['text'])
                ax.set_ylabel("Longueur (bp)", color=CHART['text'])
                ax.legend(facecolor=CHART['legend_bg'], edgecolor=CHART['legend_edge'],
                          labelcolor=CHART['legend_text'])
                ax.tick_params(colors=CHART['text'])
                ax.spines[:].set_color(CHART['spine'])
                plt.tight_layout()
                st.pyplot(fig)
                plt.close()

            fasta_out = '\n'.join(f">contig_{i + 1}\n{c}" for i, c in enumerate(contigs))
            st.download_button("Exporter FASTA", fasta_out, "contigs.fasta", "text/plain")

with tab3:
    if 'bloom' not in st.session_state:
        st.info("Construisez le filtre de Bloom d'abord.")
    else:
        mem = compare_memory(st.session_state['kmer_counts_asm'], st.session_state['bloom'])
        c1, c2, c3 = st.columns(3)
        c1.metric("k-mers distincts", f"{mem['n_distinct_kmers']:,}")
        c2.metric("Dictionnaire Python (Mo)", f"{mem['dict_memory_mb']:.2f}")
        c3.metric("Filtre de Bloom (Mo)", f"{mem['bloom_memory_mb']:.2f}",
                  delta=f"÷ {mem['compression_ratio']:.0f}")

        fig, ax = plt.subplots(figsize=(5, 3), facecolor=CHART['bg'])
        ax.set_facecolor(CHART['face'])
        vals = [mem['dict_memory_mb'], mem['bloom_memory_mb']]
        ax.bar(['Dictionnaire', 'Bloom'], vals, color=[CHART['error'], LOT_COLOR[3]], width=0.45)
        ax.set_ylabel("Mémoire (Mo)", color=CHART['text'])
        ax.tick_params(colors=CHART['text'])
        ax.spines[:].set_color(CHART['spine'])
        plt.tight_layout()
        st.pyplot(fig)
        plt.close()
