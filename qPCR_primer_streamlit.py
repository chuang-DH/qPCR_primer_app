import streamlit as st
import math
import time

# ====================================================
# 序列處理與檢查函式
# ====================================================
def sanitize_seq(seq):
    if seq is None:
        return ""
    valid = {"A", "T", "C", "G"}
    return ''.join([ch for ch in seq.upper() if ch in valid])

def calculate_tm(seq):
    seq = sanitize_seq(seq)
    length = len(seq)
    if length == 0:
        return 0.0
    A = seq.count('A')
    T = seq.count('T')
    G = seq.count('G')
    C = seq.count('C')
    if length < 14:
        tm = 4 * (G + C) + 2 * (A + T)
    else:
        tm = 64.9 + 41 * (G + C - 16.4) / (A + T + G + C)
    return round(tm, 2)

def calculate_gc_content(seq):
    seq = sanitize_seq(seq)
    if len(seq) == 0:
        return 0.0
    gc = seq.count('G') + seq.count('C')
    return round(gc / len(seq) * 100, 2)

def reverse_complement(seq):
    seq = sanitize_seq(seq)
    if not seq:
        return ""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp.get(base, base) for base in reversed(seq))

def is_self_complementary(seq, min_match=4):
    seq = sanitize_seq(seq)
    length = len(seq)
    if length < min_match + 1:
        return False
    rev_comp = reverse_complement(seq)
    for shift in range(1, length):
        match = 0
        for i in range(length - shift):
            if seq[i] == rev_comp[i + shift]:
                match += 1
                if match >= min_match:
                    return True
            else:
                match = 0
    return False

def gc_rich_3prime(seq, threshold=4):
    seq = sanitize_seq(seq)
    if len(seq) == 0:
        return False
    tail = seq[-5:]
    count_gc = tail.count('G') + tail.count('C')
    return count_gc >= threshold

def check_hairpin(seq, min_stem=4, loop_min=3, loop_max=8):
    seq = sanitize_seq(seq)
    length = len(seq)
    if length < (2 * min_stem + loop_min):
        return False
    for i in range(length - min_stem - loop_min):
        for loop_size in range(loop_min, loop_max + 1):
            if i + min_stem + loop_size >= length:
                continue
            stem1 = seq[i:i + min_stem]
            stem2 = seq[i + min_stem + loop_size:i + 2 * min_stem + loop_size]
            if stem2 == reverse_complement(stem1):
                return True
    return False

# ====================================================
# Primer candidate generation & pairing
# ====================================================
def find_primer_candidates(seq, min_len, max_len, gc_min, gc_max, tm_target, tm_tol):
    seq = sanitize_seq(seq)
    n = len(seq)
    candidates = []
    for L in range(min_len, max_len + 1):
        if L > n:
            continue
        for i in range(0, n - L + 1):
            s = seq[i:i+L]
            gc = calculate_gc_content(s)
            if gc < gc_min or gc > gc_max:
                continue
            tm = calculate_tm(s)
            if abs(tm - tm_target) > tm_tol:
                continue
            warnings = []
            if is_self_complementary(s):
                warnings.append("self-dimer")
            if gc_rich_3prime(s):
                warnings.append("3'-GC-rich")
            if check_hairpin(s):
                warnings.append("hairpin")
            candidates.append({
                "pos": i,
                "seq": s,
                "len": L,
                "tm": tm,
                "gc": gc,
                "warnings": warnings
            })
    return candidates

def design_qpcr_primers(seq,
                        primer_min=18, primer_max=22,
                        gc_min=40, gc_max=60,
                        tm_target=60.0, tm_tol=6.0,
                        amp_min=70, amp_max=200,
                        top_n=10,
                        max_candidates_per_side=8000):
    t0 = time.time()
    seq_clean = sanitize_seq(seq)
    n = len(seq_clean)
    if n == 0:
        return []
    if n > 20000:
        seq_clean = seq_clean[:20000]
        n = len(seq_clean)

    f_candidates = find_primer_candidates(seq_clean, primer_min, primer_max, gc_min, gc_max, tm_target, tm_tol)
    r_raw = find_primer_candidates(seq_clean, primer_min, primer_max, gc_min, gc_max, tm_target, tm_tol)

    r_by_pos = {}
    for rc in r_raw:
        pos = rc["pos"]
        r_by_pos.setdefault(pos, []).append(rc)

    def candidate_score_for_prune(c):
        return abs(c["gc"] - 50) + abs(c["tm"] - tm_target)
    f_candidates.sort(key=candidate_score_for_prune)
    r_sorted_all = sorted(r_raw, key=candidate_score_for_prune)

    if len(f_candidates) > max_candidates_per_side:
        f_candidates = f_candidates[:max_candidates_per_side]
    if len(r_sorted_all) > max_candidates_per_side:
        r_sorted_all = r_sorted_all[:max_candidates_per_side]
        r_by_pos = {}
        for rc in r_sorted_all:
            r_by_pos.setdefault(rc["pos"], []).append(rc)

    r_positions = sorted(r_by_pos.keys())
    pairs = []

    for f in f_candidates:
        i = f["pos"]
        for r_pos in r_positions:
            if r_pos <= i:
                continue
            rlist = r_by_pos.get(r_pos, [])
            for rc in rlist:
                amp_len = (rc["pos"] + rc["len"]) - i
                if amp_len < amp_min or amp_len > amp_max:
                    continue
                rev_primer_seq = reverse_complement(seq_clean[rc["pos"]:rc["pos"]+rc["len"]])
                score = 0.0
                score += abs(f["gc"] - 50) + abs(rc["gc"] - 50)
                score += 1.5 * (abs(f["tm"] - tm_target) + abs(rc["tm"] - tm_target))
                score += (len(f["warnings"]) + len(rc["warnings"])) * 20.0
                if gc_rich_3prime(f["seq"]):
                    score += 8.0
                if gc_rich_3prime(rev_primer_seq):
                    score += 8.0
                mid = (amp_min + amp_max) / 2.0
                score += abs(amp_len - mid) / 5.0
                pairs.append({
                    "forward": f,
                    "reverse": {
                        "pos": rc["pos"],
                        "seq": rev_primer_seq,
                        "len": rc["len"],
                        "tm": rc["tm"],
                        "gc": rc["gc"],
                        "warnings": rc["warnings"]
                    },
                    "amp_len": amp_len,
                    "score": score
                })
    pairs.sort(key=lambda x: x["score"])
    duration = time.time() - t0
    return pairs[:top_n]

# ====================================================
# Streamlit Web GUI
# ====================================================
st.set_page_config(page_title="qPCR Primer Designer", layout="wide")
st.title("🧬 qPCR Primer 設計工具")
st.markdown("輸入DNA序列，系統會自動搜尋符合條件的primer pair。")

seq_input = st.text_area("請輸入DNA序列（最多15000 bp）", height=200)
col1, col2, col3, col4, col5, col6 = st.columns(6)
with col1:
    primer_min = st.number_input("Primer min", 10, 30, 18)
with col2:
    primer_max = st.number_input("Primer max", 10, 40, 22)
with col3:
    gc_min = st.number_input("GC min %", 30, 70, 40)
with col4:
    gc_max = st.number_input("GC max %", 30, 70, 60)
with col5:
    tm_target = st.number_input("Tm target °C", 40.0, 70.0, 60.0)
with col6:
    tm_tol = st.number_input("Tm tolerance", 1.0, 10.0, 6.0)

col7, col8 = st.columns(2)
with col7:
    amp_min = st.number_input("Amplicon min", 50, 500, 70)
with col8:
    amp_max = st.number_input("Amplicon max", 50, 500, 200)

top_n = st.slider("顯示結果筆數", 5, 100, 20)

if st.button("🚀 開始設計Primer"):
    if not seq_input:
        st.warning("請先輸入DNA序列！")
    else:
        with st.spinner("正在設計primer，請稍候..."):
            pairs = design_qpcr_primers(seq_input,
                                        primer_min, primer_max,
                                        gc_min, gc_max,
                                        tm_target, tm_tol,
                                        amp_min, amp_max,
                                        top_n)
        if not pairs:
            st.error("沒有找到符合條件的primer對。請嘗試放寬GC/Tm或amplicon範圍。")
        else:
            data = []
            for idx, p in enumerate(pairs, start=1):
                f = p["forward"]; r = p["reverse"]
                warn = ";".join(f["warnings"] + r["warnings"])
                data.append({
                    "#": idx,
                    "F_seq": f["seq"],
                    "F_pos": f["pos"]+1,
                    "F_len": f["len"],
                    "F_Tm": f["tm"],
                    "F_GC": f["gc"],
                    "R_seq": r["seq"],
                    "R_pos": r["pos"]+1,
                    "R_len": r["len"],
                    "R_Tm": r["tm"],
                    "R_GC": r["gc"],
                    "Amp_len": p["amp_len"],
                    "Score": round(p["score"], 2),
                    "Warnings": warn
                })
            st.success(f"找到 {len(data)} 組候選primer對！")
            st.dataframe(data, use_container_width=True)
