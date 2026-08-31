//! Canonical 3′ (right) alignment of RA short indels.
//!
//! Bowtie2 CIGAR placement in tandem repeats is not unique: `INS 1999 TA` and
//! `INS 2000 AT` are the same haplotype. Genome Diff / `gdtools SUBTRACT` keys
//! are exact, so we slide INS/DEL to the 3′-most equivalent coordinate — the
//! placement breseq writes on the mutation line (RA evidence may stay 5′).

/// Insert `ins` after 1-based reference position `pos`. Slide 3′ while the
/// first inserted base equals the next reference base, rotating the oligo.
pub fn right_align_ins(ref_seq: &[u8], pos: u64, ins: &[u8]) -> (u64, Vec<u8>) {
    if ins.is_empty() || pos == 0 {
        return (pos, ins.to_vec());
    }
    let mut p = pos;
    let mut s: Vec<u8> = ins.iter().map(|b| b.to_ascii_uppercase()).collect();
    while (p as usize) < ref_seq.len() {
        let next = ref_seq[p as usize].to_ascii_uppercase();
        if s[0] != next {
            break;
        }
        s.rotate_left(1);
        p += 1;
    }
    (p, s)
}

/// Deletion of `size` bases starting at 1-based `pos`. Slide 3′ while the
/// first deleted base equals the base immediately after the deleted span.
pub fn right_align_del(ref_seq: &[u8], pos: u64, size: u64) -> u64 {
    if pos == 0 || size == 0 {
        return pos;
    }
    let n = size as usize;
    let mut p = pos;
    loop {
        let first = p as usize - 1;
        let after = first.saturating_add(n);
        if after >= ref_seq.len() {
            break;
        }
        if !ref_seq[first].eq_ignore_ascii_case(&ref_seq[after]) {
            break;
        }
        p += 1;
    }
    p
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Place `window` so that its first base is 1-based `start_1`.
    fn pad_at(start_1: u64, window: &[u8]) -> Vec<u8> {
        let start_0 = start_1 as usize - 1;
        let mut s = vec![b'N'; start_0 + window.len()];
        s[start_0..start_0 + window.len()].copy_from_slice(window);
        s
    }

    #[test]
    fn synth_ins_1999_ta_matches_breseq_2000_at() {
        // Seed-42 `generate.sh` window 1-based 1989–2013: …G T C…
        let ref_seq = pad_at(1989, b"CTTAGTGTTCGTCTCCGCTATTCTC");
        assert_eq!(ref_seq[1998], b'G');
        assert_eq!(ref_seq[1999], b'T');
        assert_eq!(ref_seq[2000], b'C');
        let (p, s) = right_align_ins(&ref_seq, 1999, b"TA");
        assert_eq!((p, s), (2000, b"AT".to_vec()));
        let (p, s) = right_align_ins(&ref_seq, 2000, b"AT");
        assert_eq!((p, s), (2000, b"AT".to_vec()));
    }

    #[test]
    fn synth_del_5600_shifts_to_breseq_5601() {
        // Seed-42 window 1-based 5589–5613: …A T A A…
        let ref_seq = pad_at(5589, b"GCTGGTAAACCATAACTGTCGCAGC");
        assert_eq!(&ref_seq[5599..5603], b"ATAA");
        assert_eq!(right_align_del(&ref_seq, 5600, 2), 5601);
        assert_eq!(right_align_del(&ref_seq, 5601, 2), 5601);
    }

    #[test]
    fn ins_already_rightmost_unchanged() {
        let mut ref_seq = vec![b'G'; 20];
        ref_seq[9] = b'C';
        let (p, s) = right_align_ins(&ref_seq, 10, b"AT");
        assert_eq!((p, s.as_slice()), (10, &b"AT"[..]));
    }

    #[test]
    fn ins_a_slides_to_end_of_homopolymer() {
        // 1-based 1–3 G, 4–7 A, 8–9 C.
        let ref_seq = b"GGGAAAACC";
        let (p, s) = right_align_ins(ref_seq, 3, b"A");
        assert_eq!((p, s), (7, b"A".to_vec()));
    }

    #[test]
    fn del_shifts_through_matching_flank() {
        // 1-based 5–8 = A T A T; DEL of AT at 5 ≡ DEL at 7.
        let ref_seq = b"GGGGATATGG";
        assert_eq!(right_align_del(ref_seq, 5, 2), 7);
    }

    #[test]
    fn del_one_bp_slides_to_end_of_homopolymer() {
        let ref_seq = b"GGGAAAACC";
        assert_eq!(right_align_del(ref_seq, 4, 1), 7);
    }
}
