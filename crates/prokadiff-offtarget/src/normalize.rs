pub fn revcomp_dna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'U' => b'A',
            other => other,
        })
        .collect()
}

pub fn revcomp_iupac(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' | b'U' => b'A',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'D' => b'H',
            b'H' => b'D',
            b'V' => b'B',
            b'N' => b'N',
            other => other,
        })
        .collect()
}

pub fn iupac_match(pattern: u8, target: u8) -> bool {
    let p = pattern.to_ascii_uppercase();
    let t = match target.to_ascii_uppercase() {
        b'U' => b'T',
        other => other,
    };
    match p {
        b'A' => t == b'A',
        b'C' => t == b'C',
        b'G' => t == b'G',
        b'T' | b'U' => t == b'T',
        b'R' => t == b'A' || t == b'G',
        b'Y' => t == b'C' || t == b'T',
        b'S' => t == b'G' || t == b'C',
        b'W' => t == b'A' || t == b'T',
        b'K' => t == b'G' || t == b'T',
        b'M' => t == b'A' || t == b'C',
        b'B' => t == b'C' || t == b'G' || t == b'T',
        b'D' => t == b'A' || t == b'G' || t == b'T',
        b'H' => t == b'A' || t == b'C' || t == b'T',
        b'V' => t == b'A' || t == b'C' || t == b'G',
        b'N' => t == b'A' || t == b'C' || t == b'G' || t == b'T',
        _ => p == t,
    }
}

pub fn iupac_match_slice(pattern: &[u8], target: &[u8]) -> bool {
    if pattern.len() != target.len() {
        return false;
    }
    pattern
        .iter()
        .zip(target.iter())
        .all(|(&p, &t)| iupac_match(p, t))
}

pub fn hamming(a: &[u8], b: &[u8]) -> u32 {
    a.iter()
        .zip(b.iter())
        .filter(|(&x, &y)| {
            let x_up = if x.eq_ignore_ascii_case(&b'U') {
                b'T'
            } else {
                x.to_ascii_uppercase()
            };
            let y_up = if y.eq_ignore_ascii_case(&b'U') {
                b'T'
            } else {
                y.to_ascii_uppercase()
            };
            x_up != y_up
        })
        .count() as u32
}

pub fn clean_nucleotide_string(s: &str) -> String {
    s.trim()
        .chars()
        .filter(|c| !c.is_whitespace() && *c != '-')
        .map(|c| c.to_ascii_uppercase())
        .collect()
}
