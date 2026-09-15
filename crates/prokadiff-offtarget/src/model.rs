use std::fmt;
use std::str::FromStr;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    pub fn as_char(self) -> char {
        match self {
            Self::Plus => '+',
            Self::Minus => '-',
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_char())
    }
}

impl FromStr for Strand {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim() {
            "+" | "plus" | "1" | "F" | "FORWARD" => Ok(Self::Plus),
            "-" | "minus" | "-1" | "R" | "REVERSE" => Ok(Self::Minus),
            other => Err(format!("unknown strand: {other}")),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BulgeType {
    None,
    Dna,
    Rna,
}

impl fmt::Display for BulgeType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::None => write!(f, "none"),
            Self::Dna => write!(f, "DNA"),
            Self::Rna => write!(f, "RNA"),
        }
    }
}

impl FromStr for BulgeType {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "" | "none" | "." | "x" => Ok(Self::None),
            "dna" | "d" => Ok(Self::Dna),
            "rna" | "r" => Ok(Self::Rna),
            other => Err(format!("unknown bulge type: {other}")),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PamSide {
    ThreePrime,
    FivePrime,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NucleaseProfile {
    pub name: String,
    pub spacer_len: usize,
    pub pam_pattern: String,
    pub pam_side: PamSide,
}

impl NucleaseProfile {
    pub fn spcas9() -> Self {
        Self {
            name: "SpCas9".to_string(),
            spacer_len: 20,
            pam_pattern: "NGG".to_string(),
            pam_side: PamSide::ThreePrime,
        }
    }

    pub fn cas12a() -> Self {
        Self {
            name: "Cas12a".to_string(),
            spacer_len: 20,
            pam_pattern: "TTTV".to_string(),
            pam_side: PamSide::FivePrime,
        }
    }

    pub fn custom(
        name: impl Into<String>,
        spacer_len: usize,
        pam_pattern: impl Into<String>,
        pam_side: PamSide,
    ) -> Self {
        Self {
            name: name.into(),
            spacer_len,
            pam_pattern: pam_pattern.into(),
            pam_side,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct OffTargetSite {
    pub site_id: String,
    pub seq_id: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,

    pub guide: String,
    pub target_seq: String,
    pub pam: String,

    pub mismatches: u32,
    pub bulge_type: BulgeType,
    pub bulge_size: u32,

    pub search_backend: String,

    pub cfd_score: Option<f64>,
    pub hsu_score: Option<f64>,
}

impl OffTargetSite {
    pub fn tsv_header() -> &'static str {
        "site_id\tseq_id\tstart\tend\tstrand\tguide\ttarget_seq\tpam\tmismatches\tbulge_type\tbulge_size\tsearch_backend\tcfd_score\thsu_score"
    }

    pub fn to_tsv_row(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.site_id,
            self.seq_id,
            self.start,
            self.end,
            self.strand,
            self.guide,
            self.target_seq,
            self.pam,
            self.mismatches,
            self.bulge_type,
            self.bulge_size,
            self.search_backend,
            self.cfd_score
                .map(|s| format!("{s:.4}"))
                .unwrap_or_else(|| "NA".to_string()),
            self.hsu_score
                .map(|s| format!("{s:.4}"))
                .unwrap_or_else(|| "NA".to_string())
        )
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MutationOffTargetLink {
    pub mutation_id: String,
    pub site_id: String,
    pub mutation_type: String,
    pub mutation_position: u64,
    pub site_start: u64,
    pub site_end: u64,
    pub distance_to_site: u64,
    pub mismatches: u32,
    pub pam: String,
    pub cfd_score: Option<f64>,
    pub association_window: u64,
}

impl MutationOffTargetLink {
    pub fn tsv_header() -> &'static str {
        "mutation_id\tsite_id\tmutation_type\tmutation_position\tsite_start\tsite_end\tdistance_to_site\tmismatches\tpam\tcfd_score\tassociation_window"
    }

    pub fn to_tsv_row(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.mutation_id,
            self.site_id,
            self.mutation_type,
            self.mutation_position,
            self.site_start,
            self.site_end,
            self.distance_to_site,
            self.mismatches,
            self.pam,
            self.cfd_score
                .map(|s| format!("{s:.4}"))
                .unwrap_or_else(|| "NA".to_string()),
            self.association_window
        )
    }
}
