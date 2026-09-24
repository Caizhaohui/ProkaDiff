use std::str::FromStr;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SchemaVersion {
    V1,
    #[default]
    V2,
}

impl SchemaVersion {
    pub const fn is_v2(self) -> bool {
        matches!(self, Self::V2)
    }
}

impl FromStr for SchemaVersion {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "v1" | "1" => Ok(Self::V1),
            "v2" | "2" => Ok(Self::V2),
            _ => Err("schema version must be v1 or v2".to_string()),
        }
    }
}
