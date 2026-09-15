use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::model::{MutationOffTargetLink, OffTargetSite};

/// Writes offtarget_sites.tsv to destination path.
pub fn write_offtarget_sites_tsv(
    sites: &[OffTargetSite],
    path: impl AsRef<Path>,
) -> std::io::Result<()> {
    let mut file = BufWriter::new(File::create(path)?);
    writeln!(file, "{}", OffTargetSite::tsv_header())?;
    for site in sites {
        writeln!(file, "{}", site.to_tsv_row())?;
    }
    file.flush()?;
    Ok(())
}

/// Writes mutation_offtarget_links.tsv to destination path.
pub fn write_mutation_offtarget_links_tsv(
    links: &[MutationOffTargetLink],
    path: impl AsRef<Path>,
) -> std::io::Result<()> {
    let mut file = BufWriter::new(File::create(path)?);
    writeln!(file, "{}", MutationOffTargetLink::tsv_header())?;
    for link in links {
        writeln!(file, "{}", link.to_tsv_row())?;
    }
    file.flush()?;
    Ok(())
}
