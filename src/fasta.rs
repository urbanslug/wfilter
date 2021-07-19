use seq_io::fasta::Reader;
use std::io::Read;

pub struct Fasta {
    pub header: Vec<u8>,
    pub seq: Vec<u8>
}

pub type FastaFile = Vec<Fasta>;

// TODO: have the from methods return Fasta and not FastaFile
impl Fasta {
    fn parse_fasta<R>(reader: Reader<R>) -> FastaFile
    where R: Read
    {
        reader.into_records().map(|r| {
            let r = r.unwrap();
            Fasta {
                header: r.head.clone(),
                seq: r.seq.clone()
            }
        }).collect()
    }

    // TODO: handle the header
    pub fn from_str(data: &str) -> FastaFile {
        let bytestring = &data.as_bytes()[..];
        let reader = Reader::new(bytestring);

        Self::parse_fasta(reader)
    }

    pub fn from_path(fp: &str) -> FastaFile {
        let reader = Reader::from_path(fp).unwrap();
        Self::parse_fasta(reader)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static SPECIES_X: &str = ">species_x\n\
                         TCTATACTGCGCGTTTATCTAGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\
                         TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT";

    // with newlines
    static SPECIES_X_NL: &str = ">species_x\n\
                              TCTATACTGCGCGTTTATCTAGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\n\
                              TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT";

    // with newlines
    static MULTI_FASTA: &str = ">species_x\n\
                         TCTATACTGCGCGTTTATCTAGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\n\
                         TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\n\
                         >species_y\n\
                         ATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAACTAAAATAGT\n\
                         TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAACTAAAATAG\n";
    #[test]
    fn test_with_newlines() {
        let multi = Fasta::from_str(MULTI_FASTA);
        let with_newlines = Fasta::from_str(SPECIES_X_NL);

        let foo = |f: FastaFile| -> Vec<String> {
            f.into_iter()
                .map(|x: Fasta| {
                    String::from_utf8(x.seq).unwrap()
                })
                .collect()
        };

        assert_eq!(vec!["TCTATACTGCGCGTTTATCTAGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\
                         TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT"],
                   foo(with_newlines));
        assert_eq!(vec!["TCTATACTGCGCGTTTATCTAGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\
                         TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT",
                        "ATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAACTAAAATAGT\
                         TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAACTAAAATAG"],
                   foo(multi));
    }

    #[test]
    fn test_without_newlines() {
        let species_x = Fasta::from_str(SPECIES_X);

        let foo = |f: FastaFile| -> Vec<String> {
            f.into_iter()
                .map(|x: Fasta| {
                    String::from_utf8(x.seq).unwrap()
                })
                .collect()
        };

        assert_eq!(vec!["TCTATACTGCGCGTTTATCTAGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT\
                         TCTATACTGCGCGTTTGGAGAAATAACTATCAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT"],
                   foo(species_x));
    }

    #[test]
    fn test_handle_headers() {
        let multi = Fasta::from_str(MULTI_FASTA);
        let species_x = Fasta::from_str(SPECIES_X);
        let with_newlines = Fasta::from_str(SPECIES_X_NL);

        let foo = |f: FastaFile| -> Vec<String> {
            f.into_iter()
                .map(|x: Fasta| {
                    String::from_utf8(x.header).unwrap()
                })
                .collect()
        };

        assert_eq!(vec!["species_x", "species_y"], foo(multi));
        assert_eq!(vec!["species_x"], foo(species_x));
        assert_eq!(vec!["species_x"], foo(with_newlines));
    }
}
