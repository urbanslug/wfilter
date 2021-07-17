use std::fs;
use seq_io::fasta::{Reader};

pub fn read_file(fp: &str) -> Vec<String> {
    let paf_byte_vector: Vec<u8> = fs::read(fp).expect("Error reading file");
    let x: String = String::from_utf8_lossy(&paf_byte_vector)
        .parse()
        .expect("Error decoding file. Does it contain non ASCII characters?");

    x.lines().map(|x| x.to_string()).collect::<Vec<String>>()
}

pub fn read_fasta(fp: &str) -> String {
    let mut reader = Reader::from_path(fp).unwrap();
    let mut x =  String::new();

    while let Some(record) = reader.next() {
        let record = record.unwrap();

        record.seq_lines().for_each(|s| {
            x.push_str(std::str::from_utf8(s).unwrap());
        });
    }

    // println!("{}", x);
    x
}
