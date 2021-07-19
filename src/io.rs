use std::fs;

pub fn read_file(fp: &str) -> Vec<String> {
    let paf_byte_vector: Vec<u8> = fs::read(fp).expect("Error reading file");
    let x: String = String::from_utf8_lossy(&paf_byte_vector)
        .parse()
        .expect("Error decoding file. Does it contain non ASCII characters?");

    x.lines().map(|x| x.to_string()).collect::<Vec<String>>()
}
